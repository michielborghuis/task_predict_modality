"""Orchestrate existing OpenProblems components; never calculate metric values."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import yaml


def read_yaml(path):
    return yaml.safe_load(Path(path).read_text())


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def resolve(repo, path):
    return (Path(repo) / Path(path).expanduser()).resolve()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def discover_metrics(repo, build):
    """Use workflow dependencies as the metric registry; do not hard-code names."""
    workflow = read_yaml(repo / 'src/workflows/run_benchmark/config.vsh.yaml')
    result = []
    for dependency in workflow['dependencies']:
        name = dependency['name'] if isinstance(dependency, dict) else dependency
        if not name.startswith('metrics/'):
            continue
        folder = build / 'nextflow' / name
        config = read_yaml(folder / '.config.vsh.yaml')
        require(config['info']['type'] == 'metric', f'Not a metric: {name}')
        arguments = config.get('arguments', []) + [
            a for group in config.get('argument_groups', []) for a in group.get('arguments', [])
        ]
        required = {a['name'] for a in arguments if a.get('required')}
        require(required <= {'--input_prediction', '--input_test_mod2', '--output'},
                f'{name} requires extra arguments: {required}. Add an explicit adapter first.')
        ids = [m['name'] for m in config['info']['metrics']]
        require(ids and len(ids) == len(set(ids)), f'Invalid metric IDs: {name}')
        require((folder / 'main.nf').is_file(), f'Missing built component: {folder}')
        result.append({'name': name.split('/')[-1], 'folder': str(folder), 'metric_ids': ids})
    require(result, 'No metric dependencies found in run_benchmark.')
    require(len({m['name'] for m in result}) == len(result), 'Duplicate metric components.')
    all_ids = [i for m in result for i in m['metric_ids']]
    require(len(set(all_ids)) == len(all_ids), 'Duplicate metric IDs across components.')
    return result


def inspect_pair(target, prediction, normalization):
    """Read HDF5 metadata and axis labels without loading expression matrices."""
    import h5py

    def text(value):
        return value.decode() if isinstance(value, bytes) else str(value)

    def inspect(path):
        with h5py.File(path, 'r') as f:
            uns = f['uns']
            info = {key: text(uns[key][()]) for key in ('dataset_id', 'method_id', 'normalization_id') if key in uns}
            layer = f['layers/normalized']
            info['shape'] = list(layer.shape if isinstance(layer, h5py.Dataset) else layer.attrs['shape'])
            labels = []
            for axis in ('obs', 'var'):
                group = f[axis]
                key = text(group.attrs.get('_index', '_index'))
                values = [text(v) for v in group[key][()]]
                require(len(values) == len(set(values)), f'{path}: duplicate {axis} names')
                labels.append(values)
            require(info['shape'] == [len(v) for v in labels], f'{path}: inconsistent axes')
        return info, labels

    truth, tlabels = inspect(target)
    pred, plabels = inspect(prediction)
    require(truth.get('dataset_id') and truth.get('dataset_id') == pred.get('dataset_id'), 'Dataset IDs differ or are missing.')
    require(pred.get('method_id'), 'Prediction has no method_id.')
    require(truth['shape'] == pred['shape'], 'Prediction and target shapes differ.')
    for axis, a, b in zip(('cells', 'genes'), tlabels, plabels):
        require(a == b, f'Prediction and target {axis} IDs/order differ. Stop and establish alignment; do not relabel blindly.')
    require(truth.get('normalization_id') == normalization, 'Target normalization differs from the configured scale.')
    require(pred.get('normalization_id', normalization) == normalization, 'Prediction normalization metadata conflicts with target.')
    return {'dataset_id': truth['dataset_id'], 'method_id': pred['method_id'],
            'shape': [int(x) for x in truth['shape']], 'normalization': normalization,
            'prediction_scale_metadata_present': 'normalization_id' in pred}


def score_rows(records, run_id, expected_dataset=None, expected_method=None, expected_metrics=None):
    require(isinstance(records, list) and records, 'Empty or invalid score list.')
    rows, seen = [], set()
    for record in records:
        dataset, method = record['dataset_id'], record['method_id']
        require(expected_dataset is None or dataset == expected_dataset, 'Score dataset mismatch.')
        require(expected_method is None or method == expected_method, 'Score method mismatch.')
        ids, values = record['metric_ids'], record['metric_values']
        require(len(ids) == len(values) and ids, 'Missing or unequal metric IDs and values.')
        for metric, value in zip(ids, values):
            key = (dataset, method, metric)
            require(key not in seen, f'Duplicate score: {key}')
            seen.add(key)
            require(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value), f'Invalid score: {metric}')
            rows.append({'run': run_id, 'dataset_id': dataset, 'method_id': method, 'metric': metric, 'value': value})
    if expected_metrics is not None:
        require({r['metric'] for r in rows} == set(expected_metrics), 'Missing or unexpected metric outputs.')
    return rows


def runtime_text(python, conda_root, overrides):
    # Same image overrides and reticulate binding as the successful baseline notebook.
    def quote(value):
        require(not any(c in str(value) for c in "'\n\r"), 'Unsupported quote/newline in runtime path.')
        return "'" + str(value) + "'"
    lines = ["process {", "  withName: '.*' { errorStrategy = 'terminate' }"]
    for name, image in overrides.items():
        require(re.fullmatch(r'[A-Za-z0-9_]+', name), f'Invalid component name: {name}')
        lines += [f"  withName: '.*{name}_process' {{", f'    container = {quote(image)}']
        if name == 'correlation':
            require(not any(c.isspace() for c in str(python) + str(conda_root)), 'Reticulate paths cannot contain whitespace.')
            options = f'--env RETICULATE_PYTHON={python} --bind {conda_root}:{conda_root}'
            lines += [f'    containerOptions = {quote(options)}']
        lines += ['  }']
    return '\n'.join(lines + ['}', ''])


def command(component, config, stage, arguments):
    cmd = ['nextflow', '-log', str(stage / 'nextflow.log'), 'run', str(component),
           '-profile', 'singularity', '-c', str(config), '-work-dir', str(stage / 'work'),
           '-with-trace', str(stage / 'trace.tsv'), '--publish_dir', str(stage / 'publish')]
    for key, value in arguments.items():
        cmd += ['--' + key, str(value)]
    return cmd


def prepare(repo, settings):
    """Validate inputs, freeze commands and provenance, but run no evaluation."""
    repo = Path(repo).resolve()
    settings = json.loads(json.dumps(settings))
    build = resolve(repo, settings['build_dir'])
    output = resolve(repo, settings['output_dir'])
    require(not output.exists(), f'Output already exists: {output}. Choose a new run directory; existing scores are never overwritten.')
    require(shutil.which('nextflow'), 'nextflow is not on PATH. Use the working thesis environment.')
    require(shutil.which('singularity'), 'singularity is not on PATH. Use the same Apptainer/Singularity setup as notebook 03.')
    prefix = os.environ.get('CONDA_PREFIX')
    require(prefix, 'Start Jupyter from the activated thesis Conda environment.')
    prefix = Path(prefix).resolve()
    python = Path(sys.executable).resolve()
    require(python.is_relative_to(prefix), 'Notebook kernel does not belong to the activated Conda environment.')
    conda_root = prefix.parent.parent if prefix.parent.name == 'envs' else prefix
    metrics = discover_metrics(repo, build)
    extractors = sorted(build.glob('dependencies/**/nextflow/utils/extract_uns_metadata/main.nf'))
    require(len(extractors) == 1, f'Expected one metadata exporter, found {len(extractors)}.')
    exporter = extractors[0].parent
    jobs, labels = [], set()
    for job in settings['predictions']:
        label = job['run']
        require(re.fullmatch(r'[A-Za-z0-9_-]+', label) and label not in labels, 'Run labels must be unique safe folder names.')
        labels.add(label)
        target, pred = resolve(repo, job['target']), resolve(repo, job['prediction'])
        require(job.get('scale_confirmed') is True, f'{label}: confirm prediction scale from the prediction code/run settings first.')
        pair = inspect_pair(target, pred, job['normalization'])
        jobs.append(dict(job, target=str(target), prediction=str(pred), metadata=pair))
    require(jobs, 'No prediction jobs configured.')
    require(len({j['metadata']['dataset_id'] for j in jobs}) == 1,
            'Use one evaluation plan per dataset; analyse different datasets separately.')
    existing = []
    for item in settings.get('existing_scores', []):
        require(item['run'] not in labels, 'Duplicate run label.')
        labels.add(item['run'])
        path = resolve(repo, item['path'])
        records = read_yaml(path)
        require(item['dataset_id'] == jobs[0]['metadata']['dataset_id'], 'Existing scores belong to a different dataset.')
        rows = score_rows(records, item['run'], item['dataset_id'], item['method_id'],
                          [i for m in metrics for i in m['metric_ids']])
        existing.append(dict(item, path=str(path), records=records, rows=rows))
    def capture(cmd):
        return subprocess.run(cmd, cwd=repo, check=True, text=True, capture_output=True).stdout.strip()
    git_commit = capture(['git', 'rev-parse', 'HEAD'])
    upstream = capture(['git', 'rev-parse', 'upstream/main'])
    folders = ['src/metrics/' + m['name'] for m in metrics]
    difference = capture(['git', 'diff', 'upstream/main', '--'] + folders)
    require(not difference, 'Metric source differs from locally stored upstream/main. Review before evaluation.')
    environment = os.environ.copy()
    environment.update(NXF_SYNTAX_PARSER='v1', NXF_VER=settings['nextflow_version'])
    version = subprocess.run(['nextflow', '-version'], cwd=repo, env=environment, check=True,
                             text=True, capture_output=True).stdout.strip()
    tracked = [repo / 'src/workflows/run_benchmark/config.vsh.yaml', Path(__file__).resolve()]
    for folder in [Path(m['folder']) for m in metrics] + [exporter]:
        tracked += [folder / n for n in ('main.nf', '.config.vsh.yaml', 'nextflow.config')]
    for name in folders:
        tracked += [p for p in (repo / name).rglob('*') if p.is_file()]
    inputs = [fingerprint(p) for j in jobs for p in (j['target'], j['prediction'])]
    inputs += [fingerprint(item['path']) for item in existing]
    provenance = {'created_utc': datetime.now(timezone.utc).isoformat(), 'git_commit': git_commit,
                  'upstream_commit': upstream, 'git_status': capture(['git', 'status', '--short']),
                  'nextflow_version_output': version, 'python': str(python), 'conda_prefix': str(prefix),
                  'python_packages': {name: importlib.metadata.version(name) for name in ('h5py', 'PyYAML', 'pandas')},
                  'files': [fingerprint(p) for p in sorted(set(tracked))], 'inputs': inputs,
                  'environment': {k: os.environ.get(k) for k in ('NXF_SINGULARITY_CACHEDIR','APPTAINER_CACHEDIR','SINGULARITY_CACHEDIR')},
                  'limitations': ['Container tags are mutable; retain container caches for exact replay.',
                                 'Prediction scale confirmation is supplied by the user; labels alone cannot prove scale.',
                                 'Existing score files retain their original run provenance; they are not re-evaluated here.']}
    config_text = runtime_text(python, conda_root, settings['container_overrides'])
    output.mkdir(parents=True)
    config = output / 'runtime.config'
    config.write_text(config_text)
    provenance['files'].append(fingerprint(config))
    (output / 'settings.json').write_text(json.dumps(settings, indent=2))
    (output / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    steps = []
    for job in jobs:
        for metric in metrics:
            base = output / job['run'] / metric['name']
            score = base / 'metric/publish/score.h5ad'
            cutoff = max(10, len(metric['metric_ids']))
            for kind, component, args, result in (
                ('metric', Path(metric['folder']), {'input_prediction': job['prediction'], 'input_test_mod2': job['target'], 'output': 'score.h5ad'}, score),
                ('export', exporter, {'input': str(score), 'output': 'metadata.yaml', 'uns_length_cutoff': cutoff}, base / 'export/publish/metadata.yaml'),
            ):
                stage = base / kind
                stage.mkdir(parents=True)
                cmd = command(component, config, stage, args)
                steps.append({'run': job['run'], 'metric': metric['name'], 'kind': kind,
                              'cwd': str(stage), 'command': cmd, 'output': str(result)})
    plan = {'repo': str(repo), 'output': str(output), 'metrics': metrics, 'jobs': jobs,
            'existing': existing, 'steps': steps, 'nextflow_version': settings['nextflow_version']}
    (output / 'plan.json').write_text(json.dumps(plan, indent=2))
    script = '#!/usr/bin/env bash\nset -euo pipefail\n' + shlex.join(
        [str(python), str(Path(__file__).resolve()), '--execute-plan', str(output / 'plan.json')]) + '\n'
    (output / 'run_command.sh').write_text(script)
    (output / 'run_command.sh').chmod(0o755)
    return plan


def execute(plan):
    """Run frozen component commands and collect exported values unchanged."""
    output = Path(plan['output'])
    require(not (output / 'status.json').exists(), 'This plan has already been attempted. Use a new output directory for a fresh attempt.')
    provenance = json.loads((output / 'provenance.json').read_text())
    for expected in provenance['inputs'] + provenance['files']:
        require(fingerprint(expected['path']) == expected, f'Input/code changed after planning: {expected["path"]}')
    env = os.environ.copy()
    env.update(NXF_SYNTAX_PARSER='v1', NXF_VER=plan['nextflow_version'])
    def status(state, **extra):
        (output / 'status.json').write_text(json.dumps({'state': state, **extra}, indent=2))
    status('running')
    try:
        for step in plan['steps']:
            print(f"{step['run']} / {step['metric']} / {step['kind']}", flush=True)
            stage = Path(step['cwd'])
            with (stage / 'console.log').open('w') as log:
                subprocess.run(step['command'], cwd=stage, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            require(Path(step['output']).is_file(), f'Missing output: {step["output"]}; inspect {stage}/console.log')
        all_rows = []
        for job in plan['jobs']:
            records = []
            for metric in plan['metrics']:
                path = output / job['run'] / metric['name'] / 'export/publish/metadata.yaml'
                record = read_yaml(path)['uns']
                pair = job['metadata']
                rows = score_rows([record], job['run'], pair['dataset_id'], pair['method_id'], metric['metric_ids'])
                records.append(record)
                all_rows.extend(dict(row, score_source=str(path), evaluation='new') for row in rows)
            (output / job['run'] / 'score_uns.yaml').write_text(yaml.safe_dump(records, sort_keys=False))
        for item in plan['existing']:
            all_rows.extend(dict(row, score_source=item['path'], evaluation='existing') for row in item['rows'])
        import pandas as pd
        frame = pd.DataFrame(all_rows)
        require(frame['dataset_id'].nunique() == 1, 'Comparison contains different datasets; analyse separately.')
        frame.to_csv(output / 'scores_long.csv', index=False)
        frame.pivot(index=['run', 'method_id'], columns='metric', values='value').to_csv(output / 'scores_wide.csv')
        status('complete')
        return frame
    except BaseException as error:
        status('failed', error=str(error))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute-plan', required=True, type=Path)
    arguments = parser.parse_args()
    execute(json.loads(arguments.execute_plan.read_text()))
