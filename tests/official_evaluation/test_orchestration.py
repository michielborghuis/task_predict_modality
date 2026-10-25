import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import official_evaluation as ev


class EvaluationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.record = {'dataset_id': 'dataset', 'method_id': 'method',
                       'metric_ids': ['rmse'], 'metric_values': [0.125]}

    def h5ad(self, name, genes=('g1', 'g2'), sparse=False):
        path = self.root / name
        with h5py.File(path, 'w') as f:
            for key, value in [('dataset_id', 'dataset'), ('method_id', 'method'), ('normalization_id', 'log_cp10k')]:
                f.create_dataset('uns/' + key, data=value)
            for axis, values in [('obs', ['c1', 'c2']), ('var', list(genes))]:
                group = f.create_group(axis)
                group.attrs['_index'] = '_index'
                group.create_dataset('_index', data=values, dtype=h5py.string_dtype())
            if sparse:
                f.create_group('layers/normalized').attrs['shape'] = np.array([2, 2])
            else:
                f.create_dataset('layers/normalized', data=np.ones((2, 2)))
        return path

    def test_input_dense_and_sparse_metadata(self):
        pair = ev.inspect_pair(self.h5ad('target.h5ad', sparse=True), self.h5ad('prediction.h5ad'), 'log_cp10k')
        self.assertEqual(pair['shape'], [2, 2])

    def test_reordered_genes_rejected(self):
        with self.assertRaisesRegex(ValueError, 'genes IDs/order'):
            ev.inspect_pair(self.h5ad('target.h5ad'), self.h5ad('prediction.h5ad', genes=('g2', 'g1')), 'log_cp10k')

    def test_duplicate_labels_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate var'):
            ev.inspect_pair(self.h5ad('target.h5ad'), self.h5ad('prediction.h5ad', genes=('g1', 'g1')), 'log_cp10k')

    def test_wrong_scale_rejected(self):
        with self.assertRaisesRegex(ValueError, 'normalization'):
            ev.inspect_pair(self.h5ad('target.h5ad'), self.h5ad('prediction.h5ad'), 'counts')

    def test_exported_values_unchanged(self):
        rows = ev.score_rows([self.record], 'run', 'dataset', 'method', ['rmse'])
        self.assertEqual(rows[0]['value'], self.record['metric_values'][0])

    def test_missing_or_duplicate_scores_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Missing or unexpected'):
            ev.score_rows([self.record], 'run', expected_metrics=['rmse', 'mae'])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            ev.score_rows([self.record, self.record], 'run')

    def test_wrong_dataset_rejected(self):
        with self.assertRaisesRegex(ValueError, 'dataset mismatch'):
            ev.score_rows([self.record], 'run', expected_dataset='different')

    def test_future_component_discovered(self):
        wf = self.root / 'src/workflows/run_benchmark'
        wf.mkdir(parents=True)
        (wf / 'config.vsh.yaml').write_text(yaml.safe_dump({'dependencies': [{'name': 'metrics/future_metric'}]}))
        folder = self.root / 'build/nextflow/metrics/future_metric'
        folder.mkdir(parents=True)
        (folder / 'main.nf').write_text('// fixture')
        (folder / '.config.vsh.yaml').write_text(yaml.safe_dump({
            'info': {'type': 'metric', 'metrics': [{'name': 'future_value'}]},
            'arguments': [{'name': n, 'required': True} for n in ('--input_prediction', '--input_test_mod2', '--output')]}))
        self.assertEqual(ev.discover_metrics(self.root, self.root / 'build')[0]['metric_ids'], ['future_value'])

    def plan(self):
        output = self.root / 'evaluation'
        stage = output / 'run/mse/export'
        stage.mkdir(parents=True)
        (output / 'provenance.json').write_text(json.dumps({'inputs': [], 'files': []}))
        return {'output': str(output), 'nextflow_version': '26.04.4',
                'steps': [{'run': 'run', 'metric': 'mse', 'kind': 'export', 'cwd': str(stage),
                           'command': ['mocked-component'], 'output': str(stage / 'publish/metadata.yaml')}],
                'jobs': [{'run': 'run', 'metadata': {'dataset_id': 'dataset', 'method_id': 'method'}}],
                'metrics': [{'name': 'mse', 'metric_ids': ['rmse']}], 'existing': []}

    def test_missing_output_fails_even_with_zero_exit(self):
        plan = self.plan()
        with patch.object(ev.subprocess, 'run'):
            with self.assertRaisesRegex(ValueError, 'Missing output'):
                ev.execute(plan)
        self.assertEqual(json.loads((Path(plan['output']) / 'status.json').read_text())['state'], 'failed')
        self.assertFalse((Path(plan['output']) / 'scores_wide.csv').exists())

    def test_export_collection_and_single_use(self):
        plan = self.plan()
        def simulated_export(*args, **kwargs):
            destination = Path(plan['steps'][0]['output'])
            destination.parent.mkdir()
            destination.write_text(yaml.safe_dump({'uns': self.record}))
        with patch.object(ev.subprocess, 'run', side_effect=simulated_export):
            frame = ev.execute(plan)
        self.assertEqual(frame['value'].iloc[0], 0.125)
        saved = ev.read_yaml(Path(plan['output']) / 'run/score_uns.yaml')
        self.assertEqual(saved, [self.record])
        with self.assertRaisesRegex(ValueError, 'already been attempted'):
            ev.execute(plan)


if __name__ == '__main__':
    unittest.main()
