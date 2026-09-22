# BABEL test run

21 September 2026

Check whether BABEL still works in the updated OpenProblems repository before using the framework for my T-ChIC data

Dataset used: BMMC multiome test dataset (`openproblems_neurips2021/bmmc_multiome/swap`). Test of the benchmark setup, not a T-ChIC experiment

Full workflow working: training, prediction, correlation, MSE, score extraction and publishing

## Setup

Official OpenProblems code at `e634c9c`, on `thesis/upstream-validation` branch. Branch at `09d422c` after validation. No changes to upstream source code

Software used:

- Viash 0.9.7
- Nextflow 26.04.4
- Apptainer 1.5.2-1.el8

Machine used: `binfgpu9`

Build command:

```bash
viash ns build -t target_benchmark_complete
```

All 58 configurations built successfully

Generated files in `target_benchmark_complete/`, not committed to Git

## Input data

Four files from:

`resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/`

- `train_mod1.h5ad`
- `train_mod2.h5ad`
- `test_mod1.h5ad`
- `test_mod2.h5ad`

`resources_test` in the updated worktree is a symlink to the original worktree

## Configuration

Four config files used:

- `gpu.config` — Apptainer and GPU access
- `cached_containers.config` — paths to available container images
- `numba.config` — writable Numba cache
- `correlation_python.config` — Python environment for the correlation metric

Original files in:

`outputs/babel_openproblems_upstream_20260921/configuration/`

Identical copies in `docs/thesis/babel_validation/config/`, checked with `cmp`

Configs contain paths specific to my account and the university servers

## Run command

Final successful run from the root of the updated repository on GPU9

```bash
export NXF_SYNTAX_PARSER=v1
export NXF_ANSI_LOG=false

RUN=outputs/babel_openproblems_upstream_20260921

nextflow run target_benchmark_complete/nextflow/workflows/run_benchmark \
  -profile singularity \
  -c "$RUN/configuration/gpu.config" \
  -c "$RUN/configuration/correlation_python.config" \
  -c "$RUN/configuration/cached_containers.config" \
  -c "$RUN/configuration/numba.config" \
  -work-dir /linuxhome/tmp/michiel/babel_upstream_20260921_work \
  -resume \
  -with-trace "$RUN/trace_resume_full_conda.txt" \
  --input_train_mod1 resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad \
  --input_train_mod2 resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad \
  --input_test_mod1 resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod1.h5ad \
  --input_test_mod2 resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod2.h5ad \
  --method_ids babel \
  --output_scores score_uns.yaml \
  --output_method_configs method_configs.yaml \
  --output_metric_configs metric_configs.yaml \
  --output_dataset_info dataset_uns.yaml \
  --output_task_info task_info.yaml \
  --publish_dir "$RUN/benchmark" \
  > "$RUN/nextflow_resume_full_conda.log" 2>&1
```

This is the command used for the successful run

The config paths above point to the original experiment files in `outputs/`. Copies also available in this README's `config/` directory

## Execution

Final Nextflow exit code: `0`

- BABEL training — completed earlier on September 21, reused from cache
- BABEL prediction — completed earlier on September 21, reused from cache
- Correlation — executed successfully, exit code 0
- MSE — executed successfully, exit code 0
- Score extraction — both tasks completed, exit code 0
- Results and completion state — published successfully

Training and prediction from the September 21 run, not reused from the older September 8 experiment

No completely fresh run with an empty Nextflow cache tested yet

## Results

Scores from:

`outputs/babel_openproblems_upstream_20260921/benchmark/score_uns.yaml`

| Metric | Value |
|---|---:|
| RMSE | 1.8147 |
| MAE | 1.7124 |
| Mean Pearson per cell | 0.3447 |
| Mean Spearman per cell | 0.2403 |
| Mean Pearson per gene | -0.0184 |
| Mean Spearman per gene | -0.0129 |
| Overall Pearson | 0.2955 |
| Overall Spearman | 0.2153 |

Results only for the OpenProblems test dataset, not T-ChIC

## Output files

Experiment directory:

`outputs/babel_openproblems_upstream_20260921/`

Main files:

- `benchmark/score_uns.yaml` — metric values
- `benchmark/output.run_benchmark.state.yaml` — published output references
- `model/` — saved BABEL model
- `prediction/` — saved RNA predictions
- `nextflow_resume_full_conda.log` — successful run log
- `trace_resume_full_conda.txt` — task execution details

Model and prediction copied from the Nextflow work directory and compared with the originals using `cmp`

Outputs ignored by Git

## Problems and fixes

| Problem | Fix |
|---|---|
| Nextflow parser error | `NXF_SYNTAX_PARSER=v1` |
| Development container tags unavailable | Explicit paths to cached `build_main` images |
| Scanpy / Numba import error | `NUMBA_CACHE_DIR=/tmp/numba_cache` |
| Correlation unable to initialize Python | Mount entire Miniforge installation instead of only the `thesis` environment |

Correlation issue details:

- September 8 correlation task worked with the older setup
- September 21 task failed while reading the solution file
- Newer Nextflow launch used `--no-home`
- Mounting only the Conda environment did not fix it
- Mounting `/home/michiel/miniforge3` allowed the solution file to load
- Same change in `correlation_python.config` followed by a successful `-resume` run

Old `resolved_nextflow.config` in the experiment directory was generated before the final fixes, not the configuration of the successful run

## Still to do

- Record checksums of input datasets and container images
- Record exact identities of the utility images downloaded by Nextflow
- Test from an empty cache if needed
- Make the configuration less dependent on paths specific to my account
- Connect the T-ChIC data and run the first baseline
