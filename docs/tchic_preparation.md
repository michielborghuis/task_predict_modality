# T-ChIC preparation and benchmark

This route prepares H3K4me3 counts and unspliced RNA for the OpenProblems
modality prediction benchmark. It currently supports H3K4me3 only.

## Data preparation

The preparation workflow:
- Selects cells and genes using the filtered dataset and the existing plate split.
- Reads raw counts and adds cell annotations.
- Normalizes each modality with the official `log_cp` component: CP10k followed by log1p.
- Creates four benchmark files: `train_mod1.h5ad`, `train_mod2.h5ad`,
  `test_mod1.h5ad` and `test_mod2.h5ad`.

Modality 1 is H3K4me3; modality 2 is RNA. The `test` files contain validation
cells. Final-test cells are excluded. This workflow does not generate a split.

Annotations are matched by sample, mark and barcode. The fixed prefix
`PZ-MB-TChIC-` is removed from annotation sample IDs before matching.

## Build

Run from the repository root with the `thesis` environment active.
Viash 0.9.7, Nextflow and Apptainer must be available.

```bash
python scripts/build_tchic_prepare.py --target target_tchic_prepare
```

This builds the preparation components and processor container. It uses the
datasets dependency pinned to commit `106edbdeae9775d71bdfb4e2dcca67e187dc843b`.
Build records are saved in `target_tchic_prepare/runtime/`.

## Prepare data

Edit the settings and input paths at the top of `scripts/run_tchic_prepare.sh`.
Choose an unused `PREPARATION_NUMBER`, then run:

```bash
bash scripts/run_tchic_prepare.sh
```

The four files are saved under:
`outputs/tchic/prepared/<DATASET>/preparation_<PREPARATION_NUMBER>/published/`.

Container paths are set in `src/workflows/prepare_tchic/runtime.config`.
These paths refer to the current local installation.

## Run the benchmark

Edit `scripts/run_benchmark/run_tchic_local.sh`. Select the same `DATASET`
and `PREPARATION_NUMBER`, and choose an unused `RUN_NUMBER`.

```bash
bash scripts/run_benchmark/run_tchic_local.sh
```

The current setup runs `mean_per_gene` and the official correlation and error
metrics. It uses the existing `target_benchmark_complete` build and
`scripts/run_benchmark/mean_per_gene.config`. That configuration also points
the correlation metric to Python in the local `thesis` environment.
Other methods may need different runtime settings.

Results are saved under:
`outputs/tchic/runs/<DATASET>/<METHOD>_<RUN_NUMBER>/publish/`.

The five outputs are `score_uns.yaml`, `dataset_uns.yaml`,
`method_configs.yaml`, `metric_configs.yaml` and `task_info.yaml`.
Nextflow also publishes a workflow state file.

Both scripts protect existing output directories and save run details and logs.

## Checked execution

- 2026-11-05 - `preparation_001`: all 6 processes completed successfully.
- 2026-11-05 - `mean_per_gene_004`: all 9 processes completed successfully.

This shows success for the local route with this dataset and baseline. Earlier runs and checks are recorded under
`outputs/tchic/runs/real_20261023_124600/verification/`.
