# T-ChIC H3K4me3 to GEX preparation

## Verified status

The synthetic preparation passed on 2026-10-23. The full, non-stub Nextflow
workflow completed all six processes with zero failures and published four
H5AD files. Their content was checked against the synthetic input fixture.
No real T-ChIC count matrices were processed.

The workflow uses the fixed plate split. It sends retained validation cells to
the task `test_mod1` and `test_mod2` files and excludes the final `test`
partition. It does not search for or regenerate a split.

## Component responsibilities

- `prepare_tchic_input` reads filtered metadata for cell and feature selection,
  joins annotations on `(sample, mark, barcode)`, rejects missing or duplicate
  join keys, and extracts aligned `chic_counts` and `unspliced` count layers
  from the raw file. The final test partition is excluded.
- The official OpenProblems `normalization/log_cp` component is run separately
  for H3K4me3 and GEX, with `n_cp=10000`, output layer `normalized`, and
  normalization ID `log_cp10k`. Each modality gets its own `size_factors`.
- The existing `process_dataset` component creates the four task files,
  preserves sparse `counts` and normalized layers and the requested `obs` and
  `var` metadata, and applies its default cap of 1,000 test cells. The
  development workflow's test slots contain validation cells, not final-test
  cells. HVG columns are preserved only when supplied; this workflow does not
  compute or insert them.

The workflow source is
`src/workflows/prepare_tchic/`; the generated workflow used for verification
is `target_tchic_verified_20261023/nextflow/workflows/prepare_tchic`.

## Pinned dependency and build

The datasets source is pinned to
`openproblems-bio/datasets@106edbdeae9775d71bdfb4e2dcca67e187dc843b`.
Verify the checkout from the repository root:

```bash
git -C .dependencies/datasets-106edbde rev-parse HEAD
```

The verified runner target is
`.dependencies/datasets-106edbde/target_viash097/`. The build helper checks
and reuses it. Run the bootstrap script only if this target does not exist:

```bash
bash scripts/bootstrap_tchic_datasets.sh \
  .dependencies/datasets-106edbde/target_viash097
```

During the component build,
Viash's local repository binding points to the pinned checkout root,
`.dependencies/datasets-106edbde/`, so it can resolve the component's build
metadata; the generated runner remains in `target_viash097/`.

Build the preparation components and isolated runtime with Viash 0.9.7 and
Apptainer:

```bash
/home/michiel/miniforge3/envs/thesis/bin/python scripts/build_tchic_prepare.py \
  --target target_tchic_verified_20261023
```

The helper checks the source pin and base image checksums, reuses the matching
verified runtime/dependency target, builds the generated wrappers and records
source and artifact hashes. The leading build is
`target_tchic_verified_20261023/`; its runtime identity is
`target_tchic_verified_20261023/runtime/runtime.json`. The isolated processor
image is
`target_tchic_verified_20261023/runtime/process_dataset.sif`.

Recorded image identities:

| Image | SHA-256 |
|---|---|
| OpenProblems Python base SIF | `c7a010cf9a282565c9401e74de3438ab33587d770bc52112b4981b6091fdd6fb` |
| OpenProblems R base SIF | `d9406ebde369d5b6e4dff7441ec9ba4866d2c99a8a5d9d35b0699cb6f4fb04da` |
| Isolated `process_dataset` SIF | `53d0e61475f7d501b70ca9335bde286eb2479890afee734d337a17ba9fdd65b6` |

The runtime lock also records source hashes for the extractor
(`4f59dfd59347429a2f036217f47df17f44b330b41393d5a16a7a32ed9de8fd92`),
processor (`277aaa9138b5024a87efb5383f477b0130463ec3ae0a54d2fc10aaa868652f44`),
workflow (`e0335141953bfe7da7a8407d4320113782a9f869d6c37498e77cb185c26569e7`),
and runtime configuration
(`f080a014bed940e2b0e6baaef327a808a756b7c465c88db09ee9662e8bbc8293`).
It also records artifact hashes. The pinned `log_cp` executable is
`65424554ccb4a71c4a450f52f0fd7952d2ff2e969a4a4540e16e981654c703dc`;
its generated Nextflow runner is
`d7bd4a19e1fd818699879ebecf76df8ee5b6bbd7638bda6fbdad60a84bf9b170`.

## Full synthetic verification

This is the command used for the passing, non-stub run. The `--output-dir`
must be a new directory; the recorded directory already exists.

```bash
APPTAINER_CACHEDIR=/home/michiel/.apptainer/cache \
NXF_SINGULARITY_CACHEDIR=/home/michiel/.apptainer/cache/nextflow_tchic \
NXF_SYNTAX_PARSER=v1 \
  /home/michiel/miniforge3/envs/thesis/bin/python \
  tests/run_tchic_synthetic_workflow.py \
  --target target_tchic_verified_20261023 \
  --output-dir outputs/tchic/verification_20261023/full_synthetic_20261023
```

The successful harness checked cell and feature selection/order, annotation
mapping and expected failures for missing/duplicate join keys, exact fixture
counts, exclusion of final-test cells, and independent official `log_cp`
normalization and size factors (including zero rows). It also checked task
metadata and output contracts; the OpenProblems validator reported all four
schemas passing. Processor checks covered the default limit of 1,000 from
1,002 validation cells, `max_test_cells=-1` retaining all 1,002, and HVG
columns both absent and preserved when present. The synthetic fixture has
different source row orders for raw, filtered, and annotation inputs.

Evidence and outputs:

- `outputs/tchic/verification_20261023/full_synthetic_20261023/result.json`
- `outputs/tchic/verification_20261023/full_synthetic_20261023/commands.json`
- `outputs/tchic/verification_20261023/full_synthetic_20261023/logs/`
- `outputs/tchic/verification_20261023/full_synthetic_20261023/nextflow.log`
- `outputs/tchic/verification_20261023/full_synthetic_20261023/published/`

The published files were independently rechecked against the original fixture:

| File | Cells | SHA-256 |
|---|---:|---|
| `train_mod1.h5ad` | 2 | `04452e5b16f6b2eeacbef2ed19a7cde8e213370212d023ccd10b99ef9e24194a` |
| `train_mod2.h5ad` | 2 | `b11e7da59b816af89683fced1a714a126698db5aa72615650b10adb0f579106d` |
| `test_mod1.h5ad` | 1,002 | `8c0d3b3d58539e9e9b32b8e7af1f9d255545f7af63f5ee17d92d43216101756f` |
| `test_mod2.h5ad` | 1,002 | `8ca26bf75b392bb2fb2ffa99751c641198c82081891a94cf6c0a54e9ce372c2a` |

The output layers remain sparse; `counts` match the original fixture values
exactly. Both modalities have their own normalized values and size factors.

## Real-data status and future invocation

Real-data preparation has **not** been run. The paths below are placeholders;
confirm the real input files, split manifest, annotation keys, and output
directory before using them. In particular, verify Python `anndata` is
available to reticulate in the selected processor image and that the real
files have the expected count layers and metadata.

From the repository root, with the image paths from the verified runtime lock:

```bash
export TCHIC_PYTHON_IMAGE=/home/michiel/.apptainer/cache/nextflow_tchic/openproblems-base_python-1.img
export TCHIC_PROCESSOR_IMAGE=/home/michiel/thesis/benchmark/task_predict_modality_upstream_validation/target_tchic_verified_20261023/runtime/process_dataset.sif
export NXF_SYNTAX_PARSER=v1
export APPTAINER_CACHEDIR=/home/michiel/.apptainer/cache
export NXF_SINGULARITY_CACHEDIR=/home/michiel/.apptainer/cache/nextflow_tchic

nextflow run target_tchic_verified_20261023/nextflow/workflows/prepare_tchic \
  -profile singularity \
  -c src/workflows/prepare_tchic/runtime.config \
  --publish_dir outputs/tchic/prepared/h3k4me3_genelevel_dev \
  --input_raw /path/to/raw_genelevel.h5ad \
  --input_filtered /path/to/filtered_genelevel.h5ad \
  --input_annotations /path/to/gast_H3K4me3_2000_bins_ctype.h5ad \
  --input_split outputs/tchic/splits/h3k4me3_plate_split.csv \
  --output_train_mod1 train_mod1.h5ad \
  --output_train_mod2 train_mod2.h5ad \
  --output_test_mod1 test_mod1.h5ad \
  --output_test_mod2 test_mod2.h5ad \
  --mark H3K4me3 \
  --annotation_sample_prefix PZ-MB-TChIC- \
  --max_test_cells -1
```

The workflow direction is fixed: mod1 is H3K4me3 and mod2 is GEX. The declared
argument names are `--output_train_mod1`, `--output_train_mod2`,
`--output_test_mod1`, and `--output_test_mod2`; there is no `--output_dir` or
`--swap` argument. The publish directory must be new. These settings and
examples document the intended invocation only; they are not evidence of a
real-data run.
