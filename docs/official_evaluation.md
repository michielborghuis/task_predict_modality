# Evaluate saved predictions with OpenProblems

## Scope

This addition contains a small notebook and a shared Python orchestration module.
It changes no files under `src/metrics`, no generated components and no model.
Metric values come only from the existing OpenProblems components. The existing
`utils/extract_uns_metadata` component exports their `.uns` data. Python collects
these records into the same list-of-records YAML structure used by run_benchmark;
it never recomputes the scores.

The initial notebook evaluates BABEL 20 epochs and reads the existing official
mean_per_gene scores. It does **not** claim to have re-evaluated that baseline.
Add any verified prediction/target pair, including a baseline, to `predictions`
to evaluate it through the same code. Full new model runs should continue to use
the existing `run_benchmark` workflow, which already includes evaluation.

## Install without overwriting existing work

Upload `tchic_official_evaluation.zip` to `/home/michiel/` through MobaXterm.
From the benchmark repository root:

```bash
unzip -n /home/michiel/tchic_official_evaluation.zip
```

`-n` preserves existing files. If these new files already exist from an earlier
installation, inspect their differences before replacing them.

Open `notebooks/05_official_evaluation.ipynb` using the same thesis environment
and Jupyter setup as notebook 03. Run only its first three code cells initially:
imports, settings, and input checks. Review their output before proceeding.

## Inputs and scale

The default target is the development dataset's `test_mod2.h5ad`; those slots
contain validation cells, not the final holdout. Exact axis labels/order and
unique labels are checked through HDF5 without loading expression matrices.
If a prediction has generic labels, stop and establish its actual order from
the prediction implementation. Do not silently replace or reorder labels.

The prediction scale must be verified from the run/model prediction code.
`layers['normalized']` and matching dimensions alone do not establish the scale.
Set `scale_confirmed=True` only after this review. No rescaling is performed.

## Components and environment

The metric set is read from `metrics/` dependencies in
`src/workflows/run_benchmark/config.vsh.yaml`. It is not hard-coded in the
orchestration module. A new metric must be registered in the upstream-style
workflow (both dependencies and its metrics list), built with Viash, and expose
the standard prediction/target/output interface and declared metric IDs.
Metrics requiring additional inputs need an explicit adapter; they fail clearly.

The runtime settings reuse the images and reticulate/Conda binding documented in
the successful notebook 03. Nextflow 26.04.4 and the v1 parser are requested to
match that observed run. No tool upgrade is performed by this package, although
the Nextflow launcher may need to obtain that requested version if absent.
The existing `singularity` command/compatibility setup must be on PATH.

The standalone generated components have their own entry workflow, so a new
Nextflow workflow is unnecessary. They are called sequentially, followed by the
official metadata exporter. `errorStrategy='terminate'` is set at runtime and
the helper independently checks each expected output and all metric IDs.

## Outputs and reproducibility

Planning creates a fresh output directory containing `settings.json`,
`runtime.config`, `plan.json`, `provenance.json` and `run_command.sh`.
Provenance includes Git HEAD/upstream/status, metric-source and generated-code
SHA256 hashes, prediction/target and existing-score hashes, relevant Python
package versions, Nextflow version output and selected cache environment paths.
Inputs and code are checked again against those hashes before execution.

Each component has separate work/publish folders, Nextflow log, trace and
`console.log`. Each newly evaluated prediction gets an exported `score_uns.yaml`.
`scores_long.csv` preserves run labels and score-source paths; `scores_wide.csv`
is a comparison view. Plot output is `development_scores.png`.

The original baseline YAML is never modified. Its original provenance remains
with its original run: this package does not certify its historical container
contents or whether the target has changed since that earlier run.

Plans are single-use. On failure, inspect the saved logs and use a new output
directory after fixing the cause. This intentionally avoids ambiguous reuse of
partial/stale outputs. Reproducibility does not depend on terminal history:
the saved `run_command.sh` executes the frozen plan outside Jupyter.

**Remaining reproducibility limit:** reused container tags are mutable and this
package does not lock container bytes/digests or the full Conda environment.
Retain the image caches and environment for faithful replay; complete locking
belongs in the shared benchmark runtime setup. Recording tag names is not an
immutable environment specification.

## Validation status

The source files, generated component interfaces, successful notebook 03 and
original baseline YAML were inspected. Local tests cover input alignment,
score completeness/identity, registry discovery and orchestration failure
handling using synthetic data and mocked component execution.

Actual Nextflow/container execution on the T-ChIC inputs is **not yet tested**.
No new BABEL scores have been generated locally. The first server run is the
remaining integration test, after the input/scale review.

Run the local helper tests from the repository root:

```bash
python -m unittest discover -s tests/official_evaluation -v
```

After a successful server run, commit the notebook, helper, tests and this note
in the benchmark repository. Keep result files in their normal outputs location.
