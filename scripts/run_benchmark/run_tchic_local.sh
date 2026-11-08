#!/usr/bin/env bash

set -euo pipefail # bash settings

# Settings: run this script from within the repository.
REPO=$(git rev-parse --show-toplevel)
SCRIPT=$(realpath "${BASH_SOURCE[0]}")

METHOD="mean_per_gene"
RUN_NUMBER="004"
DATASET="h3k4me3_genelevel_dev"

PREPARATION_NUMBER="001"

PREPARED="$REPO/outputs/tchic/prepared/$DATASET/preparation_$PREPARATION_NUMBER/published" # 4 input files
RUNTIME="$REPO/scripts/run_benchmark/mean_per_gene.config" # runtime configs
WORKFLOW="$REPO/target_benchmark_complete/nextflow/workflows/run_benchmark" # workflow to run

# Derived output directory
RUN="$REPO/outputs/tchic/runs/$DATASET/${METHOD}_${RUN_NUMBER}"

# Check required files.
for FILE in \
    "$PREPARED/train_mod1.h5ad" \
    "$PREPARED/train_mod2.h5ad" \
    "$PREPARED/test_mod1.h5ad" \
    "$PREPARED/test_mod2.h5ad" \
    "$RUNTIME" \
    "$WORKFLOW/main.nf"; do
    test -f "$FILE" || {
        echo "Missing file: $FILE" >&2
        exit 1
    }
done

# Create a new run directory; refuse to overwrite an existing run.
mkdir -p "$(dirname "$RUN")"
mkdir "$RUN"

# Preserve the script, configuration and version information.
cp "$SCRIPT" "$RUN/run_command.sh"
cp "$RUNTIME" "$RUN/runtime.config"
git -C "$REPO" rev-parse HEAD > "$RUN/source_commit.txt"
cp "$WORKFLOW/.config.vsh.yaml" "$RUN/workflow_build_config.yaml"

cd "$RUN"
export NXF_SYNTAX_PARSER=v1 # nextflow syntax parser
export NXF_ANSI_LOG=false # plain text terminal output

nextflow run "$WORKFLOW" \
  -profile singularity \
  -c "$RUN/runtime.config" \
  -work-dir "$RUN/work" \
  -with-trace "$RUN/trace.txt" \
  --method_ids "$METHOD" \
  --input_train_mod1 "$PREPARED/train_mod1.h5ad" \
  --input_train_mod2 "$PREPARED/train_mod2.h5ad" \
  --input_test_mod1 "$PREPARED/test_mod1.h5ad" \
  --input_test_mod2 "$PREPARED/test_mod2.h5ad" \
  --output_scores score_uns.yaml \
  --output_method_configs method_configs.yaml \
  --output_metric_configs metric_configs.yaml \
  --output_dataset_info dataset_uns.yaml \
  --output_task_info task_info.yaml \
  --publish_dir "$RUN/publish" \
  2>&1 | tee "$RUN/run.log"
