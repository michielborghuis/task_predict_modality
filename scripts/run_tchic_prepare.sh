#!/usr/bin/env bash
set -euo pipefail  # Bash settings

# Settings
REPO=$(git rev-parse --show-toplevel)
SCRIPT=$(realpath "${BASH_SOURCE[0]}")

DATASET="h3k4me3_genelevel_dev"
PREPARATION_NUMBER="001"
MARK="H3K4me3"
MAX_TEST_CELLS="-1"

# Input files
SOURCE="/home/vivek/NOBINFBACKUP/projects/mGast_tchic/adata_files"
RAW="$SOURCE/allsamples_genelevel_RNAandChIC_counts.h5ad"
FILTERED="$SOURCE/allsamples_genelevel_RNAandChIC_counts.filt-norm.h5ad"
ANNOTATIONS="/home/michiel/thesis/data/gastruloid_tchic_counts/gast_H3K4me3_2000_bins_ctype.h5ad"
SPLIT="$REPO/outputs/tchic/splits/h3k4me3_plate_split.csv"

# Build and runtime
TARGET="$REPO/target_tchic_prepare"
RUNTIME="$REPO/src/workflows/prepare_tchic/runtime.config"  # Runtime configs
WORKFLOW="$TARGET/nextflow/workflows/prepare_tchic"        # Workflow to run
BUILD_INFO="$TARGET/runtime/runtime.json"

# Derived output directory
RUN="$REPO/outputs/tchic/prepared/$DATASET/preparation_$PREPARATION_NUMBER"
PREPARED="$RUN/published"  # 4 output files

# Check required files.
for FILE in \
    "$RAW" "$FILTERED" "$ANNOTATIONS" "$SPLIT" \
    "$RUNTIME" "$WORKFLOW/main.nf" "$BUILD_INFO"; do
    test -f "$FILE" || {
        echo "Missing file: $FILE" >&2
        exit 1
    }
done

# Create a new preparation directory; refuse to overwrite an existing one.
mkdir -p "$(dirname "$RUN")"
mkdir "$RUN"

# Preserve the script, split, configuration and version information.
cp "$SCRIPT" "$RUN/run_command.sh"
cp "$SPLIT" "$RUN/plate_split.csv"
cp "$RUNTIME" "$RUN/runtime.config"
cp "$BUILD_INFO" "$RUN/build_runtime.json"
cp "$WORKFLOW/.config.vsh.yaml" "$RUN/workflow_build_config.yaml"
git -C "$REPO" rev-parse HEAD > "$RUN/source_commit.txt"


export NXF_SYNTAX_PARSER=v1  # Nextflow syntax parser
export NXF_ANSI_LOG=false   # Plain-text terminal output

cd "$RUN"

nextflow run "$WORKFLOW" \
  -profile singularity \
  -c "$RUN/runtime.config" \
  -work-dir "$RUN/work" \
  -with-trace "$RUN/trace.txt" \
  --input_raw "$RAW" \
  --input_filtered "$FILTERED" \
  --input_annotations "$ANNOTATIONS" \
  --input_split "$RUN/plate_split.csv" \
  --mark "$MARK" \
  --max_test_cells "$MAX_TEST_CELLS" \
  --output_train_mod1 train_mod1.h5ad \
  --output_train_mod2 train_mod2.h5ad \
  --output_test_mod1 test_mod1.h5ad \
  --output_test_mod2 test_mod2.h5ad \
  --publish_dir "$PREPARED" \
  2>&1 | tee "$RUN/run.log"
