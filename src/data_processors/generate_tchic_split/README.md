# T-ChIC plate split generator

Splits whole plates into train, validation and final test using cell metadata.
Development benchmarks use validation as test data.

## Generate a split

Run from the repository root with the Python environment active.
Set FILTERED to the input file.

```bash
FILTERED="/path/to/allsamples_genelevel_RNAandChIC_counts.filt-norm.h5ad"

viash run src/data_processors/generate_tchic_split/config.vsh.yaml \
  --engine native -- \
  --input_metadata "$FILTERED" \
  --input_settings src/data_processors/generate_tchic_split/h3k4me3_split_settings.json \
  --output_split outputs/tchic/splits/h3k4me3_plate_split.csv \
  --output_split_info outputs/tchic/splits/h3k4me3_plate_split_info.json
```

Existing outputs are not overwritten. Use the same split CSV for later
preparations to keep the train/validation/test assignment unchanged.

## Compare splits (optional)

Set the two CSV paths. Row order and extra columns are ignored.

```bash
REFERENCE="outputs/tchic/splits/h3k4me3_plate_split.csv"
NEW_SPLIT="outputs/tchic/splits/h3k4me3_plate_split_new.csv"

python src/data_processors/generate_tchic_split/compare_splits.py \
  "$REFERENCE" "$NEW_SPLIT"
```

## Tests

Run the small tests:

```bash
python src/data_processors/generate_tchic_split/tests/test_generate_tchic_split.py
```
