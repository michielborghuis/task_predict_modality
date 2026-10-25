# T-ChIC plate split generator

This Viash component reads only `obs` metadata from a filtered H5AD (opened in
backed/read-only mode), aggregates one row per `sample` plate, and searches
candidate plate assignments. It does not read expression/count matrices.

The settings in `h3k4me3_split_settings.json` reproduce the current H3K4me3
experiment. `test` in the generated manifest means the **final holdout**;
downstream development preparation must use `validation` for the benchmark
`test_mod1` and `test_mod2` slots. Daily preparation must reuse a reviewed
manifest and must not run this search again.

Run explicitly when creating or reviewing a split (example command; paths are
placeholders and this command has not been run):

```bash
viash run src/data_processors/generate_tchic_split/config.vsh.yaml \
  --engine native -- \
  --input_metadata /path/to/allsamples_genelevel_RNAandChIC_counts.filt-norm.h5ad \
  --input_settings src/data_processors/generate_tchic_split/h3k4me3_split_settings.json \
  --output_manifest outputs/tchic/splits/generated_h3k4me3_plate_split.csv \
  --output_provenance outputs/tchic/splits/generated_h3k4me3_plate_split.provenance.json
```

The component refuses to overwrite either output. Keep the existing reviewed
CSV, inspect the new manifest and provenance, then compare assignments without
depending on CSV order or formatting:

```bash
python3 src/data_processors/generate_tchic_split/compare_splits.py \
  outputs/tchic/splits/h3k4me3_plate_split.csv \
  outputs/tchic/splits/generated_h3k4me3_plate_split.csv
```

The current settings retain the notebook procedure: `sample` is the plate key;
metadata is grouped in stable sorted plate order; day and plate constraints
are explicit; candidates use `numpy.random.default_rng(seed=42)`; and equal
scores retain the first candidate (`<`, not `<=`). The provenance records
versions, metadata fingerprint, settings, score terms, number of valid
candidates, RNG bit generator, and the manifest hash.

The stored notebook environment reported Python at
`/home/michiel/miniforge3/envs/thesis/bin/python`, AnnData 0.11.4, NumPy 2.2.6,
and pandas 2.3.3. The original NumPy/pandas versions are therefore known from
the notebook's saved output; the actual run provenance is still emitted for
each new generation.

Tests are in `tests/test_generate_tchic_split.py`. Run the focused Python tests
with the task's Python environment, and run the Viash component checks with
`viash test src/data_processors/generate_tchic_split/config.vsh.yaml`. Do not
use the full-data split search merely as a routine component check.
