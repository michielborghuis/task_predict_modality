from pathlib import Path

## VIASH START
meta = {
    "config": "src/data_processors/generate_tchic_split/config.vsh.yaml",
}
## VIASH END


config_text = Path(meta["config"]).read_text(encoding="utf-8")
required_entries = [
    "--input_metadata",
    "--input_settings",
    "--output_manifest",
    "--output_provenance",
]
missing = [entry for entry in required_entries if entry not in config_text]
if missing:
    raise AssertionError(f"Viash config is missing expected entries: {missing}")
