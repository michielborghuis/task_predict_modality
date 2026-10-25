"""Data preparation helpers for the Full Run Benchmark Notebook."""

import numpy as np
import scipy.sparse as sp
import pandas as pd
import h5py
from anndata.io import read_elem
from anndata.io import sparse_dataset
import anndata as ad
import hashlib
import json
from pathlib import Path


def log_cp10k(counts):
    """Normalize each cell to 10,000 counts, then apply log1p."""
    normalized = sp.csr_matrix(counts, dtype=np.float64, copy=True)
    totals = np.asarray(normalized.sum(axis=1)).ravel()

    scale = np.divide(
        10_000,
        totals,
        out=np.zeros_like(totals),
        where=totals > 0,
    )

    normalized = normalized.multiply(scale[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized


def load_plate_split(path, histone_mark):
    """Load and validate plate assignments for one histone mark."""
    split = pd.read_csv(path, dtype=str)

    required = {"sample", "mark", "split"}
    if not required.issubset(split.columns):
        raise ValueError(f"Split file must contain columns: {sorted(required)}")

    split = split.loc[split["mark"] == histone_mark].copy()

    if split.empty:
        raise ValueError(f"No plates found for {histone_mark}.")
    if split["sample"].isna().any() or split["sample"].duplicated().any():
        raise ValueError("Each plate must have one unique, non-missing sample ID.")
    if set(split["split"]) != {"train", "validation", "test"}:
        raise ValueError("Split must contain train, validation and test assignments.")

    return split


def select_cells_and_genes(filtered_path, plate_split, histone_mark):
    """Read retained metadata and assign cells to the saved plate split."""
    with h5py.File(filtered_path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])

    obs = obs.loc[obs["mark"].astype(str) == histone_mark].copy()

    if obs.empty or var.empty:
        raise ValueError("No cells or genes found for this selection.")
    if not obs.index.is_unique or not var.index.is_unique:
        raise ValueError("Cell and gene IDs must be unique.")

    assignments = plate_split.set_index("sample")["split"]
    partitions = obs["sample"].astype(str).map(assignments)

    if partitions.isna().any():
        missing = obs.loc[partitions.isna(), "sample"].unique()
        raise ValueError(f"Cells belong to plates missing from the split: {missing}")

    cells = {
        part: obs.index[partitions == part].tolist()
        for part in ("train", "validation", "test")
    }

    if any(len(ids) == 0 for ids in cells.values()):
        raise ValueError("Every partition must contain retained cells.")

    return obs, var, cells


def read_raw_counts(raw_path, layer, cell_ids, gene_ids):
    """Read selected raw counts in the requested cell and gene order."""
    with h5py.File(raw_path, "r") as handle:
        raw_obs = read_elem(handle["obs"])
        raw_var = read_elem(handle["var"])

        if not raw_obs.index.is_unique or not raw_var.index.is_unique:
            raise ValueError("Raw cell and gene IDs must be unique.")

        row_indices = raw_obs.index.get_indexer(cell_ids)
        col_indices = raw_var.index.get_indexer(gene_ids)

        if (row_indices < 0).any() or (col_indices < 0).any():
            raise ValueError("Selected cells or genes are missing from the raw file.")

        if layer not in handle["layers"]:
            raise ValueError(f"Raw layer not found: {layer}")

        matrix = sparse_dataset(handle["layers"][layer])
        counts = matrix[row_indices, :][:, col_indices].tocsr()

    values = counts.data
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Counts must be finite and non-negative.")
    if not np.allclose(values, np.rint(values), rtol=0, atol=1e-6):
        raise ValueError("The selected layer does not contain integer counts.")
    if values.size and values.max() > np.iinfo(np.int32).max:
        raise ValueError("Counts exceed the supported integer range.")

    counts.data = np.rint(values)
    return counts.astype(np.int32)


def make_benchmark_adata(
    raw_path, raw_layer, cell_ids, obs, var, modality,
    dataset_id, normalization="log_cp10k",
):
    """Build one benchmark input with raw and normalized counts."""
    if normalization != "log_cp10k":
        raise ValueError(f"Unsupported normalization: {normalization}")

    counts = read_raw_counts(raw_path, raw_layer, cell_ids, var.index)

    cell_metadata = obs.loc[cell_ids].copy()
    cell_metadata["batch"] = cell_metadata["sample"].astype(str)

    gene_metadata = var.copy()
    gene_metadata["gene_ids"] = gene_metadata.index.astype(str)
    gene_metadata["hvg"] = False
    gene_metadata["hvg_score"] = 0.0

    metadata = {
        "dataset_id": dataset_id,
        "modality": modality,
        "normalization_id": normalization,
        "dataset_name": dataset_id,
        "dataset_summary": f"Modality prediction dataset: {dataset_id}.",
        "dataset_description": (
            f"Dataset {dataset_id}. "
            f"Values are provided as raw counts and {normalization}."
        ),
    }

    return ad.AnnData(
        X=None,
        obs=cell_metadata,
        var=gene_metadata,
        layers={
            "counts": counts,
            "normalized": log_cp10k(counts),
        },
        uns=metadata,
    )

def check_prepared_file(path, cell_ids, gene_ids, dataset_id, modality, normalization):
    """Check cell order, gene order and required benchmark metadata."""
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        metadata = read_elem(handle["uns"])

        if list(obs.index) != list(cell_ids):
            raise ValueError(f"Cell IDs or order do not match: {path}")
        if list(var.index) != list(gene_ids):
            raise ValueError(f"Gene IDs or order do not match: {path}")

        expected = {
            "dataset_id": dataset_id,
            "modality": modality,
            "normalization_id": normalization,
        }
        for key, value in expected.items():
            if metadata.get(key) != value:
                raise ValueError(f"Unexpected {key} in {path}")

        for layer in ("counts", "normalized"):
            if layer not in handle["layers"]:
                raise ValueError(f"Missing layer '{layer}' in {path}")
            matrix = sparse_dataset(handle["layers"][layer])
            if matrix.shape != (len(cell_ids), len(gene_ids)):
                raise ValueError(f"Unexpected shape for '{layer}' in {path}")


def check_prepared_values(path, raw_path, raw_layer, cell_ids, gene_ids):
    """Compare all stored values with their expected source and normalization."""
    with h5py.File(path, "r") as handle:
        stored_counts = sparse_dataset(handle["layers"]["counts"])
        stored_normalized = sparse_dataset(handle["layers"]["normalized"])

        for start in range(0, len(cell_ids), 512):
            stop = min(start + 512, len(cell_ids))

            expected_counts = read_raw_counts(
                raw_path, raw_layer, cell_ids[start:stop], gene_ids
            )
            actual_counts = stored_counts[start:stop].astype(np.float64)
            count_difference = actual_counts - expected_counts

            if not np.isfinite(actual_counts.data).all() or count_difference.nnz:
                raise ValueError(f"Stored counts differ from the raw source: {path}")

            expected_normalized = log_cp10k(expected_counts)
            actual_normalized = stored_normalized[start:stop]
            difference = actual_normalized - expected_normalized

            if (
                not np.isfinite(actual_normalized.data).all()
                or not np.allclose(difference.data, 0, rtol=0, atol=1e-6)
            ):
                raise ValueError(f"Stored normalization does not match log_cp10k: {path}")


def file_signature(path):
    """Identify a local file by its path, size and modification time."""
    path = Path(path).resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "size": stat.st_size,
        "modified_ns": stat.st_mtime_ns,
    }


def preparation_record(data_settings, split_path, prepared_dir, dataset_id):
    """Describe the inputs, settings and outputs of a preparation."""
    prepared_dir = Path(prepared_dir)
    filenames = [
        f"{part}_{mod}.h5ad"
        for part in ("train", "test")
        for mod in ("mod1", "mod2")
    ]

    return {
        "preparation_version": 1,
        "dataset_id": dataset_id,
        "settings": json.loads(json.dumps(data_settings, default=str)),
        "sources": {
            key: file_signature(data_settings[key])
            for key in ("raw_path", "filtered_path")
        },
        "split_sha256": hashlib.sha256(
            Path(split_path).read_bytes()
        ).hexdigest(),
        "outputs": {
            name: file_signature(prepared_dir / name)
            for name in filenames
        },
    }


def prepare_benchmark_data(data_settings, split_path, prepared_dir, dataset_id):
    """Prepare benchmark inputs or reuse a matching validated preparation."""
    prepared_dir = Path(prepared_dir)
    manifest_path = prepared_dir / "preparation_manifest.json"

    if data_settings["normalization"] != "log_cp10k":
        raise ValueError("Only log_cp10k is currently supported.")
    if data_settings["evaluation_split"] != "validation":
        raise ValueError("This development workflow evaluates validation cells only.")

    split = load_plate_split(split_path, data_settings["histone_mark"])
    obs, var, cells = select_cells_and_genes(
        data_settings["filtered_path"], split, data_settings["histone_mark"]
    )

    files = [
        ("train_mod1.h5ad", "train", "histone_layer", data_settings["histone_mark"]),
        ("train_mod2.h5ad", "train", "rna_layer", "GEX"),
        ("test_mod1.h5ad", "validation", "histone_layer", data_settings["histone_mark"]),
        ("test_mod2.h5ad", "validation", "rna_layer", "GEX"),
    ]

    reuse = False
    if manifest_path.exists():
        saved = json.loads(manifest_path.read_text())
        try:
            current = preparation_record(
                data_settings, split_path, prepared_dir, dataset_id
            )
        except FileNotFoundError as error:
            raise ValueError(
                "A registered file is missing. Use a new preparation directory "
                "or investigate the existing preparation."
            ) from error

        if saved != current:
            raise ValueError(
                "Files or settings differ from the registered preparation. "
                "Use a new preparation directory for changed inputs or settings."
            )
        reuse = True

    prepared_dir.mkdir(parents=True, exist_ok=True)

    for filename, partition, layer_key, modality in files:
        path = prepared_dir / filename
        cell_ids = cells[partition]

        if not reuse and not path.exists():
            print(f"Creating: {filename}", flush=True)

            obj = make_benchmark_adata(
                raw_path=data_settings["raw_path"],
                raw_layer=data_settings[layer_key],
                cell_ids=cell_ids,
                obs=obs,
                var=var,
                modality=modality,
                dataset_id=dataset_id,
                normalization=data_settings["normalization"],
            )

            temporary = path.with_suffix(".pending.h5ad")
            obj.write_h5ad(temporary, compression="lzf")
            del obj
            temporary.replace(path)

        check_prepared_file(
            path, cell_ids, var.index, dataset_id,
            modality, data_settings["normalization"],
        )

        if not reuse:
            print(f"Validating values: {filename}", flush=True)
            check_prepared_values(
                path, data_settings["raw_path"],
                data_settings[layer_key], cell_ids, var.index,
            )

    if not reuse:
        record = preparation_record(
            data_settings, split_path, prepared_dir, dataset_id
        )
        temporary_manifest = manifest_path.with_suffix(".pending.json")
        temporary_manifest.write_text(json.dumps(record, indent=2))
        temporary_manifest.replace(manifest_path)

    print(
        "Files already exist; matching validated preparation reused."
        if reuse else "Preparation completed and validated."
    )
    print(f"Dataset: {dataset_id}")
    print(f"Prepared data: {prepared_dir}")
    print(f"Split: {Path(split_path).name}")
    print("All five files present.")

    return prepared_dir


