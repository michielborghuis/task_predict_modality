import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

## VIASH START
par = {
    "input_raw": "raw.h5ad",
    "input_filtered": "filtered.h5ad",
    "input_annotations": "annotations.h5ad",
    "input_split": "split.csv",
    "mark": "H3K4me3",
    "annotation_sample_prefix": "PZ-MB-TChIC-",
    "output_mod1": "raw_mod1.h5ad",
    "output_mod2": "raw_mod2.h5ad",
}
## VIASH END

ANNOTATION_COLUMNS = ["celltype", "celltype_general", "lineage"]
RAW_LAYERS = {"H3K4me3": ("chic_counts", "H3K4me3"), "GEX": ("unspliced", "GEX")}


def _metadata(path):
    # Backed AnnData still eagerly reads layers. Read only the two dataframe
    # groups, preserving categorical dtypes, index names and row order.
    with h5py.File(path, "r") as source:
        obs = ad.io.read_elem(source["obs"])
        var = ad.io.read_elem(source["var"])
    obs_names = pd.Index(obs.index.astype(str))
    var_names = pd.Index(var.index.astype(str))
    if not obs_names.is_unique:
        raise ValueError(f"{path}: observation IDs are not unique")
    if not var_names.is_unique:
        raise ValueError(f"{path}: feature IDs are not unique")
    return obs_names, var_names, obs, var


def _subset_csr(group, row_indices, col_indices, shape):
    n_rows, n_cols = shape
    col_map = np.full(n_cols, -1, dtype=np.int64)
    col_map[np.asarray(col_indices, dtype=np.int64)] = np.arange(len(col_indices))
    indptr = group["indptr"][:]
    out_indices = []
    out_data = []
    out_indptr = [0]
    for source_row in row_indices:
        start, stop = int(indptr[source_row]), int(indptr[source_row + 1])
        if start != stop:
            source_cols = group["indices"][start:stop]
            values = group["data"][start:stop]
            mapped_cols = col_map[source_cols]
            keep = mapped_cols >= 0
            if np.any(keep):
                out_indices.append(mapped_cols[keep].astype(np.int32, copy=False))
                out_data.append(values[keep])
                out_indptr.append(out_indptr[-1] + int(keep.sum()))
            else:
                out_indptr.append(out_indptr[-1])
        else:
            out_indptr.append(out_indptr[-1])
    indices = np.concatenate(out_indices) if out_indices else np.empty(0, dtype=np.int32)
    values = np.concatenate(out_data) if out_data else np.empty(0, dtype=np.float64)
    return sp.csr_matrix(
        (values, indices, np.asarray(out_indptr, dtype=np.int64)),
        shape=(len(row_indices), len(col_indices)),
    )


def _subset_csc(group, row_indices, col_indices, shape):
    n_rows, n_cols = shape
    row_map = np.full(n_rows, -1, dtype=np.int64)
    row_map[np.asarray(row_indices, dtype=np.int64)] = np.arange(len(row_indices))
    indptr = group["indptr"][:]
    out_rows = []
    out_cols = []
    out_data = []
    for out_col, source_col in enumerate(col_indices):
        start, stop = int(indptr[source_col]), int(indptr[source_col + 1])
        if start == stop:
            continue
        source_rows = group["indices"][start:stop]
        values = group["data"][start:stop]
        mapped_rows = row_map[source_rows]
        keep = mapped_rows >= 0
        if np.any(keep):
            mapped_rows = mapped_rows[keep]
            kept_values = values[keep]
            out_rows.extend(mapped_rows.tolist())
            out_cols.extend([out_col] * len(mapped_rows))
            out_data.extend(kept_values.tolist())
    return sp.coo_matrix(
        (np.asarray(out_data), (np.asarray(out_rows), np.asarray(out_cols))),
        shape=(len(row_indices), len(col_indices)),
    ).tocsr()


def _read_count_subset(handle, layer, row_indices, col_indices):
    obj = handle["layers"][layer]
    if isinstance(obj, h5py.Dataset):
        values = np.vstack([obj[int(row), :][col_indices] for row in row_indices])
        matrix = sp.csr_matrix(values)
    else:
        shape = tuple(int(v) for v in obj.attrs["shape"])
        encoding = obj.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8")
        if encoding == "csr_matrix":
            matrix = _subset_csr(obj, row_indices, col_indices, shape)
        elif encoding == "csc_matrix":
            matrix = _subset_csc(obj, row_indices, col_indices, shape)
        else:
            raise ValueError(f"Unsupported sparse encoding for layer {layer}: {encoding}")
    values = matrix.data
    if values.size and (not np.isfinite(values).all() or np.any(values < 0)):
        raise ValueError(f"Layer {layer} contains negative or non-finite counts")
    if values.size and not np.allclose(values, np.rint(values), rtol=0, atol=1e-6):
        raise ValueError(f"Layer {layer} contains non-integer values")
    if values.size and np.max(values) > np.iinfo(np.int32).max:
        raise ValueError(f"Layer {layer} contains values outside int32 range")
    matrix.data = np.rint(values).astype(np.int32)
    matrix.eliminate_zeros()
    return matrix


def _annotation_rows(annotation_path, prefix, mark, sample_ids, barcodes):
    ann_names, _, annotations, _ = _metadata(annotation_path)
    del ann_names
    required = {"sample", "mark", "barcode", *ANNOTATION_COLUMNS}
    missing = required.difference(annotations.columns)
    if missing:
        raise ValueError(f"Annotation file is missing columns: {sorted(missing)}")
    ann = annotations.copy()
    if ann[["sample", "mark", "barcode", *ANNOTATION_COLUMNS]].isna().any().any():
        raise ValueError("Annotation source has missing key or annotation values")
    ann["sample"] = ann["sample"].astype(str)
    if prefix:
        not_prefixed = ~ann["sample"].str.startswith(prefix)
        if not_prefixed.any():
            raise ValueError(f"Annotation sample IDs must start with {prefix!r}")
        ann["sample"] = ann["sample"].str[len(prefix):]
    ann["mark"] = ann["mark"].astype(str)
    ann["barcode"] = ann["barcode"].astype(str)
    keys = list(zip(ann["sample"], ann["mark"], ann["barcode"]))
    if len(set(keys)) != len(keys):
        raise ValueError("Annotation source contains duplicate (sample, mark, barcode) keys")
    lookup = {key: idx for idx, key in enumerate(keys)}
    requested = list(zip(sample_ids.astype(str), [mark] * len(sample_ids), barcodes.astype(str)))
    if len(set(requested)) != len(requested):
        raise ValueError("Retained cells contain duplicate (sample, mark, barcode) keys")
    missing_keys = [key for key in requested if key not in lookup]
    if missing_keys:
        raise ValueError(
            f"Missing cell annotation for {len(missing_keys)} retained cells; "
            f"first missing key is {missing_keys[0]}"
        )
    selected = ann.iloc[[lookup[key] for key in requested]][ANNOTATION_COLUMNS].reset_index(drop=True)
    for column in ANNOTATION_COLUMNS:
        if selected[column].isna().any():
            raise ValueError(f"Annotation column {column!r} contains missing values for retained cells")
    return selected


def build_inputs(params):
    mark = params["mark"]
    if mark != "H3K4me3":
        raise ValueError("This workflow currently supports only mark H3K4me3")
    if mark not in RAW_LAYERS:
        raise ValueError(f"Unsupported histone mark: {mark}")

    filtered_cells, filtered_genes, filtered_obs, filtered_var = _metadata(params["input_filtered"])
    raw_cells, raw_genes, _, _ = _metadata(params["input_raw"])
    if not {"sample", "mark", "barcode_nla"}.issubset(filtered_obs.columns):
        raise ValueError("Filtered obs must contain sample, mark, and barcode_nla")

    split = pd.read_csv(params["input_split"], dtype={"sample": str, "mark": str, "split": str})
    if not {"sample", "mark", "split"}.issubset(split.columns):
        raise ValueError("Split manifest must contain sample, mark, and split columns")
    split = split.loc[split["mark"] == mark, ["sample", "mark", "split"]].copy()
    if split.empty:
        raise ValueError(f"Split manifest has no assignments for {mark}")
    if split.duplicated(["sample", "mark"]).any():
        raise ValueError("Split manifest contains duplicate sample/mark assignments")
    if not set(split["split"]).issubset({"train", "validation", "test"}):
        raise ValueError("Split manifest labels must be train, validation, or test")

    filtered_obs = filtered_obs.copy()
    filtered_obs["sample"] = filtered_obs["sample"].astype(str)
    filtered_obs["mark"] = filtered_obs["mark"].astype(str)
    selected_mask = filtered_obs["mark"].eq(mark)
    selected_samples = set(filtered_obs.loc[selected_mask, "sample"])
    manifest_samples = set(split["sample"])
    if selected_samples != manifest_samples:
        missing = sorted(selected_samples - manifest_samples)
        extra = sorted(manifest_samples - selected_samples)
        raise ValueError(f"Split/sample mismatch; missing assignments={missing[:5]}, unused manifest samples={extra[:5]}")

    assignment = dict(zip(split["sample"], split["split"]))
    labels = filtered_obs.loc[selected_mask, "sample"].map(assignment)
    development = labels.isin(["train", "validation"])
    filtered_positions = np.flatnonzero(selected_mask.to_numpy())[development.to_numpy()]
    selected_labels = labels.loc[development].astype(str).to_numpy()
    selected_cells = filtered_cells.take(filtered_positions).astype(str)
    selected_samples = filtered_obs.iloc[filtered_positions]["sample"]
    selected_barcodes = filtered_obs.iloc[filtered_positions]["barcode_nla"]
    if selected_samples.isna().any() or selected_barcodes.isna().any():
        raise ValueError("Retained cell metadata contains missing sample/barcode values")
    selected_samples = selected_samples.astype(str)
    selected_barcodes = selected_barcodes.astype(str)

    raw_cell_lookup = {cell: idx for idx, cell in enumerate(raw_cells.astype(str))}
    raw_gene_lookup = {gene: idx for idx, gene in enumerate(raw_genes.astype(str))}
    missing_cells = [cell for cell in selected_cells if cell not in raw_cell_lookup]
    missing_genes = [gene for gene in filtered_genes.astype(str) if gene not in raw_gene_lookup]
    if missing_cells:
        raise ValueError(f"{len(missing_cells)} selected filtered cells are absent from raw; first: {missing_cells[0]}")
    if missing_genes:
        raise ValueError(f"{len(missing_genes)} retained filtered features are absent from raw; first: {missing_genes[0]}")
    raw_cell_positions = np.asarray([raw_cell_lookup[cell] for cell in selected_cells], dtype=np.int64)
    raw_gene_positions = np.asarray([raw_gene_lookup[gene] for gene in filtered_genes.astype(str)], dtype=np.int64)

    annotation = _annotation_rows(
        params["input_annotations"],
        params["annotation_sample_prefix"],
        mark,
        selected_samples,
        selected_barcodes,
    )
    obs = filtered_obs.iloc[filtered_positions].copy()
    obs.index = selected_cells
    obs["batch"] = selected_samples.to_numpy()
    obs["is_train"] = selected_labels
    obs["cell_type"] = annotation["celltype"].astype(str).to_numpy()
    for column in ("celltype_general", "lineage"):
        obs[column] = annotation[column].to_numpy()

    var = filtered_var.copy()
    var.index = filtered_genes.astype(str)
    var["feature_id"] = var.index.astype(str)
    var["gene_ids"] = var["feature_id"]
    if "Symbol" in var.columns:
        var["feature_name"] = var["Symbol"].astype(str)
    else:
        var["feature_name"] = var["feature_id"]

    dataset_id = "tchic_mgast_h3k4me3_genelevel_dev"
    uns = {
        "dataset_id": dataset_id,
        "dataset_name": "T-ChIC mouse gastruloid H3K4me3 development split",
        "dataset_summary": "Gene-level H3K4me3 and unspliced RNA from paired single-cell T-ChIC.",
        "dataset_description": "Development preparation uses the fixed plate-level training and validation partitions. The final test partition is excluded.",
        "dataset_organism": "Mus musculus",
        "normalization_id": "raw_counts",
    }
    raw_layer = RAW_LAYERS[mark][0]
    with h5py.File(params["input_raw"], "r") as raw_handle:
        mod1_counts = _read_count_subset(raw_handle, raw_layer, raw_cell_positions, raw_gene_positions)
        mod2_counts = _read_count_subset(raw_handle, "unspliced", raw_cell_positions, raw_gene_positions)

    outputs = []
    for path, modality, counts in [
        (params["output_mod1"], mark, mod1_counts),
        (params["output_mod2"], "GEX", mod2_counts),
    ]:
        data = ad.AnnData(X=None, obs=obs.copy(), var=var.copy(), layers={"counts": counts}, uns={**uns, "modality": modality})
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        data.write_h5ad(path, compression="gzip")
        outputs.append(data)
    if not outputs[0].obs_names.equals(outputs[1].obs_names):
        raise RuntimeError("Extracted modality cell IDs are misaligned")
    return outputs


def main(params):
    build_inputs(params)


if __name__ == "__main__":
    main(par)
