from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from datetime import datetime, timezone
from unittest import mock
import os
import subprocess

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
N_VALIDATION = 1002


def create_fixture(folder: Path) -> dict[str, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    cells = ["train_a", "train_b"] + [f"validation_{i:04d}" for i in range(N_VALIDATION)] + ["final_holdout"]
    sample_for = {"train_a": "sample_train", "train_b": "sample_train", "final_holdout": "sample_final"}
    sample_for.update({cell: "sample_validation" for cell in cells if cell.startswith("validation_")})
    all_barcodes = {cell: f"barcode_{cell}" for cell in cells}
    rng = np.random.default_rng(42)
    raw_cells = cells + ['raw_only']
    sample_for['raw_only'] = 'sample_extra'
    raw_order = rng.permutation(raw_cells).tolist()
    all_barcodes['train_b'] = all_barcodes['validation_0001']
    filtered_order = cells[::-1]
    raw_features = ["gene_c", "gene_a", "gene_b"]
    retained_features = ["gene_b", "gene_a"]

    chic_rows = []
    rna_rows = []
    for cell in raw_order:
        if cell == "validation_0000":
            chic_rows.append([0, 0, 0])
            rna_rows.append([0, 0, 0])
        else:
            cell_num = raw_cells.index(cell) + 1
            chic_rows.append([5, cell_num % 11 + 1, cell_num % 7 + 2])
            rna_rows.append([3, cell_num % 5 + 2, cell_num % 13 + 1])
    raw_obs = pd.DataFrame({"sample": [sample_for[c] for c in raw_order]}, index=raw_order)
    raw_var = pd.DataFrame(index=raw_features)
    raw = ad.AnnData(
        X=sp.csr_matrix((len(raw_order), len(raw_features)), dtype=np.int32),
        obs=raw_obs,
        var=raw_var,
        layers={
            "chic_counts": sp.csr_matrix(np.asarray(chic_rows, dtype=np.int32)),
            "unspliced": sp.csr_matrix(np.asarray(rna_rows, dtype=np.int32)),
        },
    )
    raw_path = folder / "raw.h5ad"
    raw.write_h5ad(raw_path, compression="gzip")

    filtered_obs = pd.DataFrame(
        {
            "sample": [sample_for[c] for c in filtered_order],
            "mark": pd.Categorical(["H3K4me3"] * len(filtered_order), categories=["H3K27me3", "H3K4me3"], ordered=True),
            "barcode_nla": [all_barcodes[c] for c in filtered_order],
            "day": pd.Categorical(["day1" if sample_for[c] == "sample_train" else "day2" for c in filtered_order], categories=["day2", "day1"], ordered=True),
            "plate_batch": [1 if sample_for[c] == "sample_train" else 2 if sample_for[c] == "sample_validation" else 3 for c in filtered_order],
        },
        index=filtered_order,
    )
    filtered_var = pd.DataFrame(
        {"Chromosome": ["chr2", "chr1"], "Symbol": ["B", "A"], "Biotype": ["protein_coding"] * 2},
        index=retained_features,
    )
    filtered = ad.AnnData(
        X=sp.csr_matrix((len(filtered_order), len(retained_features)), dtype=np.int32),
        obs=filtered_obs,
        var=filtered_var,
    )
    filtered_path = folder / "filtered.h5ad"
    filtered.write_h5ad(filtered_path, compression="gzip")

    annotation_rows = []
    for cell in cells:
        sample = sample_for[cell]
        annotation_rows.append(
            {
                "sample": "PZ-MB-TChIC-" + sample,
                "mark": "H3K4me3",
                "barcode": all_barcodes[cell],
                "celltype": "type_a" if cell == "train_a" else "type_b",
                "celltype_general": "gastruloid",
                "lineage": "mesoderm",
            }
        )
    rng.shuffle(annotation_rows)
    annotation_rows.append({**annotation_rows[0], 'mark': 'H3K27me3', 'celltype': 'wrong_mark'})
    annotation_obs = pd.DataFrame(annotation_rows, index=[f"ann_{i}" for i in range(len(annotation_rows))])
    annotation = ad.AnnData(
        X=sp.csr_matrix((len(annotation_rows), 1), dtype=np.int8),
        obs=annotation_obs,
        var=pd.DataFrame(index=["annotation_feature"]),
    )
    annotation_path = folder / "annotations.h5ad"
    annotation.write_h5ad(annotation_path, compression="gzip")

    split_path = folder / "split.csv"
    pd.DataFrame(
        [
            {"sample": "sample_train", "mark": "H3K4me3", "split": "train"},
            {"sample": "sample_validation", "mark": "H3K4me3", "split": "validation"},
            {"sample": "sample_final", "mark": "H3K4me3", "split": "test"},
        ]
    ).to_csv(split_path, index=False)
    return {"raw": raw_path, "filtered": filtered_path, "annotations": annotation_path, "split": split_path}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--components-only", action="store_true")
    parser.add_argument("--target", default="target_tchic_verified_20261023")
    parser.add_argument("--output-dir", type=Path, required=True, help="new directory; refuses to reuse outputs")
    args = parser.parse_args()
    target = (ROOT / args.target).resolve()
    scratch = args.output_dir.resolve()
    scratch.mkdir(parents=True, exist_ok=False)
    lock = json.loads((target / "runtime/runtime.json").read_text())
    dep = ROOT / ".dependencies/datasets-106edbde"
    pin = "106edbdeae9775d71bdfb4e2dcca67e187dc843b"
    assert subprocess.check_output(["git", "--no-optional-locks", "-C", str(dep), "rev-parse", "HEAD"], text=True).strip() == pin
    subprocess.run(["git", "--no-optional-locks", "-C", str(dep), "diff", "--quiet", "HEAD", "--", "src/normalization/log_cp", "src/api", "_viash.yaml"], check=True)
    for key in ["python", "processor"]:
        assert sha256(lock[key + "_image"]) == lock[key + "_sha256"], f"Changed {key} image"
    for path, expected in lock["source_sha256"].items():
        assert sha256(ROOT / path) == expected, f"Rebuild changed source: {path}"
    for path, expected in lock["artifact_sha256"].items():
        assert sha256(target / path) == expected, f"Changed generated artifact: {path}"
    (scratch / "runtime.json").write_text(json.dumps(lock, indent=2))
    (scratch / "harness.py").write_bytes(Path(__file__).read_bytes())
    env = os.environ.copy()
    env.update(NXF_ANSI_LOG="false", NXF_SYNTAX_PARSER="v1",
               NXF_WORK=str(scratch / "nextflow-work"), VIASH_TEMP=str(scratch / "viash-temp"),
               TCHIC_PYTHON_IMAGE=lock["python_image"], TCHIC_PROCESSOR_IMAGE=lock["processor_image"])
    for directory in ["cache", "cache/numba", "cache/matplotlib", "logs", "viash-temp"]:
        (scratch / directory).mkdir(parents=True, exist_ok=True)
    commands = []

    def invoke(command, expect_success=True):
        number = len(commands) + 1
        result = subprocess.run([str(x) for x in command], cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = scratch / "logs" / f"{number:02d}.log"
        log.write_text(result.stdout)
        commands.append({"command": [str(x) for x in command], "exit_code": result.returncode, "log": str(log)})
        (scratch / "commands.json").write_text(json.dumps(commands, indent=2))
        print(f"command {number}: exit {result.returncode}; {log}", flush=True)
        if expect_success and result.returncode:
            raise RuntimeError(f"Command failed: {command}\n{result.stdout}")
        return result

    def container(image, command):
        return ["apptainer", "exec", "--cleanenv", "--no-home", "--bind", str(ROOT), "--bind", str(scratch),
                "--env", f"MPLCONFIGDIR={scratch}/cache/matplotlib", "--env", f"NUMBA_CACHE_DIR={scratch}/cache/numba",
                "--env", f"XDG_CACHE_HOME={scratch}/cache", image, *[str(x) for x in command]]

    python_image, r_image = lock["python_image"], lock["processor_image"]
    paths = create_fixture(scratch / "fixture")
    source_raw = ad.read_h5ad(paths["raw"])
    source_filtered = ad.read_h5ad(paths["filtered"])
    source_ann = ad.read_h5ad(paths["annotations"]).obs
    expected_cells = [cell for cell in source_filtered.obs_names if cell != "final_holdout"]
    features = list(source_filtered.var_names)
    train_cells = [c for c in expected_cells if c.startswith("train_")]
    validation_cells = [c for c in expected_cells if c.startswith("validation_")]
    annotation_lookup = source_ann.set_index(["sample", "mark", "barcode"])

    # Verify actual metadata-reader accesses, not merely its returned shape.
    module = runpy.run_path(str(ROOT / "src/data_processors/prepare_tchic_input/script.py"))
    original_read = ad.io.read_elem
    accessed = []
    def metadata_element(element, *pos, **kw):
        accessed.append(element.name)
        assert element.name in {"/obs", "/var"}, element.name
        return original_read(element, *pos, **kw)
    with mock.patch.object(ad, "read_h5ad", side_effect=AssertionError("must not load AnnData")), mock.patch.object(ad.io, "read_elem", side_effect=metadata_element):
        _, _, obs, var = module["_metadata"](paths["filtered"])
    assert accessed == ["/obs", "/var"]
    pd.testing.assert_frame_equal(obs, source_filtered.obs)
    pd.testing.assert_frame_equal(var, source_filtered.var)

    def reference_counts(names, modality):
        layer = "chic_counts" if modality == "mod1" else "unspliced"
        rows = source_raw.obs_names.get_indexer(names)
        cols = source_raw.var_names.get_indexer(features)
        assert (rows >= 0).all() and (cols >= 0).all()
        return source_raw.layers[layer][rows, :][:, cols].toarray()

    def check_data(path, names, modality, normalized, hvg=False):
        out = ad.read_h5ad(path)
        assert list(out.obs_names) == list(names)
        assert list(out.var_names) == features
        assert "final_holdout" not in out.obs_names and "raw_only" not in out.obs_names
        expected = reference_counts(names, modality)
        np.testing.assert_array_equal(out.layers["counts"].toarray(), expected)
        assert np.issubdtype(out.layers["counts"].dtype, np.integer)
        expected_obs = source_filtered.obs.loc[names]
        for column in ["sample", "day"]:
            assert list(out.obs[column].astype(str)) == list(expected_obs[column].astype(str)), column
        assert list(out.obs["batch"].astype(str)) == list(expected_obs["sample"].astype(str))
        expected_annotations = [annotation_lookup.loc[("PZ-MB-TChIC-" + str(row["sample"]), "H3K4me3", str(row["barcode_nla"]))]
                                for _, row in expected_obs.iterrows()]
        for target_column, source_column in [("cell_type", "celltype"), ("celltype_general", "celltype_general"), ("lineage", "lineage")]:
            assert list(out.obs[target_column].astype(str)) == [str(row[source_column]) for row in expected_annotations], target_column
        for column in ["Chromosome", "Symbol", "Biotype"]:
            assert list(out.var[column].astype(str)) == list(source_filtered.var[column].astype(str)), column
        for column in ["feature_id", "gene_ids"]:
            assert list(out.var[column].astype(str)) == features
        assert list(out.var["feature_name"].astype(str)) == list(source_filtered.var["Symbol"].astype(str))
        assert out.uns["dataset_id"] == "tchic_mgast_h3k4me3_genelevel_dev"
        assert out.uns["modality"] == ("H3K4me3" if modality == "mod1" else "GEX")
        assert out.uns["normalization_id"] == ("log_cp10k" if normalized else "raw_counts")
        if hvg:
            np.testing.assert_array_equal(out.var["hvg"].to_numpy(), [True, False] if modality == "mod1" else [False, True])
            np.testing.assert_allclose(out.var["hvg_score"].to_numpy(), [9.5, 0.2] if modality == "mod1" else [0.3, 8.5])
        else:
            assert "hvg" not in out.var and "hvg_score" not in out.var
        if normalized:
            assert np.issubdtype(out.layers["normalized"].dtype, np.floating)
            row_sums = expected.sum(axis=1)
            ref = np.zeros_like(expected, dtype=float)
            nonzero = row_sums > 0
            ref[nonzero] = np.log1p(10000 * expected[nonzero] / row_sums[nonzero, None])
            np.testing.assert_allclose(out.layers["normalized"].toarray(), ref, rtol=1e-6, atol=1e-6)
            np.testing.assert_allclose(out.obs["size_factors"].to_numpy(), row_sums / 10000, rtol=1e-6, atol=1e-7)
        return out

    stages = scratch / "stages"
    stages.mkdir()
    extractor = target / "executable/data_processors/prepare_tchic_input/prepare_tchic_input"
    base_args = [extractor, "--input_raw", paths["raw"], "--input_filtered", paths["filtered"],
                 "--input_annotations", paths["annotations"], "--input_split", paths["split"],
                 "--output_mod1", stages / "raw_mod1.h5ad", "--output_mod2", stages / "raw_mod2.h5ad"]
    invoke(container(python_image, base_args))
    for modality in ["mod1", "mod2"]:
        out = check_data(stages / f"raw_{modality}.h5ad", expected_cells, modality, False)
        assert list(out.obs["is_train"].astype(str)) == ["train" if c.startswith("train_") else "validation" for c in expected_cells]
        pd.testing.assert_series_equal(out.obs["day"], source_filtered.obs.loc[expected_cells, "day"])

    log_cp = ROOT / ".dependencies/datasets-106edbde/target_viash097/executable/normalization/log_cp/log_cp"
    for modality in ["mod1", "mod2"]:
        invoke(container(python_image, [log_cp, "--input", stages / f"raw_{modality}.h5ad", "--output", stages / f"normalized_{modality}.h5ad",
               "--n_cp", "10000", "--layer_output", "normalized", "--normalization_id", "log_cp10k", "--obs_size_factors", "size_factors"]))
        check_data(stages / f"normalized_{modality}.h5ad", expected_cells, modality, True)
    assert not np.array_equal(reference_counts(expected_cells, "mod1").sum(1), reference_counts(expected_cells, "mod2").sum(1))

    # Each malformed source is a separate synthetic input, preserving the fixture.
    for kind in ["duplicate_source", "missing_source", "duplicate_receiver"]:
        bad_path = scratch / f"{kind}.h5ad"
        bad_args = base_args.copy()
        if kind == "duplicate_receiver":
            bad = source_filtered.copy()
            bad.obs["barcode_nla"] = bad.obs["barcode_nla"].astype(str)
            bad.obs.loc["validation_0002", "barcode_nla"] = bad.obs.loc["validation_0001", "barcode_nla"]
            arg_name = "--input_filtered"
            expected_error = "Retained cells contain duplicate (sample, mark, barcode)"
        else:
            ann = ad.read_h5ad(paths["annotations"])
            if kind == "duplicate_source":
                duplicate = ann[[0]].copy()
                duplicate.obs_names = ["duplicate_record"]
                bad = ad.concat([ann, duplicate])
                expected_error = "Annotation source contains duplicate (sample, mark, barcode)"
            else:
                keep = ~((ann.obs["barcode"] == "barcode_validation_0001") & (ann.obs["sample"] == "PZ-MB-TChIC-sample_validation") & (ann.obs["mark"] == "H3K4me3"))
                bad = ann[keep].copy()
                expected_error = "Missing cell annotation"
            arg_name = "--input_annotations"
        bad.write_h5ad(bad_path)
        bad_args[bad_args.index(arg_name) + 1] = bad_path
        result = invoke(container(python_image, bad_args), expect_success=False)
        assert result.returncode != 0 and expected_error in result.stdout
    print("Extractor, metadata-only reads, annotation failures and independent normalization checks passed", flush=True)
    if args.components_only:
        return

    def schemas(folder):
        code = '''from pathlib import Path
import sys, yaml, importlib.metadata
from openproblems.project.component_tests.run_and_check_output import check_format
root, outputs = map(Path, sys.argv[1:])
print("OpenProblems validator version:", importlib.metadata.version("openproblems"))
for slot in ["train_mod1", "train_mod2", "test_mod1", "test_mod2"]:
    arg = yaml.safe_load((root / "src/api" / ("file_" + slot + ".yaml")).read_text())
    arg.update(name="--" + slot, clean_name=slot, value=str(outputs / (slot + ".h5ad")))
    check_format(arg)
print("All four local task schemas passed")
'''
        invoke(container(python_image, ["python", "-c", code, ROOT, folder]))

    def outputs_check(folder, default_limit=False, hvg=False):
        tests = []
        for modality in ["mod1", "mod2"]:
            check_data(folder / f"train_{modality}.h5ad", train_cells, modality, True, hvg)
            actual_names = list(ad.read_h5ad(folder / f"test_{modality}.h5ad").obs_names)
            if default_limit:
                assert len(actual_names) == 1000 and len(set(actual_names)) == 1000
                assert set(actual_names).issubset(validation_cells)
            else:
                assert actual_names == validation_cells
            check_data(folder / f"test_{modality}.h5ad", actual_names, modality, True, hvg)
            tests.append(actual_names)
        assert tests[0] == tests[1]
        schemas(folder)

    processor = target / "executable/data_processors/process_dataset/process_dataset"
    for case in ["default_1000_no_hvg", "unlimited_with_hvg"]:
        folder = scratch / case
        folder.mkdir()
        inputs = []
        for modality in ["mod1", "mod2"]:
            path = stages / f"normalized_{modality}.h5ad"
            if case == "unlimited_with_hvg":
                data = ad.read_h5ad(path)
                data.var["hvg"] = [True, False] if modality == "mod1" else [False, True]
                data.var["hvg_score"] = [9.5, 0.2] if modality == "mod1" else [0.3, 8.5]
                path = folder / f"input_{modality}.h5ad"
                data.write_h5ad(path)
            inputs.extend(["--input_" + modality, path])
        command = [processor, *inputs, "--preserve_obs_columns", "batch,size_factors,sample,day,cell_type,celltype_general,lineage",
                   "--preserve_var_columns", "gene_ids,hvg,hvg_score,feature_id,feature_name,Chromosome,Symbol,Biotype"]
        for slot in ["train_mod1", "train_mod2", "test_mod1", "test_mod2"]:
            command.extend(["--output_" + slot, folder / (slot + ".h5ad")])
        if case == "unlimited_with_hvg":
            command.extend(["--max_test_cells", "-1"])
        invoke(container(r_image, command))
        outputs_check(folder, default_limit=case.startswith("default"), hvg=case.endswith("with_hvg"))

    published = scratch / "published"
    command = ["nextflow", "-log", scratch / "nextflow.log", "run", target / "nextflow/workflows/prepare_tchic",
               "-profile", "singularity", "-c", ROOT / "src/workflows/prepare_tchic/runtime.config", "--publish_dir", published]
    for key, path in paths.items():
        command.extend(["--input_" + key, path])
    for slot in ["train_mod1", "train_mod2", "test_mod1", "test_mod2"]:
        command.extend(["--output_" + slot, slot + ".h5ad"])
    command.extend(["--max_test_cells", "-1"])
    invoke(command)
    assert {p.name for p in published.glob("*.h5ad")} == {s + ".h5ad" for s in ["train_mod1", "train_mod2", "test_mod1", "test_mod2"]}
    outputs_check(published)
    result = {"passed": True, "finished_utc": datetime.now(timezone.utc).isoformat(), "published": str(published),
              "workflow": "non-stub", "train_cells": len(train_cells), "validation_cells": len(validation_cells),
              "processor_default_test_cells": 1000, "processor_unlimited_test_cells": 1002,
              "hvg_present_and_absent": True, "metadata_only_reader": True, "annotation_failures": 3,
              "schema_validator": "openproblems.project.component_tests.run_and_check_output.check_format",
              "files": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in published.glob("*.h5ad")}}
    (scratch / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
