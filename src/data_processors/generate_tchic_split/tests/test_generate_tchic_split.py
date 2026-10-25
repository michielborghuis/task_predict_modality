import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import anndata
import numpy as np
import pandas as pd
from scipy import sparse


SCRIPT_PATH = Path(__file__).parents[1] / "script.py"
SPEC = importlib.util.spec_from_file_location("generate_tchic_split", SCRIPT_PATH)
split = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(split)


def make_settings():
    settings_path = SCRIPT_PATH.parent / "h3k4me3_split_settings.json"
    return split.load_settings(settings_path) | {
        "n_iterations": 250,
    }


def make_plates():
    rows = []
    plate_batch = 0
    for day, number in [(0, 6), (3, 6), (4, 12), (5, 12), (6, 12), (7, 12)]:
        for within_day in range(number):
            rows.append(
                {
                    "sample": f"plate-{day}-{within_day:02}",
                    "n_cells": 20 + (plate_batch % 7),
                    "plate_batch": plate_batch,
                    "day": day,
                    "rep": f"rep{within_day % 3}",
                    "mark": "H3K4me3",
                    "plate_grp": within_day + 1,
                    "seqrun": f"run{within_day % 4}",
                    "exp_batch": "batchA" if day < 5 else "batchB",
                    "day_rep": f"{day}_rep{within_day % 3}",
                    "batch_day_rep": f"{day}_batchA_rep{within_day % 3}",
                }
            )
            plate_batch += 1
    return pd.DataFrame(rows)


class SplitGeneratorTests(unittest.TestCase):
    def test_repeatability_and_constraints(self):
        settings = make_settings()
        plates = make_plates()
        first, first_scores, first_valid = split.search_best_split(plates, settings)
        second, second_scores, second_valid = split.search_best_split(plates, settings)

        self.assertEqual(
            first[["sample", "split"]].to_dict("records"),
            second[["sample", "split"]].to_dict("records"),
        )
        self.assertEqual(first_scores, second_scores)
        self.assertEqual(first_valid, second_valid)
        self.assertEqual(first["split"].value_counts().to_dict(), {
            "train": 42,
            "validation": 9,
            "test": 9,
        })
        self.assertEqual(len(first), 60)
        self.assertTrue(first["sample"].is_unique)
        split.validate_assignment(
            plates.sort_values("sample", kind="stable").reset_index(drop=True),
            first.sort_values("sample", kind="stable")["split"].cat.codes.to_numpy(),
            settings,
        )
        self.assertGreater(first_valid, 0)

    def test_categorical_sample_order_need_not_match_candidate_day_order(self):
        settings = make_settings()
        plates = make_plates()
        reordered = pd.concat(
            [
                plates.loc[plates["day"] == day]
                for day in [3, 4, 5, 6, 7, 0]
            ],
            ignore_index=True,
        )

        result, _, valid_candidates = split.search_best_split(reordered, settings)

        self.assertEqual(
            result["split"].value_counts().to_dict(),
            {"train": 42, "validation": 9, "test": 9},
        )
        self.assertGreater(valid_candidates, 0)

    def test_missing_and_inconsistent_metadata_fail_clearly(self):
        settings = make_settings()
        obs = pd.DataFrame(
            {
                "sample": ["p1", "p1", "p2"],
                "plate_batch": [1, 2, 3],
                "day": [0, 0, 0],
                "rep": ["r1", "r1", "r1"],
                "mark": ["H3K4me3"] * 3,
                "plate_grp": [1, 1, 2],
                "seqrun": ["s1"] * 3,
                "exp_batch": ["b1"] * 3,
                "day_rep": ["d0r1"] * 3,
                "batch_day_rep": ["d0b1r1"] * 3,
            }
        )
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            split.metadata_to_plate_table(obs, settings["mark"])
        with self.assertRaisesRegex(ValueError, "missing required metadata"):
            split.metadata_to_plate_table(obs.drop(columns="seqrun"), settings["mark"])

    def test_no_valid_candidate_fails_clearly(self):
        settings = make_settings()
        plates = make_plates()
        plates["rep"] = plates["sample"]
        with self.assertRaisesRegex(RuntimeError, "No valid candidate"):
            split.search_best_split(plates, settings | {"n_iterations": 3})

    def test_mapping_comparison_ignores_row_order(self):
        import sys

        compare_path = SCRIPT_PATH.parent / "compare_splits.py"
        compare_spec = importlib.util.spec_from_file_location(
            "compare_splits", compare_path
        )
        compare = importlib.util.module_from_spec(compare_spec)
        sys.modules[compare_spec.name] = compare
        compare_spec.loader.exec_module(compare)
        reference = pd.DataFrame(
            {"mark": ["H3K4me3", "H3K4me3"], "sample": ["a", "b"],
             "split": ["train", "test"], "other": [1, 2]}
        )
        reordered = reference.iloc[::-1].drop(columns="other")
        with tempfile.TemporaryDirectory() as temp_dir:
            reference_path = Path(temp_dir) / "reference.csv"
            generated_path = Path(temp_dir) / "generated.csv"
            reference.to_csv(reference_path, index=False)
            reordered.to_csv(generated_path, index=False)
            _, changed = compare.compare_mappings(reference_path, generated_path)
        self.assertTrue(changed.empty)

    def test_h5ad_outputs_and_provenance(self):
        settings = make_settings()
        plates = make_plates()
        obs = plates.drop(columns="n_cells").copy()
        sample_categories = obs["sample"].tolist()
        obs["sample"] = pd.Categorical(
            obs["sample"],
            categories=sample_categories,
            ordered=True,
        )
        obs["day"] = pd.Categorical(
            obs["day"].astype(str),
            categories=["0", "3", "4", "5", "6", "7"],
            ordered=True,
        )
        obs.index = [f"cell-{index:03}" for index in range(len(obs))]
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            metadata_path = root / "metadata.h5ad"
            settings_path = root / "settings.json"
            manifest_path = root / "manifest.csv"
            provenance_path = root / "provenance.json"
            anndata.AnnData(
                X=sparse.csr_matrix((len(obs), 1), dtype=np.float32),
                obs=obs,
            ).write_h5ad(metadata_path)
            loaded_obs, source_shape = split.read_h5ad_obs_metadata(metadata_path)
            self.assertEqual(source_shape, [60, 1])
            self.assertEqual(
                list(loaded_obs["sample"].cat.categories),
                sample_categories,
            )
            self.assertTrue(loaded_obs["sample"].cat.ordered)
            self.assertEqual(
                list(loaded_obs["day"].cat.categories),
                ["0", "3", "4", "5", "6", "7"],
            )
            settings_path.write_text(
                json.dumps(
                    {
                        key: value
                        for key, value in settings.items()
                        if not key.startswith("_")
                    }
                ),
                encoding="utf-8",
            )
            original_par = split.par
            split.par = {
                "input_metadata": str(metadata_path),
                "input_settings": str(settings_path),
                "output_manifest": str(manifest_path),
                "output_provenance": str(provenance_path),
            }
            try:
                with mock.patch.object(
                    split.anndata,
                    "read_h5ad",
                    side_effect=AssertionError("read_h5ad must not load the matrix"),
                ):
                    split.main()
            finally:
                split.par = original_par

            result = pd.read_csv(manifest_path)
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            self.assertEqual(len(result), 60)
            self.assertTrue(result["sample"].is_unique)
            self.assertEqual(
                provenance["manifest_sha256"],
                hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                provenance["metadata_fingerprint_sha256"].__len__(), 64
            )
            self.assertIn("selected_scores", provenance["search"])
            self.assertIn("software_versions", provenance)
            self.assertIn(
                "final holdout",
                provenance["manifest_test_label"].lower(),
            )
            self.assertEqual(provenance["settings"]["seed"], 42)
            self.assertGreater(provenance["search"]["valid_candidates"], 0)


if __name__ == "__main__":
    unittest.main()
