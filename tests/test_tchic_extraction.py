import runpy
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
extract = runpy.run_path(
    str(ROOT / "src/data_processors/prepare_tchic_input/script.py")
)["build_inputs"]
create_fixture = runpy.run_path(
    str(ROOT / "tests/run_tchic_synthetic_workflow.py")
)["create_fixture"]


class ExtractionTests(unittest.TestCase):
    def test_counts_and_metadata(self):
        for encoding in ("csr", "csc", "dense"):
            with self.subTest(encoding=encoding), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                paths = create_fixture(folder / "inputs")
                raw = ad.read_h5ad(paths["raw"])
                filtered = ad.read_h5ad(paths["filtered"])
                annotations = ad.read_h5ad(paths["annotations"]).obs
                annotations = annotations.set_index(["sample", "mark", "barcode"])

                for layer in ("chic_counts", "unspliced"):
                    values = raw.layers[layer].toarray()
                    raw.layers[layer] = (
                        values if encoding == "dense"
                        else getattr(sp, f"{encoding}_matrix")(values)
                    )
                raw.write_h5ad(paths["raw"])

                extract({
                    **{f"input_{key}": str(path) for key, path in paths.items()},
                    "mark": "H3K4me3",
                    "output_mod1": str(folder / "mod1.h5ad"),
                    "output_mod2": str(folder / "mod2.h5ad"),
                })

                cells = filtered.obs_names[filtered.obs_names != "final_holdout"]
                genes = filtered.var_names
                rows = raw.obs_names.get_indexer(cells)
                cols = raw.var_names.get_indexer(genes)

                for number, layer in [(1, "chic_counts"), (2, "unspliced")]:
                    out = ad.read_h5ad(folder / f"mod{number}.h5ad")
                    self.assertEqual(list(out.obs_names), list(cells))
                    self.assertEqual(list(out.var_names), list(genes))
                    expected = raw.layers[layer][rows, :][:, cols]
                    if sp.issparse(expected):
                        expected = expected.toarray()
                    np.testing.assert_array_equal(out.layers["counts"].toarray(), expected)
                    self.assertEqual(out.layers["counts"].dtype, np.dtype("int32"))
                    pd.testing.assert_frame_equal(
                        out.obs[["sample", "day"]].astype(str),
                        filtered.obs.loc[cells, ["sample", "day"]].astype(str),
                    )
                    self.assertEqual(
                        list(out.obs["is_train"]),
                        ["train" if c.startswith("train_") else "validation" for c in cells],
                    )
                    for cell in cells:
                        source = filtered.obs.loc[cell]
                        annotation = annotations.loc[
                            ("PZ-MB-TChIC-" + str(source["sample"]),
                             "H3K4me3", str(source["barcode_nla"]))
                        ]
                        for target, origin in [
                            ("cell_type", "celltype"),
                            ("celltype_general", "celltype_general"),
                            ("lineage", "lineage"),
                        ]:
                            self.assertEqual(out.obs.loc[cell, target], annotation[origin])


if __name__ == "__main__":
    unittest.main()
