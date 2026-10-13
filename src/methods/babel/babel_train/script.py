import logging
import pickle
import sys

import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csc_matrix, csr_matrix, issparse
import skorch
import torch
from torch import nn
from torch.utils.data import Dataset

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

## VIASH START
par = {
    "input_train_mod1": "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad",
    "input_train_mod2": "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad",
    "input_test_mod1": "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod1.h5ad",
    "output": "output_model.pkl",
    "hidden_dim": 16,
    "nn_epochs": 100,
    "lr": 0.001,
    "batch_size": 64,
    "loss2_weight": 3.0,
    "early_stopping_patience": 20,
}
meta = {"name": "babel", "resources_dir": "src/methods/babel"}
## VIASH END

sys.path.append(meta["resources_dir"])

from chrom_utils import parse_chrom_groups
from exit_codes import exit_non_applicable
from model import AssymSplicedAutoEncoder
from losses import QuadLoss


def _as_csr(X):
    return X.tocsr() if issparse(X) else csr_matrix(X)


def _rna_matrix(adata):
    # Keep everything sparse; the full ATAC/RNA matrices are only densified one cell
    # at a time in PairedDataset.__getitem__ (mirrors the original BABEL data loader).
    # Densifying up front OOMs on real ATAC data (e.g. 130k cells x 229k peaks would
    # need ~222 GiB dense).
    counts = _as_csr(adata.layers.get("counts", adata.X)).astype(np.float32)
    if "normalized" in adata.layers:
        X = _as_csr(adata.layers["normalized"]).astype(np.float32)
    else:
        tmp = ad.AnnData(counts.copy())
        sc.pp.normalize_total(tmp)
        sc.pp.log1p(tmp)
        X = _as_csr(tmp.X).astype(np.float32)
    size_factors = np.asarray(counts.sum(axis=1)).reshape(-1, 1)
    med = np.median(size_factors[size_factors > 0]) if np.any(size_factors > 0) else 1.0
    size_factors = size_factors / (med if med > 0 else 1.0)
    return X, counts, size_factors.astype(np.float32)


def _atac_binarized(adata):
    counts = _as_csr(adata.layers.get("counts", adata.X))
    return (counts > 0).astype(np.float32)


def _row_tensor(mat, idx):
    """Densify a single row of a sparse matrix into a 1-D float32 tensor
    (mirrors BABEL's per-cell `utils.ensure_arr(X[i]).flatten()`)."""
    return torch.from_numpy(mat[idx].toarray().ravel().astype(np.float32))


class PairedDataset(Dataset):
    def __init__(self, x1, x2_per_chrom, y1, y2, size_factors1):
        # Store the RNA/ATAC matrices sparse (CSR); rows are densified on demand in
        # __getitem__ so the full dense matrix is never materialized.
        self.x1 = _as_csr(x1)
        self.x2_per_chrom = [_as_csr(c) for c in x2_per_chrom]
        self.y1 = _as_csr(y1)
        self.y2 = _as_csr(y2)
        self.size_factors1 = torch.from_numpy(size_factors1)

    def __len__(self):
        return self.x1.shape[0]

    def __getitem__(self, idx):
        X = {
            "x1": _row_tensor(self.x1, idx),
            "x2_per_chrom": [_row_tensor(c, idx) for c in self.x2_per_chrom],
            "size_factors1": self.size_factors1[idx],
        }
        y = torch.cat([_row_tensor(self.y1, idx), _row_tensor(self.y2, idx)])  # dummy combined target, unused directly
        return X, y


def _collate(batch):
    xs, ys = zip(*batch)
    x1 = torch.stack([x["x1"] for x in xs])
    size_factors1 = torch.stack([x["size_factors1"] for x in xs])
    n_chroms = len(xs[0]["x2_per_chrom"])
    x2_per_chrom = [torch.stack([x["x2_per_chrom"][i] for x in xs]) for i in range(n_chroms)]
    y = torch.stack(ys)
    return {"x1": x1, "x2_per_chrom": x2_per_chrom, "size_factors1": size_factors1}, y


class BabelModule(nn.Module):
    """Wraps AssymSplicedAutoEncoder so skorch's forward(**X) call signature works
    with the dict-of-tensors batch produced by PairedDataset/_collate."""

    def __init__(self, input_dim1, input_dim2, hidden_dim=16):
        super().__init__()
        self.autoencoder = AssymSplicedAutoEncoder(input_dim1, input_dim2, hidden_dim=hidden_dim)

    def forward(self, x1, x2_per_chrom, size_factors1):
        return self.autoencoder(x1, x2_per_chrom, size_factors1=size_factors1)


class BabelNet(skorch.NeuralNet):
    def __init__(self, *args, n_genes, n_peaks, **kwargs):
        self._n_genes = n_genes
        self._n_peaks = n_peaks
        super().__init__(*args, **kwargs)

    def get_loss(self, y_pred, y_true, X=None, training=False):
        y_true = y_true.to(self.device)
        preds11, preds12, preds21, preds22, _, _ = y_pred
        target1 = y_true[:, : self._n_genes]
        target2_bin = y_true[:, self._n_genes : self._n_genes + self._n_peaks]
        return self.criterion_(preds11, preds12, preds21, preds22, target1, target2_bin)


logger.info("Reading input files...")
adata_mod1_train = ad.read_h5ad(par["input_train_mod1"])
adata_mod2_train = ad.read_h5ad(par["input_train_mod2"])

modality1 = adata_mod1_train.uns.get("modality")
modality2 = adata_mod2_train.uns.get("modality")

if modality1 not in {"ATAC", "H3K4me3"} or modality2 != "GEX":
    exit_non_applicable(f"Unsupported source/target: {modality1!r}, {modality2!r}")
if not adata_mod1_train.obs_names.equals(adata_mod2_train.obs_names):
    raise ValueError("Paired training cell order differs between modalities")
if modality1 == "H3K4me3":
    for data in (adata_mod1_train, adata_mod2_train):
        if "counts" not in data.layers:
            raise ValueError("T-ChIC adaptation requires raw counts layers")
    logger.info("T-ChIC adaptation: binary gene-level H3K4me3, chromosome grouping")

adata_atac, adata_rna = adata_mod1_train, adata_mod2_train
direction = "mod1_is_atac"

logger.info("Preprocessing RNA and ATAC...")
X_rna, Y_rna_counts, size_factors = _rna_matrix(adata_rna)
X_atac_bin = _atac_binarized(adata_atac)

if modality1 == "H3K4me3":
    chrom = adata_atac.var["Chromosome"]
    if chrom.isna().any() or chrom.astype(str).str.strip().eq("").any():
        raise ValueError("Missing chromosome annotation")
    labels = chrom.astype(str).str.strip().to_numpy()
    chrom_groups = {c: np.flatnonzero(labels == c) for c in sorted(set(labels))}
    chrom_counts = [len(indices) for indices in chrom_groups.values()]
else:
    chrom_counts, chrom_groups = parse_chrom_groups(adata_atac.var_names)
chrom_order = np.concatenate(list(chrom_groups.values()))
# Slice peaks per chromosome on a CSC view (fast column indexing); each group stays
# sparse and is densified per-cell later.
X_atac_bin_csc = X_atac_bin.tocsc()
X_atac_per_chrom = [X_atac_bin_csc[:, idxs] for idxs in chrom_groups.values()]

n_genes = X_rna.shape[1]
n_peaks = X_atac_bin.shape[1]

dataset = PairedDataset(X_rna, X_atac_per_chrom, Y_rna_counts, X_atac_bin[:, chrom_order], size_factors)

logger.info("Building model (n_genes=%d, n_peaks=%d, %d chromosome groups)...", n_genes, n_peaks, len(chrom_counts))

# Internal early-stopping split uses training data only.
internal_validation_samples = []
internal_split = skorch.dataset.ValidSplit(0.1)

if modality1 == "H3K4me3":
    metadata = adata_atac.obs[["sample", "day"]].copy()
    if metadata.isna().any().any():
        raise ValueError("Missing sample/day metadata")
    metadata = metadata.astype(str)
    if (metadata.groupby("sample")["day"].nunique() != 1).any():
        raise ValueError("A plate contains multiple developmental days")

    plates = metadata.drop_duplicates().sort_values(["day", "sample"])
    rng = np.random.default_rng(42)
    for day, group in plates.groupby("day", sort=True):
        if len(group) < 2:
            raise ValueError(f"Insufficient training plates for day {day}")
        internal_validation_samples.append(
            str(rng.choice(group["sample"].to_numpy()))
        )

    is_valid = metadata["sample"].isin(internal_validation_samples).to_numpy()
    fit_indices = np.flatnonzero(~is_valid)
    valid_indices = np.flatnonzero(is_valid)

    def internal_split(dataset, y=None, **kwargs):
        return (
            torch.utils.data.Subset(dataset, fit_indices),
            torch.utils.data.Subset(dataset, valid_indices),
        )

    logger.info(
        "Internal plate split: %d fitting cells / %d validation cells; held-out plates=%s",
        len(fit_indices), len(valid_indices), internal_validation_samples,
    )

net = BabelNet(
    module=BabelModule,
    module__input_dim1=n_genes,
    module__input_dim2=chrom_counts,
    module__hidden_dim=par["hidden_dim"],
    n_genes=n_genes,
    n_peaks=n_peaks,
    criterion=QuadLoss,
    criterion__loss2_weight=par["loss2_weight"],
    optimizer=torch.optim.Adam,
    lr=par["lr"],
    max_epochs=par["nn_epochs"],
    batch_size=par["batch_size"],
    iterator_train__collate_fn=_collate,
    iterator_train__shuffle=True,
    iterator_valid__collate_fn=_collate,
    train_split=internal_split,
    callbacks=[
        skorch.callbacks.EarlyStopping(patience=par["early_stopping_patience"], load_best=True),
        skorch.callbacks.LRScheduler(
            policy=torch.optim.lr_scheduler.ReduceLROnPlateau,
            monitor="valid_loss",
            mode="min",
            factor=0.1,
            patience=10,
        ),
        skorch.callbacks.GradientNormClipping(gradient_clip_value=5),
    ],
    device="cuda" if torch.cuda.is_available() else "cpu",
)

logger.info("Training BABEL autoencoder...")
net.fit(dataset, y=None)

logger.info("Saving model bundle...")
bundle = {
    "source_modality": modality1,
    "internal_validation_samples": internal_validation_samples,
    "adaptation": "tchic_genelevel_binary_v1" if modality1 == "H3K4me3" else "none",
    "state_dict": net.module_.autoencoder.state_dict(),
    "n_genes": n_genes,
    "n_peaks": n_peaks,
    "chrom_counts": chrom_counts,
    "chrom_groups": chrom_groups,
    "hidden_dim": par["hidden_dim"],
    "direction": direction,
    "rna_var_names": list(adata_rna.var_names),
    "atac_var_names": list(adata_atac.var_names),
    "size_factor_median": float(np.median(np.asarray(Y_rna_counts.sum(axis=1)).ravel())),
}

with open(par["output"], "wb") as f:
    pickle.dump(bundle, f, protocol=4)

logger.info("Training complete. Model saved to %s", par["output"])
