cat("Loading dependencies\n")
if (file.exists("/opt/tchic-python/bin/python")) {
  Sys.setenv(RETICULATE_PYTHON = "/opt/tchic-python/bin/python",
             RETICULATE_USE_MANAGED_VENV = "no")
}
library(anndata, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod1.h5ad",
  input_mod2 = "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod2.h5ad",
  output_train_mod1 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_cite/train_mod1.h5ad",
  output_train_mod2 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_cite/train_mod2.h5ad",
  output_test_mod1 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_cite/test_mod1.h5ad",
  output_test_mod2 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_cite/test_mod2.h5ad",
  swap = TRUE,
  seed = 1L,
  max_test_cells = 1000L,
  preserve_obs_columns = "batch,size_factors",
  preserve_var_columns = "gene_ids,hvg,hvg_score"
)
# par <- list(
#   input_mod1 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/output_mod1.h5ad",
#   input_mod2 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/output_atac.h5ad",
#   output_train_mod1 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/train_mod1.h5ad",
#   output_train_mod2 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/train_mod2.h5ad",
#   output_test_mod1 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/test_mod1.h5ad",
#   output_test_mod2 = "resources_test/task_predict_modality/openproblems_neurips2021/bmmc_multiome/test_mod2.h5ad",
#   swap = TRUE,
#   seed = 1L
# )
## VIASH END

cat("Using seed ", par$seed, "\n", sep = "")
set.seed(par$seed)

cat("Reading input data\n")
ad1 <- anndata::read_h5ad(if (!par$swap) par$input_mod1 else par$input_mod2)
ad2 <- anndata::read_h5ad(if (!par$swap) par$input_mod2 else par$input_mod1)

cell_ids <- function(ad, label) {
  ids <- rownames(ad$obs)
  if (is.null(ids) || length(ids) != nrow(ad$obs)) {
    stop(label, " observation IDs could not be read")
  }
  ids <- as.character(ids)
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(label, " observation IDs must be non-missing and unique")
  }
  ids
}
ad1_cell_ids <- cell_ids(ad1, "input_mod1")
ad2_cell_ids <- cell_ids(ad2, "input_mod2")
if (!identical(ad1_cell_ids, ad2_cell_ids)) {
  stop("Input modalities must have exactly aligned cell IDs in the same order")
}
if (!all(c("counts", "normalized") %in% names(ad1$layers)) ||
    !all(c("counts", "normalized") %in% names(ad2$layers))) {
  stop("Both input modalities must contain counts and normalized layers")
}
if (("is_train" %in% colnames(ad1$obs)) != ("is_train" %in% colnames(ad2$obs))) {
  stop("is_train must be present in both modalities or neither")
}
if ("is_train" %in% colnames(ad1$obs) &&
    !identical(as.character(ad1$obs[["is_train"]]), as.character(ad2$obs[["is_train"]]))) {
  stop("is_train labels must be aligned between modalities")
}

# input checks -- used to balance the subsample of test cells further down
if (!"cell_type" %in% colnames(ad1$obs)) {
  stop("obs['cell_type'] is required but missing from the input dataset")
}

# use heuristic to determine modality
# TODO: should be removed once modality is stored in the uns
determine_modality <- function(ad, mod1 = TRUE) {
  if ("modality" %in% names(ad$uns)) {
    ad$uns[["modality"]]
  } else if ("feature_types" %in% colnames(ad$var)) {
    unique(ad$var[["feature_types"]])
  } else if (mod1) {
    "GEX"
  } else if (grepl("cite", ad$uns[["dataset_id"]])) {
    "ADT"
  } else if (grepl("multiome", ad$uns[["dataset_id"]])) {
    "ATAC"
  } else {
    stop("Could not determine modality")
  }
}
ad1_mod <- determine_modality(ad1, !par$swap)
ad2_mod <- determine_modality(ad2, par$swap)

# determine new uns
uns_vars <- c("dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism", "normalization_id")
ad1_uns <- ad1$uns[uns_vars]
ad2_uns <- ad2$uns[uns_vars]
ad1_uns$modality <- ad1_mod
ad2_uns$modality <- ad2_mod

# Create new dataset id and name depending on the modality
if (!is.null(par$dataset_id)) {
  ad1_uns[["common_dataset_id"]] <- ad2_uns[["common_dataset_id"]] <- ad1_uns$dataset_id
  ad1_uns$dataset_id <- ad2_uns$dataset_id <- par$dataset_id
}

new_dataset_name <- paste0(ad1_uns$dataset_name, " (", ad1_mod, "2", ad2_mod, ")")
ad1_uns$dataset_name <- ad2_uns$dataset_name <- new_dataset_name

# determine new obsm
ad1_obsm <- ad2_obsm <- list()

# determine new var and obs metadata
parse_columns <- function(x) {
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  unique(values[nzchar(values)])
}
preserve_var_columns <- parse_columns(par$preserve_var_columns)
preserve_obs_columns <- parse_columns(par$preserve_obs_columns)
optional_var_columns <- c("gene_ids", "hvg", "hvg_score")
missing_var_columns <- setdiff(setdiff(preserve_var_columns, optional_var_columns), colnames(ad1$var))
missing_var_columns2 <- setdiff(setdiff(preserve_var_columns, optional_var_columns), colnames(ad2$var))
if (length(missing_var_columns) || length(missing_var_columns2)) {
  stop("Requested preserve_var_columns are missing: ", paste(unique(c(missing_var_columns, missing_var_columns2)), collapse = ", "))
}
optional_obs_columns <- c("size_factors")
missing_obs_columns <- setdiff(setdiff(preserve_obs_columns, optional_obs_columns), colnames(ad1$obs))
missing_obs_columns2 <- setdiff(setdiff(preserve_obs_columns, optional_obs_columns), colnames(ad2$obs))
if (length(missing_obs_columns) || length(missing_obs_columns2)) {
  stop("Requested preserve_obs_columns are missing: ", paste(unique(c(missing_obs_columns, missing_obs_columns2)), collapse = ", "))
}
ad1_var <- ad1$var[, intersect(preserve_var_columns, colnames(ad1$var)), drop = FALSE]
ad2_var <- ad2$var[, intersect(preserve_var_columns, colnames(ad2$var)), drop = FALSE]

if (ad1_mod == "ATAC" && "gene_activity" %in% names(ad1$obsm)) {
  # copy gene activity in new object
  ad1_uns$gene_activity_var_names <- ad1$uns$gene_activity_var_names
  ad1_obsm$gene_activity <- as(ad1$obsm$gene_activity, "CsparseMatrix")
}

if (ad2_mod == "ATAC") {
  # subset to make the task computationally feasible
  if (ncol(ad2) > 10000) {
    poss_ix <- which(Matrix::colSums(ad2$layers[["normalized"]]) > 0)
    # sample by position -- sample(x, n) treats a length-1 x as seq_len(x)
    sel_ix <- sort(poss_ix[sample.int(length(poss_ix), min(10000, length(poss_ix)))])
    ad2 <- ad2[, sel_ix]$copy()
    ad2_var <- ad2_var[sel_ix, , drop = FALSE]
  }

  if ("gene_activity" %in% names(ad2$obsm)) {
    # copy gene activity in new object
    ad2_uns$gene_activity_var_names <- ad2$uns$gene_activity_var_names
    ad2_obsm$gene_activity <- as(ad2$obsm$gene_activity, "CsparseMatrix")
  }
}

cat("Creating train/test split\n")
if ("is_train" %in% colnames(ad1$obs)) {
  is_train <- which(ad1$obs[["is_train"]] == "train")
  is_test <- which(ad1$obs[["is_train"]] != "train")
} else {
  # no predefined split -- hold out a quarter of the batches, so the test cells
  # come from donors the method has not seen
  batches <- unique(as.character(ad1$obs[["batch"]]))
  if (length(batches) > 1) {
    test_batches <- sample(batches, max(1, floor(length(batches) / 4)))
    cat("No obs['is_train'], holding out batches: ", paste(test_batches, collapse = ", "), "\n", sep = "")
    in_test <- as.character(ad1$obs[["batch"]]) %in% test_batches
  } else {
    cat("No obs['is_train'] and only one batch, holding out a quarter of the cells\n")
    in_test <- seq_len(nrow(ad1)) %in% sample.int(nrow(ad1), max(1, floor(nrow(ad1) / 4)))
  }
  is_train <- which(!in_test)
  is_test <- which(in_test)
}
cat("Train cells: ", length(is_train), ", test cells: ", length(is_test), "\n", sep = "")

max_test_cells <- as.integer(par$max_test_cells)
if (is.na(max_test_cells) || max_test_cells == 0L || max_test_cells < -1L) {
  stop("max_test_cells must be -1 (unlimited) or a positive integer")
}
if (max_test_cells > 0L && length(is_test) > max_test_cells) {
  ct <- as.character(ad1$obs[["cell_type"]][is_test])
  ct_tab <- table(ct)
  ct_freq <- setNames(as.vector(ct_tab) / sum(ct_tab), names(ct_tab))
  is_test <- sample(is_test, max_test_cells, prob = sqrt(1 / ct_freq[ct]))
}

obs_columns1 <- intersect(preserve_obs_columns, colnames(ad1$obs))
obs_columns2 <- intersect(preserve_obs_columns, colnames(ad2$obs))
train_obs1 <- ad1$obs[is_train, obs_columns1, drop = FALSE]
test_obs1 <- ad1$obs[is_test, obs_columns1, drop = FALSE]
train_obs2 <- ad2$obs[is_train, obs_columns2, drop = FALSE]
test_obs2 <- ad2$obs[is_test, obs_columns2, drop = FALSE]
subset_mats <- function(li, obs_filt) {
  out <- list()
  for (n in names(li)) {
    out[[n]] <- li[[n]][obs_filt, , drop = FALSE]
  }
  out
}

# Keep layers as Python matrices: reticulate's conversion of integer CSR
# counts to an R dgRMatrix can fail because its x slot requires doubles.
# Slice without densifying or changing count dtype, then pass them directly
# to the existing AnnData constructor. R indices are converted to zero-based.
subset_layers <- function(ad, obs_filt) {
  layers <- ad$layers$.get_py_object()
  selection <- reticulate::tuple(
    reticulate::r_to_py(as.list(as.integer(obs_filt - 1L))),
    reticulate::py_eval("slice(None)", convert = FALSE)
  )
  out <- list()
  for (name in c("counts", "normalized")) {
    out[[name]] <- reticulate::py_get_item(layers$get(name), selection)
  }
  out
}

cat("Create train objects\n")
output_train_mod1 <- anndata::AnnData(
  layers = subset_layers(ad1, is_train),
  obsm = subset_mats(ad1_obsm, is_train),
  obs = train_obs1,
  var = ad1_var,
  uns = ad1_uns
)
output_train_mod2 <- anndata::AnnData(
  layers = subset_layers(ad2, is_train),
  obsm = subset_mats(ad2_obsm, is_train),
  obs = train_obs2,
  var = ad2_var,
  uns = ad2_uns
)

cat("Create test objects\n")
output_test_mod1 <- anndata::AnnData(
  layers = subset_layers(ad1, is_test),
  obsm = subset_mats(ad1_obsm, is_test),
  obs = test_obs1,
  var = ad1_var,
  uns = ad1_uns
)
output_test_mod2 <- anndata::AnnData(
  layers = subset_layers(ad2, is_test),
  obsm = subset_mats(ad2_obsm, is_test),
  obs = test_obs2,
  var = ad2_var,
  uns = ad2_uns
)

cat("Saving output files as h5ad\n")
zzz <- output_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- output_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- output_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- output_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")
