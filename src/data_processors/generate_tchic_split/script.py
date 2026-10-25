import hashlib
import json
import os
import platform
import tempfile
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input_metadata": "resources_test/tchic_split/metadata.h5ad",
    "input_settings": "src/data_processors/generate_tchic_split/h3k4me3_split_settings.json",
    "output_manifest": "outputs/tchic/splits/generated_h3k4me3_plate_split.csv",
    "output_provenance": "outputs/tchic/splits/generated_h3k4me3_plate_split.provenance.json",
}
## VIASH END

PLATE_METADATA_COLUMNS = [
    "plate_batch",
    "day",
    "rep",
    "mark",
    "plate_grp",
    "seqrun",
    "exp_batch",
    "day_rep",
    "batch_day_rep",
]


def load_settings(path):
    with open(path, encoding="utf-8") as settings_file:
        settings = json.load(settings_file)

    required = {
        "mark",
        "partition_names",
        "target_proportions",
        "target_plate_counts",
        "seed",
        "n_iterations",
        "day_order",
        "expected_plates_per_day",
        "fixed_allocations_by_plate_count",
        "special_allocation",
        "required_coverage_columns",
        "score_weights",
    }
    missing = sorted(required - settings.keys())
    if missing:
        raise ValueError(f"Settings are missing required keys: {missing}")

    names = settings["partition_names"]
    proportions = np.asarray(settings["target_proportions"], dtype=float)
    counts = np.asarray(settings["target_plate_counts"], dtype=int)
    if len(names) != 3 or names != ["train", "validation", "test"]:
        raise ValueError(
            "partition_names must be ordered as train, validation, test."
        )
    if len(proportions) != 3 or not np.isclose(proportions.sum(), 1.0):
        raise ValueError("target_proportions must contain three values summing to 1.")
    if len(counts) != 3 or np.any(counts <= 0):
        raise ValueError("target_plate_counts must contain three positive counts.")
    if settings["n_iterations"] <= 0:
        raise ValueError("n_iterations must be positive.")
    if set(settings["score_weights"]) != {
        "within_day_cell_rmse",
        "overall_cell_rmse",
        "day_rep_rmse",
        "rep_rmse",
        "seqrun_rmse",
    }:
        raise ValueError("score_weights must define all five notebook score terms.")

    settings["_target_proportions_array"] = proportions
    settings["_target_plate_counts_array"] = counts
    settings["_expected_plates_per_day_int"] = {
        int(day): int(number)
        for day, number in settings["expected_plates_per_day"].items()
    }
    return settings


def metadata_to_plate_table(obs, mark):
    required_columns = ["sample", *PLATE_METADATA_COLUMNS]
    missing_columns = [column for column in required_columns if column not in obs]
    if missing_columns:
        raise ValueError(f"obs is missing required metadata columns: {missing_columns}")
    if obs.empty:
        raise ValueError("The H5AD obs table contains no cells.")
    if obs[required_columns].isna().any().any():
        missing = obs[required_columns].columns[
            obs[required_columns].isna().any()
        ].tolist()
        raise ValueError(f"Required cell metadata contains missing values: {missing}")
    for column in required_columns:
        empty_values = np.fromiter(
            (
                isinstance(value, str) and not value.strip()
                for value in obs[column].to_numpy()
            ),
            dtype=bool,
            count=len(obs),
        )
        if empty_values.any():
            raise ValueError(f"Required cell metadata contains empty strings: {column}")

    values_per_sample = (
        obs.groupby("sample", observed=True, sort=True)[PLATE_METADATA_COLUMNS]
        .nunique(dropna=False)
    )
    inconsistent = values_per_sample.columns[
        values_per_sample.gt(1).any()
    ].tolist()
    if inconsistent:
        raise ValueError(
            "Plate-level metadata is inconsistent within sample for columns: "
            f"{inconsistent}"
        )

    plate_batches_per_sample = obs.groupby(
        "sample", observed=True, sort=True
    )["plate_batch"].nunique()
    samples_per_plate_batch = obs.groupby(
        "plate_batch", observed=True, sort=True
    )["sample"].nunique()
    if plate_batches_per_sample.gt(1).any():
        raise ValueError("A sample maps to multiple plate_batch values.")
    if samples_per_plate_batch.gt(1).any():
        raise ValueError("A plate_batch maps to multiple sample values.")

    plate_table = (
        obs.groupby("sample", observed=True, sort=True)
        .agg(
            n_cells=("sample", "size"),
            plate_batch=("plate_batch", "first"),
            day=("day", "first"),
            rep=("rep", "first"),
            mark=("mark", "first"),
            plate_grp=("plate_grp", "first"),
            seqrun=("seqrun", "first"),
            exp_batch=("exp_batch", "first"),
            day_rep=("day_rep", "first"),
            batch_day_rep=("batch_day_rep", "first"),
        )
        .reset_index()
    )
    try:
        plate_table["day"] = pd.to_numeric(plate_table["day"].astype(str))
    except (TypeError, ValueError) as error:
        raise ValueError("Metadata column 'day' must contain numeric day labels.") from error

    if not plate_table["sample"].is_unique:
        raise ValueError("Plate metadata contains duplicate sample values.")
    if not plate_table["plate_batch"].is_unique:
        raise ValueError("plate_batch must uniquely identify each sample.")

    selected = plate_table.loc[plate_table["mark"].astype(str) == mark].copy()
    if selected.empty:
        raise ValueError(f"No plates found for configured mark {mark!r}.")
    return selected.reset_index(drop=True)


def validate_split_settings(frame, settings):
    expected = settings["_expected_plates_per_day_int"]
    actual = frame["day"].value_counts().sort_index().to_dict()
    if actual != expected:
        raise ValueError(
            f"Plate counts by day do not match settings: expected {expected}, got {actual}."
        )
    if len(frame) != int(settings["_target_plate_counts_array"].sum()):
        raise ValueError(
            "Configured target plate counts do not sum to the number of selected plates."
        )

    fixed = settings["fixed_allocations_by_plate_count"]
    special = settings["special_allocation"]
    special_count = int(special["plate_count"])
    special_days = [
        day for day, number in expected.items() if number == special_count
    ]
    if len(special_days) < int(special["days_per_candidate"]):
        raise ValueError(
            "There are fewer special-count days than the candidate generator "
            "must select."
        )
    for day, number in expected.items():
        if number != special_count and str(number) not in fixed:
            raise ValueError(
                f"No fixed allocation configured for a {number}-plate day."
            )
    for allocation in [special["default"], *special["choices"]]:
        if len(allocation) != 3 or sum(allocation) != special_count:
            raise ValueError(
                f"Special-day allocation must sum to {special_count}: {allocation}."
            )
    if len(special["choices"]) != 2:
        raise ValueError("The current candidate generator requires two special choices.")
    for count, allocation in fixed.items():
        if len(allocation) != 3 or sum(allocation) != int(count):
            raise ValueError(f"Fixed allocation does not sum to its day size: {count}.")


def generate_candidate_labels(frame, settings, rng):
    labels = np.full(len(frame), fill_value=-1, dtype=np.int8)
    expected = settings["_expected_plates_per_day_int"]
    special = settings["special_allocation"]
    special_count = int(special["plate_count"])
    special_days = [
        day for day, count in expected.items() if count == special_count
    ]
    special_choice = set(
        rng.choice(
            special_days,
            size=int(special["days_per_candidate"]),
            replace=False,
        ).tolist()
    )
    validation_light_day = int(rng.choice(sorted(special_choice)))
    day_values = frame["day"].to_numpy()

    for day in settings["day_order"]:
        day_indices = np.flatnonzero(day_values == day)
        shuffled_indices = rng.permutation(day_indices)
        if len(day_indices) == special_count:
            if day not in special_choice:
                allocation = special["default"]
            elif day == validation_light_day:
                allocation = special["choices"][0]
            else:
                allocation = special["choices"][1]
        else:
            allocation = settings["fixed_allocations_by_plate_count"][
                str(len(day_indices))
            ]

        train_end = int(allocation[0])
        validation_end = train_end + int(allocation[1])
        labels[shuffled_indices[:train_end]] = 0
        labels[shuffled_indices[train_end:validation_end]] = 1
        labels[shuffled_indices[validation_end:]] = 2

    if np.any(labels < 0):
        raise RuntimeError("Candidate generation left one or more plates unassigned.")
    if not np.array_equal(
        np.bincount(labels, minlength=3),
        settings["_target_plate_counts_array"],
    ):
        raise RuntimeError("Generated candidate does not match target plate counts.")
    return labels


def category_indicator_matrix(frame, column):
    values = frame[column].astype(str).to_numpy()
    categories = np.sort(np.unique(values))
    matrix = np.vstack([values == category for category in categories]).astype(np.int8)
    return categories, matrix


def category_counts_by_split(matrix, labels):
    return np.column_stack(
        [
            matrix[:, labels == split_index].sum(axis=1)
            for split_index in range(3)
        ]
    )


def has_required_coverage(frame, labels, settings):
    for column in settings["required_coverage_columns"]:
        if column not in frame:
            raise ValueError(f"Coverage metadata column {column!r} is missing.")
        _, matrix = category_indicator_matrix(frame, column)
        if np.any(category_counts_by_split(matrix, labels) == 0):
            return False
    return True


def prepare_search_data(frame, settings):
    return {
        "cell_counts": frame["n_cells"].to_numpy(dtype=float),
        "day_masks": {
            day: frame["day"].to_numpy() == day
            for day in settings["day_order"]
        },
        "category_matrices": {
            column: category_indicator_matrix(frame, column)[1]
            for column in [
                "day_rep",
                "rep",
                "seqrun",
            ]
        },
    }


def categorical_rmse_from_matrix(matrix, labels, target_proportions):
    counts = category_counts_by_split(matrix, labels)
    totals = matrix.sum(axis=1, keepdims=True)
    observed_proportions = counts / totals
    return float(
        100
        * np.sqrt(
            np.mean((observed_proportions - target_proportions) ** 2)
        )
    )


def score_candidate(labels, prepared, settings):
    target_proportions = settings["_target_proportions_array"]
    cell_counts = prepared["cell_counts"]
    overall_cells_by_split = np.array(
        [
            cell_counts[labels == split_index].sum()
            for split_index in range(3)
        ]
    )
    overall_cell_proportions = overall_cells_by_split / overall_cells_by_split.sum()
    overall_cell_rmse = float(
        100
        * np.sqrt(np.mean((overall_cell_proportions - target_proportions) ** 2))
    )

    within_day_differences = []
    for day in settings["day_order"]:
        day_mask = prepared["day_masks"][day]
        day_labels = labels[day_mask]
        day_cells = cell_counts[day_mask]
        day_plate_counts = np.bincount(day_labels, minlength=3)
        day_plate_proportions = day_plate_counts / day_plate_counts.sum()
        day_cell_counts = np.array(
            [
                day_cells[day_labels == split_index].sum()
                for split_index in range(3)
            ]
        )
        day_cell_proportions = day_cell_counts / day_cell_counts.sum()
        within_day_differences.extend(
            day_cell_proportions - day_plate_proportions
        )

    components = {
        "within_day_cell_rmse": float(
            100 * np.sqrt(np.mean(np.square(within_day_differences)))
        ),
        "overall_cell_rmse": overall_cell_rmse,
        "day_rep_rmse": categorical_rmse_from_matrix(
            prepared["category_matrices"]["day_rep"],
            labels,
            target_proportions,
        ),
        "rep_rmse": categorical_rmse_from_matrix(
            prepared["category_matrices"]["rep"],
            labels,
            target_proportions,
        ),
        "seqrun_rmse": categorical_rmse_from_matrix(
            prepared["category_matrices"]["seqrun"],
            labels,
            target_proportions,
        ),
    }
    components["total_score"] = sum(
        settings["score_weights"][name] * value
        for name, value in components.items()
    )
    return components


def validate_assignment(frame, labels, settings):
    if not np.array_equal(
        np.bincount(labels, minlength=3),
        settings["_target_plate_counts_array"],
    ):
        raise ValueError("Assignment does not match configured plate counts.")
    counts_by_day = {
        day: np.bincount(
            labels[frame["day"].to_numpy() == day],
            minlength=3,
        )
        for day in settings["day_order"]
    }
    special = settings["special_allocation"]
    for day, day_counts in counts_by_day.items():
        if len(frame.loc[frame["day"] == day]) == int(special["plate_count"]):
            allowed = [
                special["default"],
                *special["choices"],
            ]
        else:
            allowed = [
                settings["fixed_allocations_by_plate_count"][
                    str(len(frame.loc[frame["day"] == day]))
                ]
            ]
        if tuple(day_counts.tolist()) not in [tuple(row) for row in allowed]:
            raise ValueError(
                f"Assignment violates the configured allocation for day {day}."
            )
    for column in settings["required_coverage_columns"]:
        _, matrix = category_indicator_matrix(frame, column)
        if np.any(category_counts_by_split(matrix, labels) == 0):
            raise ValueError(f"Assignment violates coverage for {column!r}.")


def search_best_split(frame, settings):
    validate_split_settings(frame, settings)
    rng = np.random.default_rng(int(settings["seed"]))
    prepared = prepare_search_data(frame, settings)
    best_labels = None
    best_scores = None
    best_total_score = np.inf
    valid_candidates = 0

    for _ in range(int(settings["n_iterations"])):
        labels = generate_candidate_labels(frame, settings, rng)
        if not has_required_coverage(frame, labels, settings):
            continue
        valid_candidates += 1
        scores = score_candidate(labels, prepared, settings)
        if scores["total_score"] < best_total_score:
            best_total_score = scores["total_score"]
            best_labels = labels.copy()
            best_scores = scores.copy()

    if best_labels is None:
        raise RuntimeError(
            "No valid candidate was found after "
            f"{settings['n_iterations']} iterations."
        )
    validate_assignment(frame, best_labels, settings)
    result = frame.copy()
    result["split"] = pd.Categorical(
        np.asarray(settings["partition_names"])[best_labels],
        categories=settings["partition_names"],
        ordered=True,
    )
    result = result.sort_values("sample", kind="stable").reset_index(drop=True)
    return result, best_scores, valid_candidates


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def metadata_fingerprint(obs):
    canonical = obs.to_json(orient="split", date_format="iso", double_precision=15)
    return sha256_text(canonical)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_h5ad_obs_metadata(input_metadata):
    required_columns = ["sample", *PLATE_METADATA_COLUMNS]
    with h5py.File(input_metadata, "r") as source:
        if "obs" not in source or "var" not in source:
            raise ValueError("Input H5AD is missing its obs or var group.")

        obs_group = source["obs"]
        missing_columns = [
            column for column in required_columns if column not in obs_group
        ]
        if missing_columns:
            raise ValueError(
                f"H5AD obs is missing required metadata columns: {missing_columns}"
            )

        obs_index_key = obs_group.attrs.get("_index", "_index")
        if isinstance(obs_index_key, bytes):
            obs_index_key = obs_index_key.decode("utf-8")
        if obs_index_key not in obs_group:
            raise ValueError(
                f"H5AD obs index element {obs_index_key!r} is missing."
            )

        obs_index = anndata.io.read_elem(obs_group[obs_index_key])
        obs = pd.DataFrame(
            {
                column: anndata.io.read_elem(obs_group[column])
                for column in required_columns
            },
            index=pd.Index(obs_index, name=obs_index_key),
        )

        var_group = source["var"]
        var_index_key = var_group.attrs.get("_index", "_index")
        if isinstance(var_index_key, bytes):
            var_index_key = var_index_key.decode("utf-8")
        if var_index_key not in var_group:
            raise ValueError(
                f"H5AD var index element {var_index_key!r} is missing."
            )
        source_shape = [len(obs_index), int(var_group[var_index_key].shape[0])]

    return obs, source_shape


def atomic_write_text(path, content):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")
    descriptor, temporary_path = tempfile.mkstemp(
        prefix=f".{output_path.name}.",
        dir=output_path.parent,
        text=True,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output:
            output.write(content)
        os.replace(temporary_path, output_path)
    except Exception:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise


def generate_from_h5ad(input_metadata, settings):
    obs, source_shape = read_h5ad_obs_metadata(input_metadata)
    table = metadata_to_plate_table(obs, settings["mark"])
    result, scores, valid_candidates = search_best_split(table, settings)
    return obs, source_shape, result, scores, valid_candidates


def main():
    output_paths = [
        Path(par["output_manifest"]),
        Path(par["output_provenance"]),
    ]
    existing = [str(path) for path in output_paths if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs: " + ", ".join(existing)
        )
    settings = load_settings(par["input_settings"])
    obs, source_shape, manifest, scores, valid_candidates = generate_from_h5ad(
        par["input_metadata"], settings
    )
    manifest_csv = manifest.to_csv(index=False, lineterminator="\n")
    manifest_hash = sha256_text(manifest_csv)
    provenance = {
        "manifest_test_label": (
            "Final holdout. Later preparation uses validation for benchmark test slots."
        ),
        "input_metadata_path": str(Path(par["input_metadata"]).resolve()),
        "input_shape_cells_features": source_shape,
        "metadata_fingerprint_sha256": metadata_fingerprint(obs),
        "manifest_sha256": manifest_hash,
        "settings": {
            key: value
            for key, value in settings.items()
            if not key.startswith("_")
        },
        "search": {
            "valid_candidates": int(valid_candidates),
            "iterations": int(settings["n_iterations"]),
            "selected_scores": {
                key: float(value) for key, value in scores.items()
            },
            "tie_break": "first candidate with the strictly lowest score",
            "rng": type(np.random.default_rng().bit_generator).__name__,
        },
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "anndata": anndata.__version__,
        },
    }
    atomic_write_text(par["output_manifest"], manifest_csv)
    atomic_write_text(
        par["output_provenance"],
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
    )
    print(f"Wrote split manifest: {par['output_manifest']}")
    print(f"Wrote split provenance: {par['output_provenance']}")
    print(f"Valid candidates: {valid_candidates}")
    print(f"Best total score: {scores['total_score']:.10f}")


if __name__ == "__main__":
    main()
