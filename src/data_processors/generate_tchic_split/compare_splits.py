import argparse
import sys

import pandas as pd


KEY_COLUMNS = ["mark", "sample"]
VALUE_COLUMN = "split"


def read_mapping(path):
    frame = pd.read_csv(path, dtype={"mark": str, "sample": str, "split": str})
    missing = set(KEY_COLUMNS + [VALUE_COLUMN]) - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    if frame.duplicated(KEY_COLUMNS).any():
        raise ValueError(f"{path} contains duplicate (mark, sample) keys.")
    return frame.set_index(KEY_COLUMNS)[VALUE_COLUMN].sort_index()


def compare_mappings(reference_path, generated_path):
    reference = read_mapping(reference_path)
    generated = read_mapping(generated_path)
    combined = pd.concat(
        [reference.rename("reference"), generated.rename("generated")],
        axis=1,
    )
    changed = combined.loc[
        combined["reference"].ne(combined["generated"])
        | combined["reference"].isna()
        | combined["generated"].isna()
    ]
    return combined, changed


def main():
    parser = argparse.ArgumentParser(
        description="Compare split assignments by (mark, sample), ignoring CSV row order."
    )
    parser.add_argument("reference_csv")
    parser.add_argument("generated_csv")
    args = parser.parse_args()

    combined, changed = compare_mappings(args.reference_csv, args.generated_csv)
    print(f"Reference plates: {combined['reference'].notna().sum()}")
    print(f"Generated plates: {combined['generated'].notna().sum()}")
    print(f"Different/missing assignments: {len(changed)}")
    if not changed.empty:
        print(changed.to_string())
        return 1
    print("Assignments are identical independent of row order or other CSV columns.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
