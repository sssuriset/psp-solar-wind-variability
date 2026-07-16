from pathlib import Path
import argparse
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "metadata" / "rotation_manifest.csv"

COLS = [
    "cr",
    "start_date",
    "end_date",
    "omni_start",
    "omni_end",
    "role",
    "magnetogram_source",
    "magnetogram_status",
    "omni_status",
    "pfss_status",
    "mapping_status",
    "notes",
    "proxy_status",
    "omni_feature_status",
    "phase_match_status",
    "ballistic_match_status",
    "lag_scan_status",
    "null_test_status",
]

ROWS = [
    {
        "cr": 2281,
        "start_date": "2024-02-14",
        "end_date": "2024-03-12",
        "omni_start": "2024-02-12",
        "omni_end": "2024-03-19",
        "role": "ambient",
        "magnetogram_source": "HMI",
        "notes": "Baseline early 2024 sample",
    },
    {
        "cr": 2283,
        "start_date": "2024-04-09",
        "end_date": "2024-05-06",
        "omni_start": "2024-04-07",
        "omni_end": "2024-05-13",
        "role": "ambient",
        "magnetogram_source": "HMI",
        "notes": "Pre-storm or transition sample",
    },
    {
        "cr": 2284,
        "start_date": "2024-05-06",
        "end_date": "2024-06-02",
        "omni_start": "2024-05-04",
        "omni_end": "2024-06-09",
        "role": "stress_test",
        "magnetogram_source": "HMI",
        "notes": "Disturbed May 2024 interval. Interpret separately from ambient sample.",
    },
    {
        "cr": 2286,
        "start_date": "2024-06-29",
        "end_date": "2024-07-27",
        "omni_start": "2024-06-27",
        "omni_end": "2024-08-03",
        "role": "ambient",
        "magnetogram_source": "HMI",
        "notes": "Post-active comparison sample",
    },
    {
        "cr": 2287,
        "start_date": "2024-07-27",
        "end_date": "2024-08-23",
        "omni_start": "2024-07-25",
        "omni_end": "2024-08-30",
        "role": "debug",
        "magnetogram_source": "HMI",
        "notes": "Pipeline test rotation. Keep outside primary final sample unless promoted as a sensitivity case.",
    },
    {
        "cr": 2290,
        "start_date": "2024-10-16",
        "end_date": "2024-11-13",
        "omni_start": "2024-10-14",
        "omni_end": "2024-11-20",
        "role": "ambient",
        "magnetogram_source": "HMI",
        "notes": "Near-maximum sample",
    },
]


def build_manifest():
    rows = []
    for row in ROWS:
        full = {col: "" for col in COLS}
        full.update(row)
        rows.append(full)
    return pd.DataFrame(rows, columns=COLS)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if OUT.exists() and not args.force:
        print(f"Manifest exists: {OUT}")
        print("Use --force to overwrite it.")
        return

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df = build_manifest()
    df.to_csv(OUT, index=False)
    print(f"Wrote manifest: {OUT}")
    print(df[["cr", "start_date", "end_date", "omni_start", "omni_end", "role", "magnetogram_source"]])


if __name__ == "__main__":
    main()
