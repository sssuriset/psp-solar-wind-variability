from pathlib import Path
import numpy as np
import pandas as pd
from cdasws import CdasWs

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"

manifest_path = base / "metadata" / "rotation_manifest.csv"
outdir = base / "data" / "raw" / "omni"
summary_path = base / "outputs" / "phase2" / "omni_availability_summary.csv"

outdir.mkdir(parents=True, exist_ok=True)
summary_path.parent.mkdir(parents=True, exist_ok=True)

dataset = "OMNI2_H0_MRG1HR"

variables = [
    "V1800",
    "BX_GSE1800",
    "BY_GSE1800",
    "BZ_GSE1800",
    "ABS_B1800",
    "N1800",
]

cdas = CdasWs()

df = pd.read_csv(manifest_path)

status_cols = [
    "magnetogram_status",
    "omni_status",
    "pfss_status",
    "mapping_status",
    "notes",
]

for col in status_cols:
    df[col] = df[col].fillna("").astype(str)

def get_values(data, name):
    arr = data[name].values
    arr = np.asarray(arr, dtype=float)
    return arr

def get_time(data, name):
    arr = data[name]

    if "Epoch_1800" in arr.coords:
        return arr.coords["Epoch_1800"].values

    if "Epoch" in arr.coords:
        return arr.coords["Epoch"].values

    first_dim = arr.dims[0]
    return arr.coords[first_dim].values

def clean_speed(values):
    arr = np.asarray(values, dtype=float)
    arr[(arr <= 0) | (arr > 3000)] = np.nan
    return arr

def clean_b(values):
    arr = np.asarray(values, dtype=float)
    arr[np.abs(arr) > 1000] = np.nan
    return arr

def clean_density(values):
    arr = np.asarray(values, dtype=float)
    arr[(arr <= 0) | (arr > 10000)] = np.nan
    return arr

rows = []

for i, row in df.iterrows():
    cr = int(row["cr"])

    start = str(row["omni_start"]) + "T00:00:00Z"
    end = str(row["omni_end"]) + "T23:59:59Z"

    print()
    print(f"Checking OMNI for CR {cr}")
    print(f"  Window: {start} to {end}")

    try:
        status, data = cdas.get_data(dataset, variables, start, end)
    except Exception as exc:
        print("  Status: fail_download_error")
        print(f"  Error: {exc}")

        df.loc[i, "omni_status"] = "fail_download_error"
        rows.append({
            "cr": cr,
            "status": "fail_download_error",
            "rows": 0,
            "speed_valid_fraction": np.nan,
            "imf_valid_fraction": np.nan,
            "notes": str(exc),
        })
        continue

    try:
        time = pd.to_datetime(get_time(data, "V1800"), utc=True, errors="coerce")

        out = pd.DataFrame({
            "time": time,
            "speed_km_s": clean_speed(get_values(data, "V1800")),
            "bx_gse_nt": clean_b(get_values(data, "BX_GSE1800")),
            "by_gse_nt": clean_b(get_values(data, "BY_GSE1800")),
            "bz_gse_nt": clean_b(get_values(data, "BZ_GSE1800")),
            "bmag_nt": clean_b(get_values(data, "ABS_B1800")),
            "density_cm3": clean_density(get_values(data, "N1800")),
        })

    except Exception as exc:
        print("  Status: fail_parse_error")
        print(f"  Error: {exc}")

        df.loc[i, "omni_status"] = "fail_parse_error"
        rows.append({
            "cr": cr,
            "status": "fail_parse_error",
            "rows": 0,
            "speed_valid_fraction": np.nan,
            "imf_valid_fraction": np.nan,
            "notes": str(exc),
        })
        continue

    out = out.dropna(subset=["time"]).sort_values("time")

    csv_path = outdir / f"omni_cr{cr}.csv"
    out.to_csv(csv_path, index=False)

    n_rows = len(out)

    if n_rows == 0:
        speed_valid_fraction = 0.0
        imf_valid_fraction = 0.0
    else:
        speed_valid_fraction = float(out["speed_km_s"].notna().mean())

        imf_valid = (
            out["bx_gse_nt"].notna()
            & out["by_gse_nt"].notna()
            & out["bz_gse_nt"].notna()
        )

        imf_valid_fraction = float(imf_valid.mean())

    start_dt = pd.to_datetime(start)
    end_dt = pd.to_datetime(end)

    expected_hours = int((end_dt - start_dt).total_seconds() / 3600) + 1
    row_fraction = n_rows / expected_hours if expected_hours > 0 else 0

    if row_fraction >= 0.80 and speed_valid_fraction >= 0.80 and imf_valid_fraction >= 0.80:
        omni_status = "pass"
    else:
        omni_status = "backup_or_review_needed"

    df.loc[i, "omni_status"] = omni_status

    old_note = df.loc[i, "notes"]
    extra_note = (
        f"OMNI check {omni_status}: rows={n_rows}, "
        f"speed_valid={speed_valid_fraction:.3f}, "
        f"imf_valid={imf_valid_fraction:.3f}."
    )

    if "OMNI check" not in old_note:
        df.loc[i, "notes"] = old_note + " | " + extra_note

    print(f"  Rows saved: {n_rows}")
    print(f"  Expected hourly rows: {expected_hours}")
    print(f"  Row fraction: {row_fraction:.3f}")
    print(f"  Speed valid fraction: {speed_valid_fraction:.3f}")
    print(f"  IMF valid fraction: {imf_valid_fraction:.3f}")
    print(f"  Saved: {csv_path}")
    print(f"  Status: {omni_status}")

    rows.append({
        "cr": cr,
        "status": omni_status,
        "rows": n_rows,
        "expected_hours": expected_hours,
        "row_fraction": row_fraction,
        "speed_valid_fraction": speed_valid_fraction,
        "imf_valid_fraction": imf_valid_fraction,
        "saved_file": str(csv_path),
    })

df.to_csv(manifest_path, index=False)

summary = pd.DataFrame(rows)
summary.to_csv(summary_path, index=False)

print()
print("OMNI availability summary:")
print(summary.to_string(index=False))

print()
print(f"Updated manifest: {manifest_path}")
print(f"Saved summary: {summary_path}")
print()
print("Status: pass")
