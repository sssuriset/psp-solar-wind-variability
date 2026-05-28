from pathlib import Path
import os
import shutil
import time

import astropy.units as u
import numpy as np
import pandas as pd
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

base = Path(__file__).resolve().parents[1]
outdir = base / "data" / "raw" / "magnetograms" / "hmi"
manifest_path = base / "metadata" / "rotation_manifest.csv"

outdir.mkdir(parents=True, exist_ok=True)

if "JSOC_EMAIL" not in os.environ:
    raise RuntimeError("Set JSOC_EMAIL first with: export JSOC_EMAIL='your registered email'")

df = pd.read_csv(manifest_path)

for col in ["magnetogram_status", "omni_status", "pfss_status", "mapping_status", "notes", "magnetogram_source"]:
    if col not in df.columns:
        df[col] = ""
    df[col] = df[col].fillna("").astype(str)

crs = sorted(df["cr"].dropna().astype(int).unique())

print(f"Project base: {base}")
print(f"JSOC email: {os.environ['JSOC_EMAIL']}")
print(f"Carrington rotations in manifest: {crs}")
print()

series = a.jsoc.Series("hmi.synoptic_mr_polfil_720s")
notify = a.jsoc.Notify(os.environ["JSOC_EMAIL"])

for cr in crs:
    raw_path = outdir / f"hmi_cr{cr}.fits"
    small_path = outdir / f"hmi_cr{cr}_resampled_180x360.fits"

    print("=" * 70)
    print(f"CR {cr}")

    if raw_path.exists() and small_path.exists():
        print("Files already exist. Verifying resampled map only.")
        test_map = sunpy.map.Map(small_path)
        print(f"Existing resampled shape: {test_map.data.shape}")

        mask = df["cr"] == cr
        df.loc[mask, "magnetogram_source"] = "HMI"
        df.loc[mask, "magnetogram_status"] = "pass_hmi_existing"
        continue

    print("Searching JSOC...")
    result = Fido.search(series, a.jsoc.PrimeKey("CAR_ROT", cr), notify)
    print(result)

    if len(result) == 0:
        print(f"No HMI synoptic result found for CR {cr}.")
        mask = df["cr"] == cr
        df.loc[mask, "magnetogram_source"] = "HMI"
        df.loc[mask, "magnetogram_status"] = "fail_hmi_missing"
        continue

    print("Fetching file...")
    files = Fido.fetch(result, path=str(outdir / "{file}"))

    if len(files) == 0:
        print(f"Fido.fetch returned no files for CR {cr}.")
        mask = df["cr"] == cr
        df.loc[mask, "magnetogram_source"] = "HMI"
        df.loc[mask, "magnetogram_status"] = "fail_hmi_fetch"
        continue

    downloaded = Path(files[0])

    if downloaded != raw_path:
        shutil.copy2(downloaded, raw_path)

    print(f"Saved HMI map: {raw_path}")

    hmi_map = sunpy.map.Map(raw_path)

    print(f"Original shape: {hmi_map.data.shape}")
    print(f"Unit: {hmi_map.unit}")
    print(f"Date: {hmi_map.date}")

    finite = np.isfinite(hmi_map.data)
    finite_fraction = finite.sum() / hmi_map.data.size
    print(f"Finite fraction: {finite_fraction:.3f}")

    if finite_fraction < 0.80:
        print("Too many missing values. Marking fail.")
        mask = df["cr"] == cr
        df.loc[mask, "magnetogram_source"] = "HMI"
        df.loc[mask, "magnetogram_status"] = "fail_hmi_missing_values"
        continue

    small_map = hmi_map.resample([360, 180] * u.pix)
    print(f"Resampled shape: {small_map.data.shape}")

    small_map.save(small_path, overwrite=True)
    print(f"Saved resampled map: {small_path}")

    mask = df["cr"] == cr
    df.loc[mask, "magnetogram_source"] = "HMI"
    df.loc[mask, "magnetogram_status"] = "pass_hmi_downloaded"

    old_note = df.loc[mask, "notes"].iloc[0]
    extra_note = "HMI synoptic magnetogram downloaded and resampled."

    if extra_note not in old_note:
        df.loc[mask, "notes"] = old_note + " | " + extra_note

    df.to_csv(manifest_path, index=False)

    time.sleep(2)

df.to_csv(manifest_path, index=False)

print()
print("=" * 70)
print("HMI download summary:")
print(df[["cr", "magnetogram_source", "magnetogram_status"]].to_string(index=False))
print()
print("Status: pass")
