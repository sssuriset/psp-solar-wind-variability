from pathlib import Path
import os
import shutil

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

cr = 2287
series = a.jsoc.Series("hmi.synoptic_mr_polfil_720s")
crot = a.jsoc.PrimeKey("CAR_ROT", cr)
notify = a.jsoc.Notify(os.environ["JSOC_EMAIL"])

print(f"Project base: {base}")
print(f"Searching JSOC HMI synoptic map for CR {cr}...")

result = Fido.search(series, crot, notify)

print(result)

if len(result) == 0:
    raise RuntimeError(f"No HMI synoptic result found for CR {cr}")

print()
print("Fetching file...")
files = Fido.fetch(result, path=str(outdir / "{file}"))

if len(files) == 0:
    raise RuntimeError("Fido.fetch returned no files.")

downloaded = Path(files[0])
target = outdir / f"hmi_cr{cr}.fits"

if downloaded != target:
    shutil.copy2(downloaded, target)

print(f"Saved HMI map: {target}")

print()
print("Loading map...")
hmi_map = sunpy.map.Map(target)

print(f"Original shape: {hmi_map.data.shape}")
print(f"Unit: {hmi_map.unit}")
print(f"Date: {hmi_map.date}")

finite = np.isfinite(hmi_map.data)
finite_fraction = finite.sum() / hmi_map.data.size
print(f"Finite fraction: {finite_fraction:.3f}")

if finite_fraction < 0.80:
    raise RuntimeError("HMI map has too many missing values.")

print()
print("Testing resample to PFSS-friendly shape...")
small_map = hmi_map.resample([360, 180] * u.pix)

print(f"Resampled shape: {small_map.data.shape}")

small_path = outdir / f"hmi_cr{cr}_resampled_180x360.fits"
small_map.save(small_path, overwrite=True)

print(f"Saved resampled map: {small_path}")

df = pd.read_csv(manifest_path)

for col in ["magnetogram_status", "omni_status", "pfss_status", "mapping_status", "notes", "magnetogram_source"]:
    df[col] = df[col].fillna("").astype(str)

mask = df["cr"] == cr
df.loc[mask, "magnetogram_source"] = "HMI"
df.loc[mask, "magnetogram_status"] = "pass_hmi_debug"

old_note = df.loc[mask, "notes"].iloc[0]
extra_note = "HMI synoptic magnetogram downloaded and resampled for debug CR."

if extra_note not in old_note:
    df.loc[mask, "notes"] = old_note + " | " + extra_note

df.to_csv(manifest_path, index=False)

print()
print("Updated manifest for CR 2287.")
print(df[df["cr"] == cr][["cr", "magnetogram_source", "magnetogram_status", "notes"]].to_string(index=False))

print()
print("Status: pass")
