from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sunpy.map
from sunkit_magex import pfss

base = Path(__file__).resolve().parents[1]

cr = 2287
nrho = 35
rss = 2.5

mag_path = base / "data" / "raw" / "magnetograms" / "hmi" / f"hmi_cr{cr}_resampled_180x360.fits"
manifest_path = base / "metadata" / "rotation_manifest.csv"

outdir = base / "data" / "processed" / "pfss" / "hmi"
figdir = base / "outputs" / "figures" / "pfss" / "hmi"

outdir.mkdir(parents=True, exist_ok=True)
figdir.mkdir(parents=True, exist_ok=True)

print(f"Project base: {base}")
print(f"Running HMI PFSS debug for CR {cr}")
print(f"Magnetogram: {mag_path}")

if not mag_path.exists():
    raise FileNotFoundError(f"Missing magnetogram: {mag_path}")

hmi_map = sunpy.map.Map(mag_path)

print()
print("Input map:")
print(f"Shape: {hmi_map.data.shape}")
print(f"Unit: {hmi_map.unit}")
print(f"Date: {hmi_map.date}")
print(f"Min Br: {np.nanmin(hmi_map.data):.3f}")
print(f"Max Br: {np.nanmax(hmi_map.data):.3f}")
print(f"Finite fraction: {np.isfinite(hmi_map.data).sum() / hmi_map.data.size:.3f}")

print()
print(f"Building PFSS input with nrho={nrho}, rss={rss}...")
pfss_in = pfss.Input(hmi_map, nrho, rss)

print("Solving PFSS...")
pfss_out = pfss.pfss(pfss_in)

ss_br = pfss_out.source_surface_br

print()
print("Source surface map:")
print(f"Shape: {ss_br.data.shape}")
print(f"Unit: {ss_br.unit}")
print(f"Min Br: {np.nanmin(ss_br.data):.6f}")
print(f"Max Br: {np.nanmax(ss_br.data):.6f}")
print(f"Mean unsigned Br: {np.nanmean(np.abs(ss_br.data)):.6f}")
print(f"Finite fraction: {np.isfinite(ss_br.data).sum() / ss_br.data.size:.3f}")

source_surface_path = outdir / f"pfss_hmi_cr{cr}_rss{rss:.1f}_source_surface_br.fits"
ss_br.save(source_surface_path, overwrite=True)
print(f"Saved source surface Br map: {source_surface_path}")

metrics_path = outdir / f"pfss_hmi_cr{cr}_rss{rss:.1f}_metrics.csv"

metrics = pd.DataFrame([{
    "cr": cr,
    "rss": rss,
    "nrho": nrho,
    "input_shape": str(hmi_map.data.shape),
    "source_surface_shape": str(ss_br.data.shape),
    "input_min_br": float(np.nanmin(hmi_map.data)),
    "input_max_br": float(np.nanmax(hmi_map.data)),
    "source_surface_min_br": float(np.nanmin(ss_br.data)),
    "source_surface_max_br": float(np.nanmax(ss_br.data)),
    "source_surface_mean_abs_br": float(np.nanmean(np.abs(ss_br.data))),
    "source_surface_finite_fraction": float(np.isfinite(ss_br.data).sum() / ss_br.data.size),
}])

metrics.to_csv(metrics_path, index=False)
print(f"Saved metrics: {metrics_path}")

fig_path = figdir / f"pfss_hmi_cr{cr}_rss{rss:.1f}_source_surface_br.png"

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(projection=ss_br)
ss_br.plot(axes=ax)
ax.set_title(f"HMI PFSS Source Surface Br, CR {cr}, rss={rss} Rs")
plt.colorbar()
plt.tight_layout()
fig.savefig(fig_path, dpi=200)
plt.close(fig)

print(f"Saved figure: {fig_path}")

df = pd.read_csv(manifest_path)

for col in ["magnetogram_status", "omni_status", "pfss_status", "mapping_status", "notes", "magnetogram_source"]:
    if col not in df.columns:
        df[col] = ""
    df[col] = df[col].fillna("").astype(str)

mask = df["cr"] == cr
df.loc[mask, "pfss_status"] = "pass_hmi_pfss_debug"

old_note = df.loc[mask, "notes"].iloc[0]
extra_note = f"HMI PFSS debug run passed for rss={rss}, nrho={nrho}."

if extra_note not in old_note:
    df.loc[mask, "notes"] = old_note + " | " + extra_note

df.to_csv(manifest_path, index=False)

print()
print("Updated manifest:")
print(df[df["cr"] == cr][["cr", "magnetogram_status", "pfss_status", "notes"]].to_string(index=False))

print()
print("Status: pass")
