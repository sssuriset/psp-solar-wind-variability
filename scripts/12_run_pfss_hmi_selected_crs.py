from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sunpy.map
from sunkit_magex import pfss

base = Path(__file__).resolve().parents[1]

nrho = 35
rss = 2.5

magdir = base / "data" / "raw" / "magnetograms" / "hmi"
manifest_path = base / "metadata" / "rotation_manifest.csv"

outdir = base / "data" / "processed" / "pfss" / "hmi"
figdir = base / "outputs" / "figures" / "pfss" / "hmi"

outdir.mkdir(parents=True, exist_ok=True)
figdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(manifest_path)

for col in ["magnetogram_status", "omni_status", "pfss_status", "mapping_status", "notes", "magnetogram_source"]:
    if col not in df.columns:
        df[col] = ""
    df[col] = df[col].fillna("").astype(str)

crs = sorted(df["cr"].dropna().astype(int).unique())

all_metrics = []

print(f"Project base: {base}")
print(f"Running HMI PFSS batch with rss={rss}, nrho={nrho}")
print(f"Carrington rotations: {crs}")
print()

for cr in crs:
    print("=" * 70)
    print(f"CR {cr}")

    mag_path = magdir / f"hmi_cr{cr}_resampled_180x360.fits"

    source_surface_path = outdir / f"pfss_hmi_cr{cr}_rss{rss:.1f}_source_surface_br.fits"
    metrics_path = outdir / f"pfss_hmi_cr{cr}_rss{rss:.1f}_metrics.csv"
    fig_path = figdir / f"pfss_hmi_cr{cr}_rss{rss:.1f}_source_surface_br.png"

    if not mag_path.exists():
        print(f"Missing resampled HMI magnetogram: {mag_path}")
        mask = df["cr"] == cr
        df.loc[mask, "pfss_status"] = "fail_missing_hmi_resampled"
        continue

    if source_surface_path.exists() and metrics_path.exists() and fig_path.exists():
        print("PFSS outputs already exist. Reading metrics.")
        existing = pd.read_csv(metrics_path)
        all_metrics.append(existing)
        mask = df["cr"] == cr
        df.loc[mask, "pfss_status"] = "pass_hmi_pfss_existing"
        continue

    hmi_map = sunpy.map.Map(mag_path)

    print(f"Input shape: {hmi_map.data.shape}")
    print(f"Input unit: {hmi_map.unit}")
    print(f"Input date: {hmi_map.date}")
    print(f"Input min Br: {np.nanmin(hmi_map.data):.3f}")
    print(f"Input max Br: {np.nanmax(hmi_map.data):.3f}")
    print(f"Input finite fraction: {np.isfinite(hmi_map.data).sum() / hmi_map.data.size:.3f}")

    print("Building PFSS input...")
    pfss_in = pfss.Input(hmi_map, nrho, rss)

    print("Solving PFSS...")
    pfss_out = pfss.pfss(pfss_in)

    ss_br = pfss_out.source_surface_br

    print(f"Source surface shape: {ss_br.data.shape}")
    print(f"Source surface unit: {ss_br.unit}")
    print(f"Source surface min Br: {np.nanmin(ss_br.data):.6f}")
    print(f"Source surface max Br: {np.nanmax(ss_br.data):.6f}")
    print(f"Source surface mean unsigned Br: {np.nanmean(np.abs(ss_br.data)):.6f}")
    print(f"Source surface finite fraction: {np.isfinite(ss_br.data).sum() / ss_br.data.size:.3f}")

    ss_br.save(source_surface_path, overwrite=True)
    print(f"Saved source surface Br map: {source_surface_path}")

    metrics = pd.DataFrame([{
        "cr": cr,
        "rss": rss,
        "nrho": nrho,
        "input_shape": str(hmi_map.data.shape),
        "source_surface_shape": str(ss_br.data.shape),
        "input_min_br": float(np.nanmin(hmi_map.data)),
        "input_max_br": float(np.nanmax(hmi_map.data)),
        "input_mean_abs_br": float(np.nanmean(np.abs(hmi_map.data))),
        "input_finite_fraction": float(np.isfinite(hmi_map.data).sum() / hmi_map.data.size),
        "source_surface_min_br": float(np.nanmin(ss_br.data)),
        "source_surface_max_br": float(np.nanmax(ss_br.data)),
        "source_surface_mean_abs_br": float(np.nanmean(np.abs(ss_br.data))),
        "source_surface_finite_fraction": float(np.isfinite(ss_br.data).sum() / ss_br.data.size),
    }])

    metrics.to_csv(metrics_path, index=False)
    all_metrics.append(metrics)
    print(f"Saved metrics: {metrics_path}")

    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(projection=ss_br)
    ss_br.plot(axes=ax)
    ax.set_title(f"HMI PFSS Source Surface Br, CR {cr}, rss={rss} Rs")
    plt.colorbar()
    plt.tight_layout()
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    print(f"Saved figure: {fig_path}")

    mask = df["cr"] == cr
    df.loc[mask, "pfss_status"] = "pass_hmi_pfss"

    old_note = df.loc[mask, "notes"].iloc[0]
    extra_note = f"HMI PFSS run passed for rss={rss}, nrho={nrho}."

    if extra_note not in old_note:
        df.loc[mask, "notes"] = old_note + " | " + extra_note

    df.to_csv(manifest_path, index=False)

if all_metrics:
    summary = pd.concat(all_metrics, ignore_index=True)
    summary_path = outdir / f"pfss_hmi_selected_crs_rss{rss:.1f}_metrics_summary.csv"
    summary.to_csv(summary_path, index=False)
    print()
    print(f"Saved batch metrics summary: {summary_path}")

df.to_csv(manifest_path, index=False)

print()
print("=" * 70)
print("PFSS batch summary:")
print(df[["cr", "magnetogram_status", "pfss_status"]].to_string(index=False))

print()
print("Status: pass")
