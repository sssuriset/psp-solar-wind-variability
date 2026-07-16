from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sunpy.map
from sunkit_magex import pfss


ROOT = Path(__file__).resolve().parents[1]
NRHO = 35
RSS_ALL = [2.0, 2.5, 3.0]

MAGDIR = ROOT / "data" / "raw" / "magnetograms" / "hmi"
MANIFEST = ROOT / "metadata" / "rotation_manifest.csv"
OUTDIR = ROOT / "data" / "processed" / "pfss" / "hmi"
FIGDIR = ROOT / "outputs" / "figures" / "pfss" / "hmi"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rss", type=float, nargs="*", default=RSS_ALL)
    ap.add_argument("--cr", type=int, nargs="*", default=None)
    ap.add_argument("--redo", action="store_true")
    return ap.parse_args()


def rss_name(rss):
    return f"rss{rss:.1f}"


def finite_frac(arr):
    return float(np.isfinite(arr).sum() / arr.size)


def prep_manifest(df):
    cols = [
        "magnetogram_status",
        "omni_status",
        "pfss_status",
        "mapping_status",
        "notes",
        "magnetogram_source",
    ]
    for col in cols:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str)
    return df


def add_note(df, cr, msg):
    mask = df["cr"] == cr
    old = df.loc[mask, "notes"].iloc[0]
    if msg not in old:
        sep = " | " if old else ""
        df.loc[mask, "notes"] = old + sep + msg


def paths(cr, rss):
    tag = rss_name(rss)
    return {
        "mag": MAGDIR / f"hmi_cr{cr}_resampled_180x360.fits",
        "ss": OUTDIR / f"pfss_hmi_cr{cr}_{tag}_source_surface_br.fits",
        "metrics": OUTDIR / f"pfss_hmi_cr{cr}_{tag}_metrics.csv",
        "fig": FIGDIR / f"pfss_hmi_cr{cr}_{tag}_source_surface_br.png",
    }


def run_one(cr, rss, df, redo=False):
    p = paths(cr, rss)
    mask = df["cr"] == cr

    print("=" * 70)
    print(f"CR {cr}, rss={rss:.1f}")

    if not p["mag"].exists():
        print(f"Missing resampled HMI magnetogram: {p['mag']}")
        df.loc[mask, "pfss_status"] = "fail_missing_hmi_resampled"
        return None

    if not redo and p["ss"].exists() and p["metrics"].exists() and p["fig"].exists():
        print("PFSS outputs already exist. Reading metrics.")
        df.loc[mask, "pfss_status"] = "pass_hmi_pfss_existing"
        return pd.read_csv(p["metrics"])

    hmi_map = sunpy.map.Map(p["mag"])

    print(f"Input shape: {hmi_map.data.shape}")
    print(f"Input unit: {hmi_map.unit}")
    print(f"Input date: {hmi_map.date}")
    print(f"Input min Br: {np.nanmin(hmi_map.data):.3f}")
    print(f"Input max Br: {np.nanmax(hmi_map.data):.3f}")
    print(f"Input finite fraction: {finite_frac(hmi_map.data):.3f}")

    print("Building PFSS input.")
    pfss_in = pfss.Input(hmi_map, NRHO, rss)

    print("Solving PFSS.")
    pfss_out = pfss.pfss(pfss_in)
    ss_br = pfss_out.source_surface_br

    print(f"Source surface shape: {ss_br.data.shape}")
    print(f"Source surface unit: {ss_br.unit}")
    print(f"Source surface min Br: {np.nanmin(ss_br.data):.6f}")
    print(f"Source surface max Br: {np.nanmax(ss_br.data):.6f}")
    print(f"Source surface mean unsigned Br: {np.nanmean(np.abs(ss_br.data)):.6f}")
    print(f"Source surface finite fraction: {finite_frac(ss_br.data):.3f}")

    ss_br.save(p["ss"], overwrite=True)
    print(f"Saved source surface Br map: {p['ss']}")

    metrics = pd.DataFrame([{
        "cr": cr,
        "rss": rss,
        "nrho": NRHO,
        "input_shape": str(hmi_map.data.shape),
        "source_surface_shape": str(ss_br.data.shape),
        "input_min_br": float(np.nanmin(hmi_map.data)),
        "input_max_br": float(np.nanmax(hmi_map.data)),
        "input_mean_abs_br": float(np.nanmean(np.abs(hmi_map.data))),
        "input_finite_fraction": finite_frac(hmi_map.data),
        "source_surface_min_br": float(np.nanmin(ss_br.data)),
        "source_surface_max_br": float(np.nanmax(ss_br.data)),
        "source_surface_mean_abs_br": float(np.nanmean(np.abs(ss_br.data))),
        "source_surface_finite_fraction": finite_frac(ss_br.data),
    }])

    metrics.to_csv(p["metrics"], index=False)
    print(f"Saved metrics: {p['metrics']}")

    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(projection=ss_br)
    ss_br.plot(axes=ax)
    ax.set_title(f"HMI PFSS Source Surface Br, CR {cr}, rss={rss:.1f} Rs")
    plt.colorbar()
    plt.tight_layout()
    fig.savefig(p["fig"], dpi=200)
    plt.close(fig)
    print(f"Saved figure: {p['fig']}")

    df.loc[mask, "pfss_status"] = "pass_hmi_pfss"
    add_note(df, cr, f"HMI PFSS run passed for rss={rss:.1f}, nrho={NRHO}.")

    return metrics


def run_batch(rss, crs, df, redo=False):
    print()
    print(f"Running HMI PFSS batch with rss={rss:.1f}, nrho={NRHO}")
    print(f"Carrington rotations: {crs}")
    print()

    all_metrics = []
    for cr in crs:
        metrics = run_one(cr, rss, df, redo=redo)
        if metrics is not None:
            all_metrics.append(metrics)
        df.to_csv(MANIFEST, index=False)

    if all_metrics:
        summary = pd.concat(all_metrics, ignore_index=True)
        summary_path = OUTDIR / f"pfss_hmi_selected_crs_{rss_name(rss)}_metrics_summary.csv"
        summary.to_csv(summary_path, index=False)
        print()
        print(f"Saved batch metrics summary: {summary_path}")


def main():
    args = parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(MANIFEST)
    df = prep_manifest(df)

    if args.cr is None:
        crs = sorted(df["cr"].dropna().astype(int).unique())
    else:
        crs = sorted(args.cr)

    print(f"Project base: {ROOT}")

    for rss in args.rss:
        run_batch(float(rss), crs, df, redo=args.redo)

    df.to_csv(MANIFEST, index=False)

    print()
    print("=" * 70)
    print("PFSS batch summary:")
    print(df[["cr", "magnetogram_status", "pfss_status"]].to_string(index=False))
    print()
    print("Status: pass")


if __name__ == "__main__":
    main()
