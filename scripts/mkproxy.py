from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[1]
RSS_ALL = [2.0, 2.5, 3.0]

MANIFEST = ROOT / "metadata" / "rotation_manifest.csv"
PFSSDIR = ROOT / "data" / "processed" / "pfss" / "hmi"
OUTBASE = ROOT / "data" / "processed" / "pfss_longitude_proxies"
FIGBASE = ROOT / "outputs" / "figures" / "pfss_longitude_proxies"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rss", type=float, nargs="*", default=RSS_ALL)
    ap.add_argument("--cr", type=int, nargs="*", default=None)
    return ap.parse_args()


def rss_name(rss):
    return f"rss{rss:.1f}"


def finite_frac(a):
    return float(np.isfinite(a).sum() / a.size)


def sign0(a):
    return np.where(a > 0, 1, np.where(a < 0, -1, 0))


def grids(shape):
    nlat, nlon = shape
    lat = np.linspace(-90.0, 90.0, nlat)
    lon = np.linspace(0.0, 360.0, nlon, endpoint=False)
    return lat, lon


def read_br(path):
    with fits.open(path) as hdul:
        data = hdul[0].data
        if data is None and len(hdul) > 1:
            data = hdul[1].data
    arr = np.array(data, dtype=float)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D source-surface Br map, got shape {arr.shape}")
    return arr


def band_stats(arr, lat, mask):
    sub = arr[mask, :]
    signed = np.nanmean(sub, axis=0)
    unsigned = np.nanmean(np.abs(sub), axis=0)
    return signed, unsigned


def profile_for(cr, rss):
    tag = rss_name(rss)
    src = PFSSDIR / f"pfss_hmi_cr{cr}_{tag}_source_surface_br.fits"
    if not src.exists():
        raise FileNotFoundError(src)

    arr = read_br(src)
    lat, lon = grids(arr.shape)

    eq = np.abs(lat) <= 20.0
    mid = np.abs(lat) <= 45.0
    all_lat = np.isfinite(lat)

    eq_s, eq_a = band_stats(arr, lat, eq)
    mid_s, mid_a = band_stats(arr, lat, mid)
    glob_s, glob_a = band_stats(arr, lat, all_lat)

    prof = pd.DataFrame({
        "cr": cr,
        "rss": rss,
        "longitude_deg": lon,
        "equator_signed_br": eq_s,
        "equator_abs_br": eq_a,
        "midlat_signed_br": mid_s,
        "midlat_abs_br": mid_a,
        "global_signed_br": glob_s,
        "global_abs_br": glob_a,
        "equator_polarity": sign0(eq_s),
        "global_polarity": sign0(glob_s),
    })

    summ = pd.DataFrame([{
        "cr": cr,
        "rss": rss,
        "nlat": arr.shape[0],
        "nlon": arr.shape[1],
        "finite_fraction": finite_frac(arr),
        "equator_abs_mean": float(np.nanmean(eq_a)),
        "equator_abs_std": float(np.nanstd(eq_a)),
        "equator_abs_min": float(np.nanmin(eq_a)),
        "equator_abs_max": float(np.nanmax(eq_a)),
        "equator_abs_range": float(np.nanmax(eq_a) - np.nanmin(eq_a)),
        "midlat_abs_mean": float(np.nanmean(mid_a)),
        "midlat_abs_std": float(np.nanstd(mid_a)),
        "midlat_abs_min": float(np.nanmin(mid_a)),
        "midlat_abs_max": float(np.nanmax(mid_a)),
        "midlat_abs_range": float(np.nanmax(mid_a) - np.nanmin(mid_a)),
        "global_abs_mean": float(np.nanmean(glob_a)),
        "global_abs_std": float(np.nanstd(glob_a)),
        "global_abs_min": float(np.nanmin(glob_a)),
        "global_abs_max": float(np.nanmax(glob_a)),
        "global_abs_range": float(np.nanmax(glob_a) - np.nanmin(glob_a)),
    }])

    return prof, summ


def update_manifest(crs, rss):
    if not MANIFEST.exists():
        return
    df = pd.read_csv(MANIFEST)
    for col in ["proxy_status", "notes"]:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str)

    for cr in crs:
        mask = df["cr"].astype(int) == int(cr)
        if mask.any():
            df.loc[mask, "proxy_status"] = "pass"
            note = f"Longitude proxies built for rss={rss:.1f}."
            old = df.loc[mask, "notes"].iloc[0]
            if note not in old:
                df.loc[mask, "notes"] = old + (" | " if old else "") + note

    df.to_csv(MANIFEST, index=False)


def plot_profiles(all_prof, rss, figdir):
    fig, ax = plt.subplots(figsize=(10, 5))
    for cr, g in all_prof.groupby("cr"):
        ax.plot(g["longitude_deg"], g["equator_abs_br"], label=f"CR {cr}")
    ax.set_xlabel("Carrington longitude [deg]")
    ax.set_ylabel("Equator |Br|")
    ax.set_title(f"PFSS longitude proxy, rss={rss:.1f} Rs")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    out = figdir / f"pfss_hmi_longitude_proxy_equator_abs_{rss_name(rss)}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved figure: {out}")


def run_rss(rss, crs):
    tag = rss_name(rss)
    outdir = OUTBASE / tag
    figdir = FIGBASE / tag
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    profiles = []
    summaries = []

    print()
    print("=" * 70)
    print(f"Building longitude proxies for rss={rss:.1f}")
    print(f"Carrington rotations: {crs}")

    for cr in crs:
        print(f"CR {cr}")
        prof, summ = profile_for(cr, rss)
        prof_path = outdir / f"pfss_hmi_cr{cr}_{tag}_longitude_profile.csv"
        prof.to_csv(prof_path, index=False)
        print(f"Saved profile: {prof_path}")
        profiles.append(prof)
        summaries.append(summ)

    all_prof = pd.concat(profiles, ignore_index=True)
    all_summ = pd.concat(summaries, ignore_index=True)

    all_prof_path = outdir / f"pfss_hmi_selected_crs_{tag}_longitude_profiles_all.csv"
    summ_path = outdir / f"pfss_hmi_selected_crs_{tag}_longitude_proxy_summary.csv"

    all_prof.to_csv(all_prof_path, index=False)
    all_summ.to_csv(summ_path, index=False)

    print(f"Saved all profiles: {all_prof_path}")
    print(f"Saved summary: {summ_path}")

    plot_profiles(all_prof, rss, figdir)
    update_manifest(crs, rss)


def main():
    args = parse_args()

    if args.cr is None:
        df = pd.read_csv(MANIFEST)
        crs = sorted(df["cr"].dropna().astype(int).unique())
    else:
        crs = sorted(args.cr)

    for rss in args.rss:
        run_rss(float(rss), crs)

    print()
    print("Status: pass")


if __name__ == "__main__":
    main()
