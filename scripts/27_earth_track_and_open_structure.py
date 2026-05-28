from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

b = Path.cwd()

heights = [2.0, 2.5, 3.0]

outtab = b / "outputs" / "tables"
outfig = b / "outputs" / "poster_figures"
outres = b / "outputs" / "results"

outtab.mkdir(parents=True, exist_ok=True)
outfig.mkdir(parents=True, exist_ok=True)
outres.mkdir(parents=True, exist_ok=True)

track_rows = []
structure_rows = []

def zero_crossings(y):
    y = np.array(y, dtype=float, copy=True)
    y[~np.isfinite(y)] = np.nan
    s = np.sign(y)
    ok = np.isfinite(s)
    if ok.sum() < 2:
        return np.nan
    s = s[ok]
    s[s == 0] = np.nan
    s = pd.Series(s).ffill().bfill().to_numpy()
    return int(np.sum(s != np.roll(s, 1)))

def strong_patch_count(y, threshold):
    y = np.array(y, dtype=float, copy=True)
    mask = np.isfinite(y) & (y >= threshold)
    if mask.sum() == 0:
        return 0
    return int(np.sum(mask & ~np.roll(mask, 1)))

for rss in heights:
    h = f"rss{rss:.1f}"

    mf = b / "data" / "processed" / "comparison" / h / "pfss_omni_ballistic_matched_rows.csv"
    pf = b / "data" / "processed" / "pfss" / "hmi_longitude_proxies" / f"pfss_hmi_selected_crs_rss{rss:.1f}_longitude_profiles_all.csv"
    lf = b / "data" / "processed" / "comparison" / h / "pfss_omni_ballistic_lag_summary.csv"

    if not mf.exists():
        raise RuntimeError(f"Missing matched rows: {mf}")
    if not pf.exists():
        raise RuntimeError(f"Missing longitude profiles: {pf}")

    m = pd.read_csv(mf)
    p = pd.read_csv(pf)
    lag = pd.read_csv(lf) if lf.exists() else pd.DataFrame()

    for c in m.columns:
        if c != "time":
            m[c] = pd.to_numeric(m[c], errors="coerce")

    for c in p.columns:
        p[c] = pd.to_numeric(p[c], errors="coerce")

    m["rss"] = rss
    p["rss"] = rss

    if "time" in m.columns:
        m["time_dt"] = pd.to_datetime(m["time"], errors="coerce")
    else:
        m["time_dt"] = pd.NaT

    if "ballistic_phase10_deg" not in m.columns:
        m["ballistic_phase10_deg"] = np.floor(np.mod(m["ballistic_lon_deg"], 360.0) / 10.0) * 10.0

    for cr in sorted(m["cr"].dropna().astype(int).unique()):
        q = m[m["cr"].astype(int) == cr].copy()

        occupied_bins = q["ballistic_phase10_deg"].dropna().nunique()
        total_bins = 36
        coverage = occupied_bins / total_bins

        row = {
            "rss": rss,
            "cr": cr,
            "rows": len(q),
            "occupied_10deg_bins": occupied_bins,
            "coverage_fraction_10deg": coverage,
            "mean_speed_km_s": q["speed_km_s"].mean() if "speed_km_s" in q.columns else np.nan,
            "median_speed_km_s": q["speed_km_s"].median() if "speed_km_s" in q.columns else np.nan,
            "mean_omni_bmag_nt": q["bmag_nt"].mean() if "bmag_nt" in q.columns else np.nan,
            "mean_ballistic_lon_deg": q["ballistic_lon_deg"].mean() if "ballistic_lon_deg" in q.columns else np.nan,
            "min_ballistic_lon_deg": q["ballistic_lon_deg"].min() if "ballistic_lon_deg" in q.columns else np.nan,
            "max_ballistic_lon_deg": q["ballistic_lon_deg"].max() if "ballistic_lon_deg" in q.columns else np.nan,
        }

        if len(lag) and "cr" in lag.columns:
            lq = lag[lag["cr"].astype(int) == cr]
            if len(lq):
                for col in ["lag_mean_days", "lag_min_days", "lag_max_days", "shift_mean_deg"]:
                    if col in lq.columns:
                        row[col] = float(lq[col].iloc[0])

        track_rows.append(row)

    needed = ["cr", "lon_deg", "equator_abs_br", "equator_signed_br", "global_abs_br", "global_signed_br"]
    miss = [c for c in needed if c not in p.columns]
    if miss:
        raise RuntimeError(f"Missing profile columns for rss {rss}: {miss}")

    for cr in sorted(p["cr"].dropna().astype(int).unique()):
        q = p[p["cr"].astype(int) == cr].sort_values("lon_deg").copy()

        eq_abs = q["equator_abs_br"].to_numpy(dtype=float)
        glob_abs = q["global_abs_br"].to_numpy(dtype=float)
        eq_signed = q["equator_signed_br"].to_numpy(dtype=float)
        glob_signed = q["global_signed_br"].to_numpy(dtype=float)

        eq_thresh = np.nanpercentile(eq_abs, 75)
        glob_thresh = np.nanpercentile(glob_abs, 75)

        structure_rows.append({
            "rss": rss,
            "cr": cr,
            "equator_abs_mean": np.nanmean(eq_abs),
            "equator_abs_std": np.nanstd(eq_abs),
            "equator_abs_max": np.nanmax(eq_abs),
            "equator_abs_peak_lon_deg": float(q.iloc[int(np.nanargmax(eq_abs))]["lon_deg"]),
            "equator_strong_patch_count_p75": strong_patch_count(eq_abs, eq_thresh),
            "equator_sector_boundary_count": zero_crossings(eq_signed),
            "global_abs_mean": np.nanmean(glob_abs),
            "global_abs_std": np.nanstd(glob_abs),
            "global_abs_max": np.nanmax(glob_abs),
            "global_abs_peak_lon_deg": float(q.iloc[int(np.nanargmax(glob_abs))]["lon_deg"]),
            "global_strong_patch_count_p75": strong_patch_count(glob_abs, glob_thresh),
            "global_sector_boundary_count": zero_crossings(glob_signed),
        })

track = pd.DataFrame(track_rows)
structure = pd.DataFrame(structure_rows)

track_summary = track.groupby("rss", as_index=False).agg(
    total_rows=("rows", "sum"),
    median_occupied_10deg_bins=("occupied_10deg_bins", "median"),
    mean_coverage_fraction_10deg=("coverage_fraction_10deg", "mean"),
    mean_lag_days=("lag_mean_days", "mean"),
    mean_shift_deg=("shift_mean_deg", "mean"),
)

structure_summary = structure.groupby("rss", as_index=False).agg(
    mean_equator_abs_br=("equator_abs_mean", "mean"),
    mean_equator_sector_boundaries=("equator_sector_boundary_count", "mean"),
    mean_equator_strong_patches=("equator_strong_patch_count_p75", "mean"),
    mean_global_abs_br=("global_abs_mean", "mean"),
    mean_global_sector_boundaries=("global_sector_boundary_count", "mean"),
    mean_global_strong_patches=("global_strong_patch_count_p75", "mean"),
)

track_path = outtab / "pfss_omni_earth_connected_track_summary_by_cr.csv"
track_height_path = outtab / "pfss_omni_earth_connected_track_summary_by_height.csv"
structure_path = outtab / "pfss_open_structure_summary_by_cr.csv"
structure_height_path = outtab / "pfss_open_structure_summary_by_height.csv"

track.to_csv(track_path, index=False)
track_summary.to_csv(track_height_path, index=False)
structure.to_csv(structure_path, index=False)
structure_summary.to_csv(structure_height_path, index=False)

plt.figure(figsize=(8, 5))
for rss in heights:
    q = track[track["rss"] == rss]
    plt.plot(q["cr"], q["coverage_fraction_10deg"], marker="o", linewidth=2, label=f"rss {rss:.1f}")
plt.xlabel("Carrington rotation")
plt.ylabel("Earth-connected 10 deg longitude-bin coverage")
plt.title("Ballistic Earth-connected source-surface track coverage")
plt.legend()
plt.tight_layout()
fig_track = outfig / "fig_earth_connected_track_coverage.png"
plt.savefig(fig_track, dpi=300)
plt.close()

plt.figure(figsize=(8, 5))
for rss in heights:
    q = structure[structure["rss"] == rss]
    plt.plot(q["cr"], q["equator_sector_boundary_count"], marker="o", linewidth=2, label=f"rss {rss:.1f}")
plt.xlabel("Carrington rotation")
plt.ylabel("Equatorial source-surface polarity boundaries")
plt.title("PFSS source-surface sector structure")
plt.legend()
plt.tight_layout()
fig_sector = outfig / "fig_pfss_sector_structure_by_cr.png"
plt.savefig(fig_sector, dpi=300)
plt.close()

txt = []
txt.append("Earth-connected source-surface track and open-structure summary")
txt.append("")
txt.append("Earth-connected track definition:")
txt.append("Each OMNI sample was mapped back to a PFSS source-surface longitude using the speed-dependent ballistic lag already used in the correlation analysis.")
txt.append("")
txt.append("Track summary by source-surface height:")
txt.append(track_summary.to_string(index=False))
txt.append("")
txt.append("Open/source-surface structure definition:")
txt.append("Structure was summarized from longitude profiles of source-surface Br using unsigned-field strength, strong-field patch counts, and signed-Br sector-boundary counts.")
txt.append("")
txt.append("Structure summary by source-surface height:")
txt.append(structure_summary.to_string(index=False))
txt.append("")
txt.append("Track summary by CR:")
txt.append(track.to_string(index=False))
txt.append("")
txt.append("Structure summary by CR:")
txt.append(structure.to_string(index=False))

txt_path = outres / "pfss_omni_earth_track_open_structure_summary.txt"
txt_path.write_text("\n".join(txt))

print("Earth-connected track summary by height:")
print(track_summary.to_string(index=False))
print()
print("Open/source-surface structure summary by height:")
print(structure_summary.to_string(index=False))
print()
print("Saved tables:")
print(track_path)
print(track_height_path)
print(structure_path)
print(structure_height_path)
print()
print("Saved figures:")
print(fig_track)
print(fig_sector)
print()
print("Saved text summary:")
print(txt_path)
print()
print("Status: pass")
