from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits

b = Path.cwd()

heights = [2.0, 2.5, 3.0]

outtab = b / "outputs" / "tables"
outfig = b / "outputs" / "figures" / "diag"
outres = b / "outputs" / "results"

outtab.mkdir(parents=True, exist_ok=True)
outfig.mkdir(parents=True, exist_ok=True)
outres.mkdir(parents=True, exist_ok=True)

def read_2d_fits(path):
    with fits.open(path) as hdul:
        for hdu in hdul:
            data = hdu.data
            if data is None:
                continue
            arr = np.asarray(data)
            arr = np.squeeze(arr)
            if arr.ndim == 2:
                arr = arr.astype(float)
                arr[~np.isfinite(arr)] = np.nan
                return arr
    raise RuntimeError(f"No 2D image data found in {path}")

def interp_to_360(lon, values):
    lon = np.asarray(lon, dtype=float)
    values = np.asarray(values, dtype=float)

    ok = np.isfinite(lon) & np.isfinite(values)
    if ok.sum() < 2:
        return np.full(360, np.nan)

    x = np.mod(lon[ok], 360.0)
    y = values[ok]

    s = np.argsort(x)
    x = x[s]
    y = y[s]

    xx = np.concatenate([x - 360.0, x, x + 360.0])
    yy = np.concatenate([y, y, y])

    return np.interp(np.arange(360, dtype=float), xx, yy)

def magnetogram_profile(arr):
    nlat, nlon = arr.shape

    lon = np.linspace(0.0, 360.0, nlon, endpoint=False)
    lat = np.linspace(-90.0 + 90.0 / nlat, 90.0 - 90.0 / nlat, nlat)

    eq = np.abs(lat) <= 5.0
    mid = (np.abs(lat) >= 20.0) & (np.abs(lat) <= 40.0)

    eq_abs = np.nanmean(np.abs(arr[eq, :]), axis=0)
    mid_abs = np.nanmean(np.abs(arr[mid, :]), axis=0)
    glob_abs = np.nanmean(np.abs(arr), axis=0)

    eq_signed = np.nanmean(arr[eq, :], axis=0)
    mid_signed = np.nanmean(arr[mid, :], axis=0)
    glob_signed = np.nanmean(arr, axis=0)

    return pd.DataFrame({
        "lon_bin": np.arange(360, dtype=int),
        "photo_equator_abs_br": interp_to_360(lon, eq_abs),
        "photo_midlat_abs_br": interp_to_360(lon, mid_abs),
        "photo_global_abs_br": interp_to_360(lon, glob_abs),
        "photo_equator_signed_br": interp_to_360(lon, eq_signed),
        "photo_midlat_signed_br": interp_to_360(lon, mid_signed),
        "photo_global_signed_br": interp_to_360(lon, glob_signed),
    })

def candidate_manifest_paths(cr):
    mf = b / "metadata" / "rotation_manifest.csv"
    if not mf.exists():
        return []

    m = pd.read_csv(mf)
    if "cr" not in m.columns:
        return []

    m["cr_num"] = pd.to_numeric(m["cr"], errors="coerce")
    q = m[m["cr_num"] == int(cr)]

    paths = []
    for col in q.columns:
        cl = col.lower()
        if not any(k in cl for k in ["path", "file", "fits", "magnetogram", "hmi"]):
            continue

        for val in q[col].dropna().astype(str):
            if len(val.strip()) == 0:
                continue
            p = Path(val)
            if not p.is_absolute():
                p = b / p
            paths.append(p)

    return paths

def find_hmi_fits(cr):
    candidates = []

    for p in candidate_manifest_paths(cr):
        candidates.append(p)

    roots = [
        b / "data" / "raw" / "magnetograms" / "hmi",
        b / "data" / "processed" / "magnetograms" / "hmi",
        b / "data" / "raw",
        b / "data" / "processed",
        b / "data",
    ]

    for root in roots:
        if not root.exists():
            continue

        for p in root.rglob("*"):
            if not p.is_file():
                continue

            name = p.name.lower()
            full = str(p).lower()

            if str(cr) not in name:
                continue
            if not (name.endswith(".fits") or name.endswith(".fit") or name.endswith(".fits.gz")):
                continue
            if "source_surface" in name:
                continue
            if "metrics" in name:
                continue
            if "longitude" in name:
                continue
            if "/pfss/" in full:
                continue

            candidates.append(p)

    checked = set()
    for p in candidates:
        p = Path(p)
        if p in checked:
            continue
        checked.add(p)

        if not p.exists():
            continue

        try:
            arr = read_2d_fits(p)
            if arr.shape[0] >= 60 and arr.shape[1] >= 120:
                return p
        except Exception:
            pass

    raise RuntimeError(f"Could not find usable HMI magnetogram FITS for CR {cr}")

def corr(x, y, method):
    q = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(q) < 3:
        return np.nan, len(q)
    if q["x"].nunique() < 2 or q["y"].nunique() < 2:
        return np.nan, len(q)
    return float(q["x"].corr(q["y"], method=method)), len(q)

photo_cache = {}
expansion_profile_rows = []
matched_rows = []
binned_rows = []
corr_rows = []

for rss in heights:
    h = f"rss{rss:.1f}"

    pf = b / "data" / "processed" / "pfss" / "hmi_longitude_proxies" / f"pfss_hmi_selected_crs_rss{rss:.1f}_longitude_profiles_all.csv"
    mf = b / "data" / "processed" / "comparison" / h / "pfss_omni_ballistic_matched_rows.csv"

    if not pf.exists():
        raise RuntimeError(f"Missing PFSS longitude profile table: {pf}")
    if not mf.exists():
        raise RuntimeError(f"Missing ballistic matched rows table: {mf}")

    p = pd.read_csv(pf)
    j = pd.read_csv(mf)

    needed_pfss = [
        "cr",
        "lon_deg",
        "equator_abs_br",
        "midlat_abs_br",
        "global_abs_br",
        "equator_signed_br",
        "global_signed_br",
    ]

    needed_match = [
        "cr",
        "ballistic_lon_deg",
        "ballistic_phase10_deg",
        "speed_km_s",
        "bmag_nt",
        "density_cm3",
        "pdyn_npa",
        "ma",
    ]

    miss_pfss = [c for c in needed_pfss if c not in p.columns]
    miss_match = [c for c in needed_match if c not in j.columns]

    if miss_pfss:
        raise RuntimeError(f"{pf} missing columns: {miss_pfss}")
    if miss_match:
        raise RuntimeError(f"{mf} missing columns: {miss_match}")

    for c in needed_pfss:
        p[c] = pd.to_numeric(p[c], errors="coerce")

    for c in needed_match:
        j[c] = pd.to_numeric(j[c], errors="coerce")

    p["cr"] = p["cr"].astype(int)
    j["cr"] = j["cr"].astype(int)

    p["rss"] = rss
    p["lon_bin"] = np.floor(np.mod(p["lon_deg"], 360.0)).astype(int)

    per_height_profiles = []

    for cr in sorted(p["cr"].dropna().astype(int).unique()):
        if cr not in photo_cache:
            hmi_path = find_hmi_fits(cr)
            arr = read_2d_fits(hmi_path)
            prof = magnetogram_profile(arr)
            prof["cr"] = cr
            prof["hmi_file"] = str(hmi_path)
            photo_cache[cr] = prof

        phot = photo_cache[cr]
        q = p[p["cr"] == cr].copy()

        e = q.merge(phot, on=["cr", "lon_bin"], how="left")

        for band in ["equator", "midlat", "global"]:
            ss = pd.to_numeric(e[f"{band}_abs_br"], errors="coerce")
            ph = pd.to_numeric(e[f"photo_{band}_abs_br"], errors="coerce")

            fproxy = ph / ((rss ** 2.0) * ss.replace(0.0, np.nan))
            e[f"radial_expansion_proxy_{band}"] = fproxy
            e[f"log_radial_expansion_proxy_{band}"] = np.log10(fproxy.where(fproxy > 0.0))
            e[f"expansion_speed_proxy_{band}"] = -e[f"log_radial_expansion_proxy_{band}"]

        per_height_profiles.append(e)

    ep = pd.concat(per_height_profiles, ignore_index=True)
    expansion_profile_rows.append(ep)

    keep_cols = [
        "rss",
        "cr",
        "lon_bin",
        "photo_equator_abs_br",
        "photo_midlat_abs_br",
        "photo_global_abs_br",
        "radial_expansion_proxy_equator",
        "radial_expansion_proxy_midlat",
        "radial_expansion_proxy_global",
        "log_radial_expansion_proxy_equator",
        "log_radial_expansion_proxy_midlat",
        "log_radial_expansion_proxy_global",
        "expansion_speed_proxy_equator",
        "expansion_speed_proxy_midlat",
        "expansion_speed_proxy_global",
    ]

    jp = ep[keep_cols].copy()

    j["rss"] = rss
    j["lon_bin"] = np.floor(np.mod(j["ballistic_lon_deg"], 360.0)).astype(int)

    jm = j.merge(jp, on=["rss", "cr", "lon_bin"], how="left")
    matched_rows.append(jm)

    g = jm.groupby(["rss", "cr", "ballistic_phase10_deg"], as_index=False).agg(
        n=("cr", "size"),
        speed_mean=("speed_km_s", "mean"),
        bmag_mean=("bmag_nt", "mean"),
        density_mean=("density_cm3", "mean"),
        pdyn_mean=("pdyn_npa", "mean"),
        ma_mean=("ma", "mean"),
        expansion_speed_proxy_equator_mean=("expansion_speed_proxy_equator", "mean"),
        expansion_speed_proxy_midlat_mean=("expansion_speed_proxy_midlat", "mean"),
        expansion_speed_proxy_global_mean=("expansion_speed_proxy_global", "mean"),
        radial_expansion_proxy_equator_mean=("radial_expansion_proxy_equator", "mean"),
        radial_expansion_proxy_midlat_mean=("radial_expansion_proxy_midlat", "mean"),
        radial_expansion_proxy_global_mean=("radial_expansion_proxy_global", "mean"),
    )

    binned_rows.append(g)

all_profiles = pd.concat(expansion_profile_rows, ignore_index=True)
all_matched = pd.concat(matched_rows, ignore_index=True)
all_binned = pd.concat(binned_rows, ignore_index=True)

xcols = [
    "expansion_speed_proxy_equator_mean",
    "expansion_speed_proxy_midlat_mean",
    "expansion_speed_proxy_global_mean",
    "radial_expansion_proxy_equator_mean",
    "radial_expansion_proxy_midlat_mean",
    "radial_expansion_proxy_global_mean",
]

ycols = [
    "speed_mean",
    "bmag_mean",
    "density_mean",
    "pdyn_mean",
    "ma_mean",
]

for rss in heights:
    q = all_binned[all_binned["rss"] == rss].copy()

    for x in xcols:
        for y in ycols:
            for method in ["spearman", "pearson"]:
                r, n = corr(q[x], q[y], method)
                corr_rows.append({
                    "rss": rss,
                    "cr": "all",
                    "x": x,
                    "y": y,
                    "method": method,
                    "r": r,
                    "n": n,
                })

    for cr in sorted(q["cr"].dropna().astype(int).unique()):
        qc = q[q["cr"].astype(int) == cr].copy()

        for x in xcols:
            for y in ["speed_mean", "bmag_mean"]:
                for method in ["spearman"]:
                    r, n = corr(qc[x], qc[y], method)
                    corr_rows.append({
                        "rss": rss,
                        "cr": cr,
                        "x": x,
                        "y": y,
                        "method": method,
                        "r": r,
                        "n": n,
                    })

corr = pd.DataFrame(corr_rows)

profiles_path = outtab / "pfss_omni_radial_expansion_proxy_profiles.csv"
matched_path = outtab / "pfss_omni_radial_expansion_proxy_matched_rows.csv"
binned_path = outtab / "pfss_omni_radial_expansion_proxy_binned.csv"
corr_path = outtab / "pfss_omni_radial_expansion_proxy_correlations.csv"

all_profiles.to_csv(profiles_path, index=False)
all_matched.to_csv(matched_path, index=False)
all_binned.to_csv(binned_path, index=False)
corr.to_csv(corr_path, index=False)

target = corr[
    (corr["cr"].astype(str) == "all")
    & (corr["x"] == "expansion_speed_proxy_equator_mean")
    & (corr["y"] == "speed_mean")
    & (corr["method"] == "spearman")
].copy()

plt.figure(figsize=(7, 5))
plt.plot(target["rss"], target["r"], marker="o", linewidth=2)
plt.axhline(0.0, linewidth=1)
plt.xlabel("PFSS source-surface height, Rs")
plt.ylabel("Spearman r")
plt.title("Expansion speed proxy vs OMNI speed")
plt.xticks(target["rss"])
plt.tight_layout()
fig_path = outfig / "fig_expansion_speed_proxy_vs_omni_speed_by_height.png"
plt.savefig(fig_path, dpi=300)
plt.close()

best = corr[(corr["cr"].astype(str) == "all") & corr["r"].notna()].copy()
best["abs_r"] = best["r"].abs()
best = best.sort_values("abs_r", ascending=False).head(20)

txt = []
txt.append("Radial expansion proxy summary")
txt.append("")
txt.append("Definition:")
txt.append("f_proxy = |B_photosphere| / (rss^2 |B_source_surface|), computed in longitude bins.")
txt.append("expansion_speed_proxy = -log10(f_proxy), so larger values correspond to smaller expansion proxy.")
txt.append("This is not a true traced WSA flux-tube expansion factor.")
txt.append("")
txt.append("Target relation:")
txt.append(target.to_string(index=False))
txt.append("")
txt.append("Best all-height, all-CR correlations:")
txt.append(best[["rss", "x", "y", "method", "r", "n"]].to_string(index=False))

txt_path = outres / "pfss_omni_radial_expansion_proxy_summary.txt"
txt_path.write_text("\n".join(txt))

print("Radial expansion proxy target relation:")
print(target.to_string(index=False))
print()
print("Best all-height, all-CR correlations:")
print(best[["rss", "x", "y", "method", "r", "n"]].to_string(index=False))
print()
print(f"Saved profiles: {profiles_path}")
print(f"Saved matched rows: {matched_path}")
print(f"Saved binned table: {binned_path}")
print(f"Saved correlations: {corr_path}")
print(f"Saved figure: {fig_path}")
print(f"Saved summary: {txt_path}")
print()
print("Status: pass")
