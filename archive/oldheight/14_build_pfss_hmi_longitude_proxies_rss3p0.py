from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sunpy.map

b = Path(__file__).resolve().parents[1]
r = 3.0

man = b / "metadata" / "rotation_manifest.csv"
pfd = b / "data" / "processed" / "pfss" / "hmi"
od = b / "data" / "processed" / "pfss" / "hmi_longitude_proxies"
fd = b / "outputs" / "figures" / "pfss" / "hmi"
td = b / "outputs" / "tables"

od.mkdir(parents=True, exist_ok=True)
fd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

d = pd.read_csv(man)

if "proxy_status" not in d.columns:
    d["proxy_status"] = ""

for c in ["magnetogram_status", "omni_status", "pfss_status", "mapping_status", "notes", "magnetogram_source", "proxy_status"]:
    if c not in d.columns:
        d[c] = ""
    d[c] = d[c].fillna("").astype(str)

crs = sorted(d["cr"].dropna().astype(int).unique())

print(f"Project base: {b}")
print(f"Building PFSS longitude proxies for rss={r}")
print(f"Carrington rotations: {crs}")
print()

rows = []
plist = []

for cr in crs:
    print("=" * 70)
    print(f"CR {cr}")

    fp = pfd / f"pfss_hmi_cr{cr}_rss{r:.1f}_source_surface_br.fits"

    if not fp.exists():
        print(f"Missing PFSS source-surface file: {fp}")
        d.loc[d["cr"] == cr, "proxy_status"] = "fail_missing_pfss"
        continue

    sm = sunpy.map.Map(fp)
    br = np.asarray(sm.data, dtype=float)

    ff = np.isfinite(br).sum() / br.size

    if ff < 0.99:
        print(f"Bad finite fraction: {ff:.3f}")
        d.loc[d["cr"] == cr, "proxy_status"] = "fail_bad_finite_fraction"
        continue

    ny, nx = br.shape

    lat = np.linspace(-90, 90, ny)
    lon = np.linspace(0, 360, nx, endpoint=False)

    em = np.abs(lat) <= 20
    mm = np.abs(lat) <= 45

    es = np.nanmean(br[em, :], axis=0)
    ea = np.nanmean(np.abs(br[em, :]), axis=0)

    ms = np.nanmean(br[mm, :], axis=0)
    ma = np.nanmean(np.abs(br[mm, :]), axis=0)

    gs = np.nanmean(br, axis=0)
    ga = np.nanmean(np.abs(br), axis=0)

    p = pd.DataFrame({
        "cr": cr,
        "rss": r,
        "lon_deg": lon,
        "equator_signed_br": es,
        "equator_abs_br": ea,
        "midlat_signed_br": ms,
        "midlat_abs_br": ma,
        "global_signed_br": gs,
        "global_abs_br": ga,
        "equator_polarity": np.sign(es),
        "global_polarity": np.sign(gs),
    })

    pp = od / f"pfss_hmi_cr{cr}_rss{r:.1f}_longitude_proxy.csv"
    p.to_csv(pp, index=False)

    print(f"Saved longitude proxy: {pp}")
    print(f"Rows: {len(p)}")
    print(f"Mean equator abs Br: {p['equator_abs_br'].mean():.6f}")
    print(f"Mean global abs Br: {p['global_abs_br'].mean():.6f}")

    rows.append({
        "cr": cr,
        "rss": r,
        "nlon": nx,
        "nlat": ny,
        "finite_fraction": ff,
        "equator_abs_mean": float(p["equator_abs_br"].mean()),
        "equator_abs_std_by_lon": float(p["equator_abs_br"].std()),
        "equator_abs_min": float(p["equator_abs_br"].min()),
        "equator_abs_max": float(p["equator_abs_br"].max()),
        "equator_abs_range": float(p["equator_abs_br"].max() - p["equator_abs_br"].min()),
        "global_abs_mean": float(p["global_abs_br"].mean()),
        "global_abs_std_by_lon": float(p["global_abs_br"].std()),
        "global_abs_min": float(p["global_abs_br"].min()),
        "global_abs_max": float(p["global_abs_br"].max()),
        "global_abs_range": float(p["global_abs_br"].max() - p["global_abs_br"].min()),
    })

    plist.append(p)

    m = d["cr"] == cr
    d.loc[m, "proxy_status"] = "pass_pfss_longitude_proxy"

    old = d.loc[m, "notes"].iloc[0]
    note = "PFSS source-surface longitude proxy created."

    if note not in old:
        d.loc[m, "notes"] = old + " | " + note

if len(rows) == 0:
    raise RuntimeError("No PFSS longitude proxies were created.")

s = pd.DataFrame(rows).sort_values("cr")

sp = od / f"pfss_hmi_selected_crs_rss{r:.1f}_longitude_proxy_summary.csv"
s.to_csv(sp, index=False)

cp = td / "pfss_hmi_longitude_proxy_summary_compact_rss3.0.csv"
s.to_csv(cp, index=False)

d.to_csv(man, index=False)

ap = pd.concat(plist, ignore_index=True)

app = od / f"pfss_hmi_selected_crs_rss{r:.1f}_longitude_profiles_all.csv"
ap.to_csv(app, index=False)

fig = fd / "pfss_hmi_equator_abs_br_by_longitude.png"

plt.figure(figsize=(10, 5))

for cr in crs:
    q = ap[ap["cr"] == cr]
    if len(q) == 0:
        continue
    plt.plot(q["lon_deg"], q["equator_abs_br"], label=f"CR {cr}")

plt.xlabel("Carrington longitude index proxy, degrees")
plt.ylabel("Mean unsigned source-surface Br near equator")
plt.title("PFSS HMI source-surface longitude proxy")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(fig, dpi=200)
plt.close()

print()
print("Longitude proxy summary:")
print(s.to_string(index=False))

print()
print(f"Saved summary: {sp}")
print(f"Saved compact table: {cp}")
print(f"Saved all profiles: {app}")
print(f"Saved figure: {fig}")

print()
print("Manifest proxy status:")
print(d[["cr", "pfss_status", "proxy_status"]].to_string(index=False))

print()
print("Status: pass")
