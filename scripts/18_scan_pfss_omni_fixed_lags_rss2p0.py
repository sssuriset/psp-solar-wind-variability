from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

b = Path.cwd()

pf = b / "data" / "processed" / "pfss" / "hmi_longitude_proxies" / "pfss_hmi_selected_crs_rss2.0_longitude_profiles_all.csv"
of = b / "data" / "processed" / "omni" / "omni_selected_crs_clean_with_phase.csv"
bf = b / "data" / "processed" / "comparison" / "rss2.0" / "pfss_omni_ballistic_correlation_summary.csv"
mf = b / "metadata" / "rotation_manifest.csv"

cd = b / "data" / "processed" / "comparison" / "rss2.0"
fd = b / "outputs" / "figures" / "comparison" / "rss2.0"
td = b / "outputs" / "tables" / "rss2.0"

cd.mkdir(parents=True, exist_ok=True)
fd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

p = pd.read_csv(pf)
o = pd.read_csv(of)
br = pd.read_csv(bf)
m = pd.read_csv(mf)

if "lag_scan_status" not in m.columns:
    m["lag_scan_status"] = ""

for c in ["lag_scan_status", "notes"]:
    if c not in m.columns:
        m[c] = ""
    m[c] = m[c].fillna("").astype(str)

pc = [
    "cr",
    "lon_deg",
    "equator_abs_br",
    "midlat_abs_br",
    "global_abs_br",
]

oc = [
    "cr",
    "time",
    "lon_deg_proxy",
    "speed_km_s",
    "bmag_nt",
    "density_cm3",
    "pdyn_npa",
    "ma",
]

pm = [c for c in pc if c not in p.columns]
om = [c for c in oc if c not in o.columns]

if len(pm) > 0:
    raise RuntimeError(f"Missing PFSS columns: {pm}")

if len(om) > 0:
    raise RuntimeError(f"Missing OMNI columns: {om}")

p = p[pc].copy()
o = o[oc].copy()

p["cr"] = pd.to_numeric(p["cr"], errors="coerce").astype(int)
o["cr"] = pd.to_numeric(o["cr"], errors="coerce").astype(int)

p["lon_bin"] = np.mod(np.rint(p["lon_deg"]).astype(int), 360)
p = p.rename(columns={"lon_deg": "pfss_lon_deg"})

for c in oc:
    if c != "time":
        o[c] = pd.to_numeric(o[c], errors="coerce")

o["time"] = pd.to_datetime(o["time"], utc=True, errors="coerce")
o = o.dropna(subset=["cr", "time", "lon_deg_proxy"]).copy()

xs = [
    "pfss_equator_abs_mean",
    "pfss_midlat_abs_mean",
    "pfss_global_abs_mean",
]

ys = [
    "speed_mean",
    "bmag_mean",
    "density_mean",
    "pdyn_mean",
    "ma_mean",
]

def cor(d, x, y, md):
    q = d[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()

    if len(q) < 3:
        return np.nan, len(q)

    if q[x].std() == 0 or q[y].std() == 0:
        return np.nan, len(q)

    return float(q[x].corr(q[y], method=md)), len(q)

gr = []
rr = []

print(f"Project base: {b}")
print("Scanning fixed longitude shifts from 0 to 120 deg")
print()

for sh in range(0, 121, 5):
    q = o.copy()

    q["scan_shift_deg"] = float(sh)
    q["scan_lon_deg"] = np.mod(q["lon_deg_proxy"] - sh, 360.0)
    q["lon_bin"] = np.floor(q["scan_lon_deg"]).astype(int).clip(0, 359)

    j = q.merge(p, on=["cr", "lon_bin"], how="left")
    j["scan_phase10_deg"] = np.floor(j["scan_lon_deg"] / 10.0) * 10.0
    j["scan_phase10_deg"] = j["scan_phase10_deg"].clip(0, 350)

    g = j.groupby(["cr", "scan_phase10_deg"], as_index=False).agg(
        n=("time", "size"),
        speed_mean=("speed_km_s", "mean"),
        bmag_mean=("bmag_nt", "mean"),
        density_mean=("density_cm3", "mean"),
        pdyn_mean=("pdyn_npa", "mean"),
        ma_mean=("ma", "mean"),
        pfss_equator_abs_mean=("equator_abs_br", "mean"),
        pfss_midlat_abs_mean=("midlat_abs_br", "mean"),
        pfss_global_abs_mean=("global_abs_br", "mean"),
    )

    g["scan_shift_deg"] = float(sh)
    gr.append(g)

    for x in xs:
        for y in ys:
            for md in ["pearson", "spearman"]:
                rv, n = cor(g, x, y, md)
                rr.append({
                    "scan_shift_deg": float(sh),
                    "cr": "all",
                    "x": x,
                    "y": y,
                    "method": md,
                    "r": rv,
                    "n": n,
                })

    for cr in sorted(g["cr"].dropna().astype(int).unique()):
        z = g[g["cr"] == cr]

        for x in xs:
            for y in ys:
                for md in ["pearson", "spearman"]:
                    rv, n = cor(z, x, y, md)
                    rr.append({
                        "scan_shift_deg": float(sh),
                        "cr": cr,
                        "x": x,
                        "y": y,
                        "method": md,
                        "r": rv,
                        "n": n,
                    })

ga = pd.concat(gr, ignore_index=True)
ra = pd.DataFrame(rr)

gp = cd / "pfss_omni_fixed_lag_scan_binned.csv"
rp = cd / "pfss_omni_fixed_lag_scan_correlations.csv"
rt = td / "pfss_omni_fixed_lag_scan_correlations_compact.csv"

ga.to_csv(gp, index=False)
ra.to_csv(rp, index=False)
ra.to_csv(rt, index=False)

qa = ra[(ra["cr"].astype(str) == "all") & (ra["r"].notna())].copy()
qa["abs_r"] = qa["r"].abs()

best = qa.sort_values("abs_r", ascending=False).head(15)

bp = cd / "pfss_omni_fixed_lag_scan_best_all.csv"
bt = td / "pfss_omni_fixed_lag_scan_best_all_compact.csv"

best.to_csv(bp, index=False)
best.to_csv(bt, index=False)

bal = br[(br["cr"].astype(str) == "all") & (br["r"].notna())].copy()
bal["abs_r"] = bal["r"].abs()

bc = bal.sort_values("abs_r", ascending=False).head(15)

for x, y, md, fn in [
    ("pfss_equator_abs_mean", "bmag_mean", "spearman", "scan_eqabs_bmag_spearman_by_shift.png"),
    ("pfss_equator_abs_mean", "bmag_mean", "pearson", "scan_eqabs_bmag_pearson_by_shift.png"),
    ("pfss_equator_abs_mean", "pdyn_mean", "spearman", "scan_eqabs_pdyn_spearman_by_shift.png"),
    ("pfss_equator_abs_mean", "density_mean", "spearman", "scan_eqabs_density_spearman_by_shift.png"),
    ("pfss_midlat_abs_mean", "bmag_mean", "spearman", "scan_midabs_bmag_spearman_by_shift.png"),
]:
    q = qa[(qa["x"] == x) & (qa["y"] == y) & (qa["method"] == md)].copy()

    plt.figure(figsize=(8, 5))
    plt.plot(q["scan_shift_deg"], q["r"], marker="o")
    plt.axhline(0, linewidth=1)
    plt.xlabel("Fixed longitude shift, degrees")
    plt.ylabel(f"{md} r")
    plt.title(f"{x} vs {y}")
    plt.tight_layout()
    plt.savefig(fd / fn, dpi=200)
    plt.close()

for cr in sorted(o["cr"].dropna().astype(int).unique()):
    ms = m["cr"].astype(int) == cr
    m.loc[ms, "lag_scan_status"] = "pass_fixed_lag_scan"

    nt = "Fixed longitude lag scan completed."
    old = m.loc[ms, "notes"].iloc[0]

    if nt not in old:
        m.loc[ms, "notes"] = old + " | " + nt

m.to_csv(mf, index=False)

print("Fixed lag scan summary:")
print(f"Shift count: {len(sorted(ra['scan_shift_deg'].unique()))}")
print(f"Binned rows: {len(ga)}")
print(f"Correlation rows: {len(ra)}")

print()
print("Best fixed-lag all-CR correlations:")
print(best[["scan_shift_deg", "x", "y", "method", "r", "n"]].to_string(index=False))

print()
print("Best ballistic all-CR correlations:")
print(bc[["x", "y", "method", "r", "n"]].to_string(index=False))

print()
print(f"Saved binned scan table: {gp}")
print(f"Saved scan correlation table: {rp}")
print(f"Saved compact scan table: {rt}")
print(f"Saved best scan table: {bp}")
print(f"Saved figures: {fd}")

print()
print("Manifest lag scan status:")
print(m[["cr", "ballistic_match_status", "lag_scan_status"]].to_string(index=False))

print()
print("Status: pass")
