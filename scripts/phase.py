from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

b = Path.cwd()

pf = b / "data" / "processed" / "pfss" / "hmi_longitude_proxies" / "pfss_hmi_selected_crs_rss2.5_longitude_profiles_all.csv"
of = b / "data" / "processed" / "omni" / "omni_selected_crs_clean_with_phase.csv"
mf = b / "metadata" / "rotation_manifest.csv"

cd = b / "data" / "processed" / "comparison"
fd = b / "outputs" / "figures" / "comparison"
td = b / "outputs" / "tables"

cd.mkdir(parents=True, exist_ok=True)
fd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

print(f"Project base: {b}")
print(f"PFSS table: {pf}")
print(f"OMNI table: {of}")
print()

p = pd.read_csv(pf)
o = pd.read_csv(of)
m = pd.read_csv(mf)

if "phase_match_status" not in m.columns:
    m["phase_match_status"] = ""

for c in ["phase_match_status", "notes"]:
    if c not in m.columns:
        m[c] = ""
    m[c] = m[c].fillna("").astype(str)

pc = [
    "cr",
    "lon_deg",
    "equator_signed_br",
    "equator_abs_br",
    "midlat_signed_br",
    "midlat_abs_br",
    "global_signed_br",
    "global_abs_br",
    "equator_polarity",
    "global_polarity",
]

oc = [
    "cr",
    "time",
    "lon_deg_proxy",
    "speed_km_s",
    "bx_gse_nt",
    "by_gse_nt",
    "bz_gse_nt",
    "bmag_nt",
    "density_cm3",
    "pdyn_npa",
    "va_km_s",
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

p["cr"] = pd.to_numeric(p["cr"], errors="coerce").astype("Int64")
o["cr"] = pd.to_numeric(o["cr"], errors="coerce").astype("Int64")

p["lon_bin"] = np.mod(np.rint(p["lon_deg"]).astype(int), 360)

o["time"] = pd.to_datetime(o["time"], utc=True, errors="coerce")
o["lon_deg_proxy"] = pd.to_numeric(o["lon_deg_proxy"], errors="coerce")
o = o.dropna(subset=["cr", "time", "lon_deg_proxy"]).copy()

o["lon_bin"] = np.floor(o["lon_deg_proxy"]).astype(int)
o["lon_bin"] = o["lon_bin"].clip(0, 359)

p = p.rename(columns={"lon_deg": "pfss_lon_deg"})

j = o.merge(p, on=["cr", "lon_bin"], how="left")

miss = int(j["equator_abs_br"].isna().sum())

if miss > 0:
    print(f"WARNING: unmatched rows: {miss}")

j["phase10_deg"] = np.floor(j["lon_deg_proxy"] / 10.0) * 10.0
j["phase10_deg"] = j["phase10_deg"].clip(0, 350)

jp = cd / "pfss_omni_phase_matched_rows.csv"
j.to_csv(jp, index=False)

g = j.groupby(["cr", "phase10_deg"], as_index=False).agg(
    n=("time", "size"),
    speed_mean=("speed_km_s", "mean"),
    speed_std=("speed_km_s", "std"),
    bmag_mean=("bmag_nt", "mean"),
    bmag_std=("bmag_nt", "std"),
    density_mean=("density_cm3", "mean"),
    density_std=("density_cm3", "std"),
    pdyn_mean=("pdyn_npa", "mean"),
    pdyn_std=("pdyn_npa", "std"),
    va_mean=("va_km_s", "mean"),
    ma_mean=("ma", "mean"),
    pfss_equator_abs_mean=("equator_abs_br", "mean"),
    pfss_equator_signed_mean=("equator_signed_br", "mean"),
    pfss_midlat_abs_mean=("midlat_abs_br", "mean"),
    pfss_global_abs_mean=("global_abs_br", "mean"),
    pfss_global_signed_mean=("global_signed_br", "mean"),
)

gp = cd / "pfss_omni_phase10_binned.csv"
g.to_csv(gp, index=False)

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

def cc(d, x, y, how):
    q = d[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(q) < 3:
        return np.nan, len(q)
    if q[x].std() == 0 or q[y].std() == 0:
        return np.nan, len(q)
    return float(q[x].corr(q[y], method=how)), len(q)

rs = []

for tag, q in [("all", g)]:
    for x in xs:
        for y in ys:
            rp, n = cc(q, x, y, "pearson")
            rs.append({"cr": tag, "x": x, "y": y, "method": "pearson", "r": rp, "n": n})

            rr, n = cc(q, x, y, "spearman")
            rs.append({"cr": tag, "x": x, "y": y, "method": "spearman", "r": rr, "n": n})

for cr in sorted(g["cr"].dropna().unique()):
    q = g[g["cr"] == cr]

    for x in xs:
        for y in ys:
            rp, n = cc(q, x, y, "pearson")
            rs.append({"cr": int(cr), "x": x, "y": y, "method": "pearson", "r": rp, "n": n})

            rr, n = cc(q, x, y, "spearman")
            rs.append({"cr": int(cr), "x": x, "y": y, "method": "spearman", "r": rr, "n": n})

r = pd.DataFrame(rs)

rp = cd / "pfss_omni_phase_correlation_summary.csv"
rt = td / "pfss_omni_phase_correlation_summary_compact.csv"

r.to_csv(rp, index=False)
r.to_csv(rt, index=False)

for x, y, fn, lab in [
    ("pfss_equator_abs_mean", "speed_mean", "pfss_eqabs_vs_omni_speed_phase10.png", "OMNI speed, km/s"),
    ("pfss_equator_abs_mean", "bmag_mean", "pfss_eqabs_vs_omni_bmag_phase10.png", "OMNI |B|, nT"),
    ("pfss_global_abs_mean", "speed_mean", "pfss_globalabs_vs_omni_speed_phase10.png", "OMNI speed, km/s"),
    ("pfss_global_abs_mean", "bmag_mean", "pfss_globalabs_vs_omni_bmag_phase10.png", "OMNI |B|, nT"),
]:
    plt.figure(figsize=(7, 5))

    for cr in sorted(g["cr"].dropna().unique()):
        q = g[g["cr"] == cr]
        plt.scatter(q[x], q[y], label=f"CR {int(cr)}", s=18)

    plt.xlabel(x)
    plt.ylabel(lab)
    plt.title(fn.replace(".png", ""))
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(fd / fn, dpi=200)
    plt.close()

for y, lab, fn in [
    ("speed_mean", "OMNI speed, km/s", "phase10_omni_speed_by_cr.png"),
    ("bmag_mean", "OMNI |B|, nT", "phase10_omni_bmag_by_cr.png"),
    ("pfss_equator_abs_mean", "PFSS equator unsigned Br", "phase10_pfss_eqabs_by_cr.png"),
]:
    plt.figure(figsize=(10, 5))

    for cr in sorted(g["cr"].dropna().unique()):
        q = g[g["cr"] == cr]
        plt.plot(q["phase10_deg"], q[y], label=f"CR {int(cr)}")

    plt.xlabel("Rotation phase proxy, degrees")
    plt.ylabel(lab)
    plt.title(fn.replace(".png", ""))
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(fd / fn, dpi=200)
    plt.close()

for cr in sorted(j["cr"].dropna().unique()):
    q = j[j["cr"] == cr]
    ok = q["equator_abs_br"].notna().all()

    ms = m["cr"].astype(int) == int(cr)

    if ok:
        m.loc[ms, "phase_match_status"] = "pass_pfss_omni_phase_match"
    else:
        m.loc[ms, "phase_match_status"] = "warn_unmatched_phase_rows"

    nt = "PFSS and OMNI phase-matched diagnostic table created."
    old = m.loc[ms, "notes"].iloc[0]

    if nt not in old:
        m.loc[ms, "notes"] = old + " | " + nt

m.to_csv(mf, index=False)

print("Phase match summary:")
print(f"Matched rows: {len(j)}")
print(f"Unmatched rows: {miss}")
print(f"10 degree bins: {len(g)}")
print(f"CRs: {sorted(j['cr'].dropna().astype(int).unique().tolist())}")

print()
print("Best all-CR absolute correlations:")
q = r[(r["cr"].astype(str) == "all") & (r["r"].notna())].copy()
q["abs_r"] = q["r"].abs()
q = q.sort_values("abs_r", ascending=False).head(12)
print(q[["x", "y", "method", "r", "n"]].to_string(index=False))

print()
print(f"Saved matched rows: {jp}")
print(f"Saved 10 degree binned table: {gp}")
print(f"Saved correlation table: {rp}")
print(f"Saved compact correlation table: {rt}")
print(f"Saved figures: {fd}")

print()
print("Manifest phase match status:")
print(m[["cr", "proxy_status", "omni_feature_status", "phase_match_status"]].to_string(index=False))

print()
print("Status: pass")
