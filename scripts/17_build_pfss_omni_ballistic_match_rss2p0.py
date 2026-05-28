from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

b = Path.cwd()

pf = b / "data" / "processed" / "pfss" / "hmi_longitude_proxies" / "pfss_hmi_selected_crs_rss2.0_longitude_profiles_all.csv"
of = b / "data" / "processed" / "omni" / "omni_selected_crs_clean_with_phase.csv"
mf = b / "metadata" / "rotation_manifest.csv"

cd = b / "data" / "processed" / "comparison" / "rss2.0"
fd = b / "outputs" / "figures" / "comparison" / "rss2.0"
td = b / "outputs" / "tables" / "rss2.0"

cd.mkdir(parents=True, exist_ok=True)
fd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

au = 149597870.7
rs = 695700.0
rss = 2.0
dsec = 86400.0

p = pd.read_csv(pf)
o = pd.read_csv(of)
m = pd.read_csv(mf)

if "ballistic_match_status" not in m.columns:
    m["ballistic_match_status"] = ""

for c in ["ballistic_match_status", "notes"]:
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
    "tday",
    "rot_frac",
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

p["cr"] = pd.to_numeric(p["cr"], errors="coerce").astype(int)
o["cr"] = pd.to_numeric(o["cr"], errors="coerce").astype(int)

p["lon_bin"] = np.mod(np.rint(p["lon_deg"]).astype(int), 360)
p = p.rename(columns={"lon_deg": "pfss_lon_deg"})

o["time"] = pd.to_datetime(o["time"], utc=True, errors="coerce")

for c in oc:
    if c not in ["time"]:
        o[c] = pd.to_numeric(o[c], errors="coerce")

o = o.dropna(subset=["cr", "time", "lon_deg_proxy", "speed_km_s"]).copy()

js = []
ls = []

print(f"Project base: {b}")
print(f"PFSS table: {pf}")
print(f"OMNI table: {of}")
print()

for cr in sorted(o["cr"].dropna().astype(int).unique()):
    print("=" * 80)
    print(f"CR {cr}")

    q = o[o["cr"] == cr].copy()
    pp = p[p["cr"] == cr].copy()

    if len(pp) == 0:
        print("No PFSS rows for this CR.")
        m.loc[m["cr"].astype(int) == cr, "ballistic_match_status"] = "fail_no_pfss_rows"
        continue

    x = q[(q["lon_deg_proxy"] > 1) & (q["tday"] > 0)].copy()
    pr = 27.2753

    if len(x) > 0:
        vv = 360.0 * x["tday"] / x["lon_deg_proxy"]
        vv = vv.replace([np.inf, -np.inf], np.nan).dropna()
        vv = vv[(vv > 20) & (vv < 40)]

        if len(vv) > 0:
            pr = float(vv.median())

    v = q["speed_km_s"].copy()
    v[v <= 0] = np.nan

    tau = ((au - rss * rs) / v) / dsec
    sd = 360.0 * tau / pr
    blon = np.mod(q["lon_deg_proxy"] - sd, 360.0)

    q["ballistic_lag_days"] = tau
    q["ballistic_shift_deg"] = sd
    q["ballistic_lon_deg"] = blon

    q = q.dropna(subset=["ballistic_lon_deg"]).copy()
    q["lon_bin"] = np.floor(q["ballistic_lon_deg"]).astype(int).clip(0, 359)

    j = q.merge(pp, on=["cr", "lon_bin"], how="left")

    miss = int(j["equator_abs_br"].isna().sum())

    j["ballistic_phase10_deg"] = np.floor(j["ballistic_lon_deg"] / 10.0) * 10.0
    j["ballistic_phase10_deg"] = j["ballistic_phase10_deg"].clip(0, 350)

    js.append(j)

    ls.append({
        "cr": cr,
        "rows": int(len(j)),
        "period_days": pr,
        "lag_mean_days": float(j["ballistic_lag_days"].mean()),
        "lag_min_days": float(j["ballistic_lag_days"].min()),
        "lag_max_days": float(j["ballistic_lag_days"].max()),
        "shift_mean_deg": float(j["ballistic_shift_deg"].mean()),
        "shift_min_deg": float(j["ballistic_shift_deg"].min()),
        "shift_max_deg": float(j["ballistic_shift_deg"].max()),
        "unmatched_rows": miss,
    })

    ms = m["cr"].astype(int) == cr

    if miss == 0:
        m.loc[ms, "ballistic_match_status"] = "pass_pfss_omni_ballistic_match"
    else:
        m.loc[ms, "ballistic_match_status"] = "warn_unmatched_ballistic_rows"

    nt = "PFSS and OMNI ballistic-lag diagnostic table created."
    old = m.loc[ms, "notes"].iloc[0]

    if nt not in old:
        m.loc[ms, "notes"] = old + " | " + nt

    print(f"Rows matched: {len(j)}")
    print(f"Unmatched rows: {miss}")
    print(f"Period used: {pr:.4f} days")
    print(f"Mean lag: {j['ballistic_lag_days'].mean():.3f} days")
    print(f"Mean longitude shift: {j['ballistic_shift_deg'].mean():.3f} deg")

if len(js) == 0:
    raise RuntimeError("No ballistic matched rows were created.")

j = pd.concat(js, ignore_index=True)
l = pd.DataFrame(ls).sort_values("cr")

jp = cd / "pfss_omni_ballistic_matched_rows.csv"
lp = cd / "pfss_omni_ballistic_lag_summary.csv"
lt = td / "pfss_omni_ballistic_lag_summary_compact.csv"

j.to_csv(jp, index=False)
l.to_csv(lp, index=False)
l.to_csv(lt, index=False)

g = j.groupby(["cr", "ballistic_phase10_deg"], as_index=False).agg(
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
    lag_mean_days=("ballistic_lag_days", "mean"),
    shift_mean_deg=("ballistic_shift_deg", "mean"),
    pfss_equator_abs_mean=("equator_abs_br", "mean"),
    pfss_equator_signed_mean=("equator_signed_br", "mean"),
    pfss_midlat_abs_mean=("midlat_abs_br", "mean"),
    pfss_global_abs_mean=("global_abs_br", "mean"),
    pfss_global_signed_mean=("global_signed_br", "mean"),
)

gp = cd / "pfss_omni_ballistic_phase10_binned.csv"
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

def cor(d, x, y, md):
    q = d[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()

    if len(q) < 3:
        return np.nan, len(q)

    if q[x].std() == 0 or q[y].std() == 0:
        return np.nan, len(q)

    return float(q[x].corr(q[y], method=md)), len(q)

rs = []

for x in xs:
    for y in ys:
        rr, n = cor(g, x, y, "pearson")
        rs.append({"cr": "all", "x": x, "y": y, "method": "pearson", "r": rr, "n": n})

        rr, n = cor(g, x, y, "spearman")
        rs.append({"cr": "all", "x": x, "y": y, "method": "spearman", "r": rr, "n": n})

for cr in sorted(g["cr"].dropna().astype(int).unique()):
    q = g[g["cr"] == cr]

    for x in xs:
        for y in ys:
            rr, n = cor(q, x, y, "pearson")
            rs.append({"cr": cr, "x": x, "y": y, "method": "pearson", "r": rr, "n": n})

            rr, n = cor(q, x, y, "spearman")
            rs.append({"cr": cr, "x": x, "y": y, "method": "spearman", "r": rr, "n": n})

r = pd.DataFrame(rs)

rp = cd / "pfss_omni_ballistic_correlation_summary.csv"
rt = td / "pfss_omni_ballistic_correlation_summary_compact.csv"

r.to_csv(rp, index=False)
r.to_csv(rt, index=False)

for x, y, fn, lab in [
    ("pfss_equator_abs_mean", "speed_mean", "ballistic_eqabs_vs_omni_speed.png", "OMNI speed, km/s"),
    ("pfss_equator_abs_mean", "bmag_mean", "ballistic_eqabs_vs_omni_bmag.png", "OMNI |B|, nT"),
    ("pfss_global_abs_mean", "speed_mean", "ballistic_globalabs_vs_omni_speed.png", "OMNI speed, km/s"),
    ("pfss_global_abs_mean", "bmag_mean", "ballistic_globalabs_vs_omni_bmag.png", "OMNI |B|, nT"),
]:
    plt.figure(figsize=(7, 5))

    for cr in sorted(g["cr"].dropna().astype(int).unique()):
        q = g[g["cr"] == cr]
        plt.scatter(q[x], q[y], label=f"CR {cr}", s=18)

    plt.xlabel(x)
    plt.ylabel(lab)
    plt.title(fn.replace(".png", ""))
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(fd / fn, dpi=200)
    plt.close()

for y, lab, fn in [
    ("speed_mean", "OMNI speed, km/s", "ballistic_phase10_omni_speed_by_cr.png"),
    ("bmag_mean", "OMNI |B|, nT", "ballistic_phase10_omni_bmag_by_cr.png"),
    ("pfss_equator_abs_mean", "PFSS equator unsigned Br", "ballistic_phase10_pfss_eqabs_by_cr.png"),
    ("lag_mean_days", "Ballistic lag, days", "ballistic_phase10_lag_by_cr.png"),
]:
    plt.figure(figsize=(10, 5))

    for cr in sorted(g["cr"].dropna().astype(int).unique()):
        q = g[g["cr"] == cr]
        plt.plot(q["ballistic_phase10_deg"], q[y], label=f"CR {cr}")

    plt.xlabel("Ballistic-shifted source-surface phase, degrees")
    plt.ylabel(lab)
    plt.title(fn.replace(".png", ""))
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(fd / fn, dpi=200)
    plt.close()

m.to_csv(mf, index=False)

print()
print("Ballistic match summary:")
print(f"Matched rows: {len(j)}")
print(f"Unmatched rows: {int(j['equator_abs_br'].isna().sum())}")
print(f"10 degree bins: {len(g)}")
print(f"CRs: {sorted(j['cr'].dropna().astype(int).unique().tolist())}")

print()
print("Lag summary:")
print(l.to_string(index=False))

print()
print("Best all-CR absolute ballistic correlations:")
q = r[(r["cr"].astype(str) == "all") & (r["r"].notna())].copy()
q["abs_r"] = q["r"].abs()
q = q.sort_values("abs_r", ascending=False).head(12)
print(q[["x", "y", "method", "r", "n"]].to_string(index=False))

print()
print(f"Saved matched rows: {jp}")
print(f"Saved lag summary: {lp}")
print(f"Saved 10 degree binned table: {gp}")
print(f"Saved correlation table: {rp}")
print(f"Saved compact correlation table: {rt}")
print(f"Saved figures: {fd}")

print()
print("Manifest ballistic match status:")
print(m[["cr", "proxy_status", "omni_feature_status", "phase_match_status", "ballistic_match_status"]].to_string(index=False))

print()
print("Status: pass")
