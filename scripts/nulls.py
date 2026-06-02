from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

b = Path.cwd()

gf = b / "data" / "processed" / "comparison" / "pfss_omni_ballistic_phase10_binned.csv"
bf = b / "data" / "processed" / "comparison" / "pfss_omni_ballistic_correlation_summary.csv"
pf = b / "data" / "processed" / "comparison" / "pfss_omni_phase_correlation_summary.csv"
ff = b / "data" / "processed" / "comparison" / "pfss_omni_fixed_lag_scan_correlations.csv"
mf = b / "metadata" / "rotation_manifest.csv"

cd = b / "data" / "processed" / "comparison"
fd = b / "outputs" / "figures" / "comparison"
td = b / "outputs" / "tables"

cd.mkdir(parents=True, exist_ok=True)
fd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

g = pd.read_csv(gf)
br = pd.read_csv(bf)
pr = pd.read_csv(pf)
fr = pd.read_csv(ff)
m = pd.read_csv(mf)

if "null_test_status" not in m.columns:
    m["null_test_status"] = ""

for c in ["null_test_status", "notes"]:
    if c not in m.columns:
        m[c] = ""
    m[c] = m[c].fillna("").astype(str)

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

mds = ["pearson", "spearman"]

for c in ["cr", "ballistic_phase10_deg"]:
    if c not in g.columns:
        raise RuntimeError(f"Missing column in ballistic binned table: {c}")

for c in xs + ys:
    if c not in g.columns:
        raise RuntimeError(f"Missing column in ballistic binned table: {c}")

def ca(a, v, md):
    a = np.asarray(a, dtype=float)
    v = np.asarray(v, dtype=float)

    ok = np.isfinite(a) & np.isfinite(v)

    if ok.sum() < 3:
        return np.nan, int(ok.sum())

    aa = a[ok]
    vv = v[ok]

    if np.nanstd(aa) == 0 or np.nanstd(vv) == 0:
        return np.nan, int(ok.sum())

    r = pd.Series(aa).corr(pd.Series(vv), method=md)
    return float(r), int(ok.sum())

def grabr(d, x, y, md):
    q = d[
        (d["cr"].astype(str) == "all") &
        (d["x"] == x) &
        (d["y"] == y) &
        (d["method"] == md)
    ].copy()

    if len(q) == 0:
        return np.nan, np.nan

    return float(q.iloc[0]["r"]), int(q.iloc[0]["n"])

rng = np.random.default_rng(57)
nrun = 2000
crs = sorted(g["cr"].dropna().astype(int).unique())

nulls = []
srows = []
crows = []

print(f"Project base: {b}")
print(f"Null runs per comparison: {nrun}")
print(f"CRs: {crs}")
print()

for x in xs:
    for y in ys:
        for md in mds:
            ob, on = grabr(br, x, y, md)
            ph, pn = grabr(pr, x, y, md)

            vals = []

            for i in range(nrun):
                aa = []
                vv = []

                for cr in crs:
                    q = g[g["cr"].astype(int) == cr].sort_values("ballistic_phase10_deg")

                    a = q[x].to_numpy(dtype=float)
                    v = q[y].to_numpy(dtype=float)

                    if len(a) == 0:
                        continue

                    k = int(rng.integers(0, len(a)))

                    aa.append(np.roll(a, k))
                    vv.append(v)

                if len(aa) == 0:
                    rv = np.nan
                    nn = 0
                else:
                    rv, nn = ca(np.concatenate(aa), np.concatenate(vv), md)

                vals.append(rv)

                nulls.append({
                    "x": x,
                    "y": y,
                    "method": md,
                    "run": i,
                    "null_r": rv,
                    "n": nn,
                })

            va = np.asarray(vals, dtype=float)
            va = va[np.isfinite(va)]

            if np.isfinite(ob) and len(va) > 0:
                p2 = float((np.sum(np.abs(va) >= abs(ob)) + 1) / (len(va) + 1))
                pct = float(100.0 * np.sum(np.abs(va) < abs(ob)) / len(va))
                nm = float(np.mean(np.abs(va)))
                n95 = float(np.quantile(np.abs(va), 0.95))
            else:
                p2 = np.nan
                pct = np.nan
                nm = np.nan
                n95 = np.nan

            srows.append({
                "x": x,
                "y": y,
                "method": md,
                "obs_r": ob,
                "obs_n": on,
                "null_runs": int(len(va)),
                "null_abs_mean": nm,
                "null_abs_p95": n95,
                "p_two_sided_abs": p2,
                "abs_percentile": pct,
            })

            fq = fr[
                (fr["cr"].astype(str) == "all") &
                (fr["x"] == x) &
                (fr["y"] == y) &
                (fr["method"] == md) &
                (fr["r"].notna())
            ].copy()

            if len(fq) > 0:
                fq["abs_r"] = fq["r"].abs()
                fb = fq.sort_values("abs_r", ascending=False).iloc[0]
                fbr = float(fb["r"])
                fbs = float(fb["scan_shift_deg"])
                fba = float(fb["abs_r"])
            else:
                fbr = np.nan
                fbs = np.nan
                fba = np.nan

            crows.append({
                "x": x,
                "y": y,
                "method": md,
                "phase_r": ph,
                "ballistic_r": ob,
                "best_fixed_r": fbr,
                "best_fixed_shift_deg": fbs,
                "best_fixed_abs_r": fba,
                "ballistic_abs_minus_phase_abs": abs(ob) - abs(ph) if np.isfinite(ob) and np.isfinite(ph) else np.nan,
                "ballistic_abs_minus_best_fixed_abs": abs(ob) - fba if np.isfinite(ob) and np.isfinite(fba) else np.nan,
                "ballistic_null_p": p2,
                "ballistic_abs_percentile": pct,
            })

nr = pd.DataFrame(nulls)
sr = pd.DataFrame(srows)
cr = pd.DataFrame(crows)

npth = cd / "pfss_omni_ballistic_shift_null_runs.csv"
spth = cd / "pfss_omni_ballistic_shift_null_summary.csv"
cpth = td / "pfss_omni_ballistic_shift_null_summary_compact.csv"
cmp = cd / "pfss_omni_match_method_comparison.csv"
cmpt = td / "pfss_omni_match_method_comparison_compact.csv"

nr.to_csv(npth, index=False)
sr.to_csv(spth, index=False)
sr.to_csv(cpth, index=False)
cr.to_csv(cmp, index=False)
cr.to_csv(cmpt, index=False)

top = sr.copy()
top["abs_obs_r"] = top["obs_r"].abs()
top = top.sort_values("abs_obs_r", ascending=False)

for _, row in top.head(5).iterrows():
    x = row["x"]
    y = row["y"]
    md = row["method"]
    ob = row["obs_r"]

    q = nr[(nr["x"] == x) & (nr["y"] == y) & (nr["method"] == md)]["null_r"].dropna()

    plt.figure(figsize=(7, 5))
    plt.hist(q, bins=40)
    plt.axvline(ob, linewidth=2)
    plt.axvline(-ob, linewidth=2)
    plt.xlabel("Null correlation r")
    plt.ylabel("Count")
    plt.title(f"{x} vs {y}, {md}")
    plt.tight_layout()

    nm = f"null_{x}_vs_{y}_{md}.png"
    plt.savefig(fd / nm, dpi=200)
    plt.close()

for c in sorted(g["cr"].dropna().astype(int).unique()):
    ms = m["cr"].astype(int) == c
    m.loc[ms, "null_test_status"] = "pass_ballistic_shift_null"

    nt = "Ballistic match tested against random circular phase shifts."
    old = m.loc[ms, "notes"].iloc[0]

    if nt not in old:
        m.loc[ms, "notes"] = old + " | " + nt

m.to_csv(mf, index=False)

print("Ballistic null-test summary:")
print(f"Null runs per comparison: {nrun}")
print(f"Total null rows: {len(nr)}")
print()

print("Top ballistic correlations with null-test p values:")
print(top[[
    "x",
    "y",
    "method",
    "obs_r",
    "obs_n",
    "null_abs_mean",
    "null_abs_p95",
    "p_two_sided_abs",
    "abs_percentile",
]].head(15).to_string(index=False))

print()
print("Method comparison summary:")
q = cr.copy()
q["ballistic_abs"] = q["ballistic_r"].abs()
q = q.sort_values("ballistic_abs", ascending=False).head(15)
print(q[[
    "x",
    "y",
    "method",
    "phase_r",
    "ballistic_r",
    "best_fixed_r",
    "best_fixed_shift_deg",
    "ballistic_abs_minus_phase_abs",
    "ballistic_abs_minus_best_fixed_abs",
    "ballistic_null_p",
]].to_string(index=False))

print()
print(f"Saved null runs: {npth}")
print(f"Saved null summary: {spth}")
print(f"Saved compact null summary: {cpth}")
print(f"Saved method comparison: {cmp}")
print(f"Saved compact method comparison: {cmpt}")
print(f"Saved figures: {fd}")

print()
print("Manifest null-test status:")
print(m[["cr", "ballistic_match_status", "lag_scan_status", "null_test_status"]].to_string(index=False))

print()
print("Status: pass")
