from pathlib import Path
import math
import numpy as np
import pandas as pd

BASE = Path.cwd()

FINAL_CRS = [2281, 2283, 2284, 2286, 2290]
AMBIENT_CRS = [2281, 2283, 2286, 2290]
DEBUG_ONLY_CRS = [2287]

HEIGHTS = {
    2.0: BASE / "data" / "processed" / "comparison" / "rss2.0" / "pfss_omni_ballistic_matched_rows.csv",
    2.5: BASE / "data" / "processed" / "comparison" / "rss2.5" / "pfss_omni_ballistic_matched_rows.csv",
    3.0: BASE / "data" / "processed" / "comparison" / "rss3.0" / "pfss_omni_ballistic_matched_rows.csv",
}

OUT_TABLE = BASE / "outputs" / "tables" / "final"
OUT_RESULT = BASE / "outputs" / "results" / "final"
OUT_DATA = BASE / "data" / "processed" / "comparison" / "final"

OUT_TABLE.mkdir(parents=True, exist_ok=True)
OUT_RESULT.mkdir(parents=True, exist_ok=True)
OUT_DATA.mkdir(parents=True, exist_ok=True)

XS = [
    "pfss_equator_abs_mean",
    "pfss_midlat_abs_mean",
    "pfss_global_abs_mean",
]

YS = [
    "speed_mean",
    "bmag_mean",
    "density_mean",
    "pdyn_mean",
    "ma_mean",
]

TARGET_X = "pfss_equator_abs_mean"
TARGET_Y = "bmag_mean"
TARGET_METHOD = "spearman"


def corr(df, x, y, method):
    q = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(q) < 3:
        return np.nan, int(len(q))
    if q[x].std() == 0 or q[y].std() == 0:
        return np.nan, int(len(q))
    return float(q[x].corr(q[y], method=method)), int(len(q))


def bin10(matched):
    g = matched.groupby(["cr", "ballistic_phase10_deg"], as_index=False).agg(
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
    return g


def corrtab(g, rss, sample_name):
    rows = []

    for x in XS:
        for y in YS:
            for method in ["pearson", "spearman"]:
                r, n = corr(g, x, y, method)
                rows.append({
                    "sample": sample_name,
                    "rss": rss,
                    "cr": "all",
                    "x": x,
                    "y": y,
                    "method": method,
                    "r": r,
                    "n": n,
                    "abs_r": abs(r) if np.isfinite(r) else np.nan,
                })

    for cr in sorted(g["cr"].dropna().astype(int).unique()):
        q = g[g["cr"].astype(int) == cr]
        for x in XS:
            for y in YS:
                for method in ["pearson", "spearman"]:
                    r, n = corr(q, x, y, method)
                    rows.append({
                        "sample": sample_name,
                        "rss": rss,
                        "cr": cr,
                        "x": x,
                        "y": y,
                        "method": method,
                        "r": r,
                        "n": n,
                        "abs_r": abs(r) if np.isfinite(r) else np.nan,
                    })

    return pd.DataFrame(rows)


def shiftnull(g, x, y, method, n_iter=2000, seed=57):
    obs, n = corr(g, x, y, method)
    if not np.isfinite(obs):
        return np.nan, np.nan, np.nan, obs, n

    rng = np.random.default_rng(seed)
    base = g[["cr", "ballistic_phase10_deg", x, y]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    base = base.sort_values(["cr", "ballistic_phase10_deg"])

    nulls = []
    for _ in range(n_iter):
        parts = []
        for cr, q in base.groupby("cr"):
            qq = q.copy()
            vals = qq[x].to_numpy()
            if len(vals) > 1:
                shift = int(rng.integers(0, len(vals)))
                qq[x] = np.roll(vals, shift)
            parts.append(qq)
        shifted = pd.concat(parts, ignore_index=True)
        rr, _ = corr(shifted, x, y, method)
        if np.isfinite(rr):
            nulls.append(rr)

    nulls = np.asarray(nulls, dtype=float)
    if len(nulls) == 0:
        return np.nan, np.nan, np.nan, obs, n

    p = (np.sum(np.abs(nulls) >= abs(obs)) + 1.0) / (len(nulls) + 1.0)
    percentile = 100.0 * np.mean(np.abs(nulls) <= abs(obs))
    p95 = float(np.percentile(np.abs(nulls), 95))
    return float(p), float(percentile), p95, obs, n


def addnull(corr_df, g, rss, sample_name):
    rows = []
    q = corr_df[(corr_df["sample"] == sample_name) & (corr_df["rss"] == rss) & (corr_df["cr"].astype(str) == "all")].copy()

    for _, row in q.iterrows():
        p, pct, p95, obs, n = shiftnull(
            g,
            row["x"],
            row["y"],
            row["method"],
            n_iter=2000,
            seed=5700 + int(rss * 10),
        )
        r = dict(row)
        r["null_p_two_sided_abs"] = p
        r["null_abs_percentile"] = pct
        r["null_abs_p95"] = p95
        rows.append(r)

    return pd.DataFrame(rows)


def runsample(rss, matched_all, crs, sample_name):
    matched = matched_all[matched_all["cr"].astype(int).isin(crs)].copy()
    g = bin10(matched)

    sample_dir = OUT_DATA / f"rss{rss:.1f}" / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)

    matched.to_csv(sample_dir / "pfss_omni_ballistic_matched_rows.csv", index=False)
    g.to_csv(sample_dir / "pfss_omni_ballistic_phase10_binned.csv", index=False)

    corr = corrtab(g, rss, sample_name)
    null = addnull(corr, g, rss, sample_name)

    return matched, g, corr, null


all_corr = []
all_null = []
all_lag = []
all_polarity = []

for rss, path in HEIGHTS.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing matched-row table for rss={rss}: {path}")

    d = pd.read_csv(path)
    d["cr"] = pd.to_numeric(d["cr"], errors="coerce").astype(int)

    present = sorted(d["cr"].dropna().astype(int).unique().tolist())
    print(f"rss={rss}: input CRs = {present}")

    final_matched, final_binned, corr, null = runsample(rss, d, FINAL_CRS, "with2284")
    ambient_matched, ambient_binned, corr_amb, null_amb = runsample(rss, d, AMBIENT_CRS, "ambient_no2284")

    all_corr.extend([corr, corr_amb])
    all_null.extend([null, null_amb])

    lag = final_matched.groupby("cr", as_index=False).agg(
        rows=("time", "size"),
        lag_mean_days=("ballistic_lag_days", "mean"),
        lag_min_days=("ballistic_lag_days", "min"),
        lag_max_days=("ballistic_lag_days", "max"),
        shift_mean_deg=("ballistic_shift_deg", "mean"),
        shift_min_deg=("ballistic_shift_deg", "min"),
        shift_max_deg=("ballistic_shift_deg", "max"),
    )
    lag["rss"] = rss
    all_lag.append(lag)

    # RTN polarity agreement, using the same GSE -> RTN approximation already used in the project.
    pol = final_matched.copy()
    pol["br_rtn_nt"] = -pd.to_numeric(pol["bx_gse_nt"], errors="coerce")
    pol["ballistic_phase10_deg"] = pd.to_numeric(pol["ballistic_phase10_deg"], errors="coerce")
    pg = pol.groupby(["cr", "ballistic_phase10_deg"], as_index=False).agg(
        br_rtn_mean=("br_rtn_nt", "mean"),
        pfss_equator_signed_mean=("equator_signed_br", "mean"),
        pfss_global_signed_mean=("global_signed_br", "mean"),
    )

    for proxy_col, label in [
        ("pfss_equator_signed_mean", "equator_signed_br"),
        ("pfss_global_signed_mean", "global_signed_br"),
    ]:
        q = pg[[proxy_col, "br_rtn_mean"]].dropna().copy()
        q["pfss_sign"] = np.sign(q[proxy_col])
        q["omni_sign"] = np.sign(q["br_rtn_mean"])
        q = q[(q["pfss_sign"] != 0) & (q["omni_sign"] != 0)]

        direct = int((q["pfss_sign"] == q["omni_sign"]).sum())
        inverted = int((-q["pfss_sign"] == q["omni_sign"]).sum())
        n = int(len(q))

        all_polarity.append({
            "sample": "with2284",
            "rss": rss,
            "pfss_polarity_proxy": label,
            "omni_polarity": "RTN Br sign",
            "bins": n,
            "direct_agreement_fraction": direct / n if n else np.nan,
            "direct_agreement_count": direct,
            "inverted_agreement_fraction": inverted / n if n else np.nan,
            "inverted_agreement_count": inverted,
            "best_agreement_fraction": max(direct, inverted) / n if n else np.nan,
            "best_orientation": "direct" if direct >= inverted else "inverted",
        })

corr_all = pd.concat(all_corr, ignore_index=True)
null_all = pd.concat(all_null, ignore_index=True)
lag_all = pd.concat(all_lag, ignore_index=True)
polarity_all = pd.DataFrame(all_polarity)

corr_all.to_csv(OUT_TABLE / "allcorr.csv", index=False)
null_all.to_csv(OUT_TABLE / "corr.csv", index=False)
lag_all.to_csv(OUT_TABLE / "lag.csv", index=False)
polarity_all.to_csv(OUT_TABLE / "pol.csv", index=False)

target = null_all[
    (null_all["x"] == TARGET_X)
    & (null_all["y"] == TARGET_Y)
    & (null_all["method"] == TARGET_METHOD)
].copy()

target = target.sort_values(["sample", "rss"])
target.to_csv(OUT_TABLE / "height.csv", index=False)

# Leave-one-rotation-out for target relation.
loo_rows = []
for drop_cr in FINAL_CRS:
    keep = [cr for cr in FINAL_CRS if cr != drop_cr]
    for rss, path in HEIGHTS.items():
        d = pd.read_csv(path)
        d["cr"] = pd.to_numeric(d["cr"], errors="coerce").astype(int)
        matched = d[d["cr"].isin(keep)].copy()
        g = bin10(matched)
        r, n = corr(g, TARGET_X, TARGET_Y, TARGET_METHOD)
        loo_rows.append({
            "dropped_cr": drop_cr,
            "kept_crs": " ".join(str(x) for x in keep),
            "rss": rss,
            "x": TARGET_X,
            "y": TARGET_Y,
            "method": TARGET_METHOD,
            "r": r,
            "n": n,
            "abs_r": abs(r) if np.isfinite(r) else np.nan,
        })

loo = pd.DataFrame(loo_rows)
loo["is_best_abs_r_for_dropout"] = False
for drop_cr, q in loo.groupby("dropped_cr"):
    idx = q["abs_r"].idxmax()
    loo.loc[idx, "is_best_abs_r_for_dropout"] = True

loo.to_csv(OUT_TABLE / "loo.csv", index=False)

headline = target[target["sample"] == "with2284"].sort_values(
    ["null_p_two_sided_abs", "abs_r"],
    ascending=[True, False],
).copy()

ambient_compare = target[target["sample"].isin(["with2284", "ambient_no2284"])].copy()

summary_lines = []
summary_lines.append("Final science sample locked")
summary_lines.append(f"Included CRs: {FINAL_CRS}")
summary_lines.append(f"Debug-only CRs excluded from final science metrics: {DEBUG_ONLY_CRS}")
summary_lines.append("")
summary_lines.append("Target relation: pfss_equator_abs_mean vs bmag_mean, Spearman")
summary_lines.append("")
summary_lines.append(ambient_compare[[
    "sample",
    "rss",
    "r",
    "n",
    "null_p_two_sided_abs",
    "null_abs_percentile",
    "null_abs_p95",
]].to_string(index=False))
summary_lines.append("")
summary_lines.append("Leave-one-rotation-out target relation")
summary_lines.append("")
summary_lines.append(loo[[
    "dropped_cr",
    "rss",
    "r",
    "n",
    "abs_r",
    "is_best_abs_r_for_dropout",
]].to_string(index=False))
summary_lines.append("")
summary_lines.append("RTN IMF polarity agreement")
summary_lines.append("")
summary_lines.append(polarity_all.to_string(index=False))

(OUT_RESULT / "summary.txt").write_text("\n".join(summary_lines))

print("")
print("Final sample target relation:")
print(target[[
    "sample",
    "rss",
    "r",
    "n",
    "null_p_two_sided_abs",
    "null_abs_percentile",
    "null_abs_p95",
]].to_string(index=False))

print("")
print("Leave-one-rotation-out target relation:")
print(loo[[
    "dropped_cr",
    "rss",
    "r",
    "n",
    "abs_r",
    "is_best_abs_r_for_dropout",
]].to_string(index=False))

print("")
print("Saved:")
print(OUT_TABLE / "height.csv")
print(OUT_TABLE / "loo.csv")
print(OUT_TABLE / "corr.csv")
print(OUT_TABLE / "lag.csv")
print(OUT_TABLE / "pol.csv")
print(OUT_RESULT / "summary.txt")
print("")
print("Status: pass")
