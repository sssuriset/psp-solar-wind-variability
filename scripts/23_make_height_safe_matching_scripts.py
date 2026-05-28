from pathlib import Path

scripts = Path("scripts")
heights = [2.0, 2.5, 3.0]

def ftag(r):
    return f"rss{r:.1f}"

def stag(r):
    return f"rss{r:.1f}".replace(".", "p")

def must_replace(text, old, new, label):
    if old not in text:
        raise RuntimeError(f"Could not find expected text for {label}:\n{old}")
    return text.replace(old, new)

for r in heights:
    rf = ftag(r)
    sf = stag(r)

    # 17: ballistic match
    t17 = (scripts / "17_build_pfss_omni_ballistic_match.py").read_text()
    t17 = must_replace(t17, "rss = 2.5", f"rss = {r:.1f}", "17 rss")
    t17 = must_replace(
        t17,
        "pfss_hmi_selected_crs_rss2.5_longitude_profiles_all.csv",
        f"pfss_hmi_selected_crs_{rf}_longitude_profiles_all.csv",
        "17 PFSS input",
    )
    t17 = must_replace(
        t17,
        'cd = b / "data" / "processed" / "comparison"',
        f'cd = b / "data" / "processed" / "comparison" / "{rf}"',
        "17 comparison folder",
    )
    t17 = must_replace(
        t17,
        'fd = b / "outputs" / "figures" / "comparison"',
        f'fd = b / "outputs" / "figures" / "comparison" / "{rf}"',
        "17 figure folder",
    )
    t17 = must_replace(
        t17,
        'td = b / "outputs" / "tables"',
        f'td = b / "outputs" / "tables" / "{rf}"',
        "17 table folder",
    )
    (scripts / f"17_build_pfss_omni_ballistic_match_{sf}.py").write_text(t17)

    # 18: fixed-lag scan
    t18 = (scripts / "18_scan_pfss_omni_fixed_lags.py").read_text()
    t18 = must_replace(
        t18,
        "pfss_hmi_selected_crs_rss2.5_longitude_profiles_all.csv",
        f"pfss_hmi_selected_crs_{rf}_longitude_profiles_all.csv",
        "18 PFSS input",
    )
    t18 = must_replace(
        t18,
        'bf = b / "data" / "processed" / "comparison" / "pfss_omni_ballistic_correlation_summary.csv"',
        f'bf = b / "data" / "processed" / "comparison" / "{rf}" / "pfss_omni_ballistic_correlation_summary.csv"',
        "18 ballistic correlation input",
    )
    t18 = must_replace(
        t18,
        'cd = b / "data" / "processed" / "comparison"',
        f'cd = b / "data" / "processed" / "comparison" / "{rf}"',
        "18 comparison folder",
    )
    t18 = must_replace(
        t18,
        'fd = b / "outputs" / "figures" / "comparison"',
        f'fd = b / "outputs" / "figures" / "comparison" / "{rf}"',
        "18 figure folder",
    )
    t18 = must_replace(
        t18,
        'td = b / "outputs" / "tables"',
        f'td = b / "outputs" / "tables" / "{rf}"',
        "18 table folder",
    )
    (scripts / f"18_scan_pfss_omni_fixed_lags_{sf}.py").write_text(t18)

    # 19: null test
    t19 = (scripts / "19_null_ballistic.py").read_text()
    t19 = must_replace(
        t19,
        'gf = b / "data" / "processed" / "comparison" / "pfss_omni_ballistic_phase10_binned.csv"',
        f'gf = b / "data" / "processed" / "comparison" / "{rf}" / "pfss_omni_ballistic_phase10_binned.csv"',
        "19 ballistic binned input",
    )
    t19 = must_replace(
        t19,
        'bf = b / "data" / "processed" / "comparison" / "pfss_omni_ballistic_correlation_summary.csv"',
        f'bf = b / "data" / "processed" / "comparison" / "{rf}" / "pfss_omni_ballistic_correlation_summary.csv"',
        "19 ballistic correlation input",
    )
    t19 = must_replace(
        t19,
        'pf = b / "data" / "processed" / "comparison" / "pfss_omni_phase_correlation_summary.csv"',
        f'pf = b / "data" / "processed" / "comparison" / "{rf}" / "pfss_omni_phase_correlation_summary.csv"',
        "19 phase correlation input",
    )
    t19 = must_replace(
        t19,
        'ff = b / "data" / "processed" / "comparison" / "pfss_omni_fixed_lag_scan_correlations.csv"',
        f'ff = b / "data" / "processed" / "comparison" / "{rf}" / "pfss_omni_fixed_lag_scan_correlations.csv"',
        "19 fixed lag input",
    )
    t19 = must_replace(
        t19,
        'cd = b / "data" / "processed" / "comparison"',
        f'cd = b / "data" / "processed" / "comparison" / "{rf}"',
        "19 comparison folder",
    )
    t19 = must_replace(
        t19,
        'fd = b / "outputs" / "figures" / "comparison"',
        f'fd = b / "outputs" / "figures" / "comparison" / "{rf}"',
        "19 figure folder",
    )
    t19 = must_replace(
        t19,
        'td = b / "outputs" / "tables"',
        f'td = b / "outputs" / "tables" / "{rf}"',
        "19 table folder",
    )
    t19 = must_replace(
        t19,
        "pr = pd.read_csv(pf)",
        'if pf.exists():\n    pr = pd.read_csv(pf)\nelse:\n    pr = pd.DataFrame(columns=["cr", "x", "y", "method", "r", "n"])',
        "19 optional phase read",
    )
    (scripts / f"19_null_ballistic_{sf}.py").write_text(t19)

summary = r'''from pathlib import Path
import numpy as np
import pandas as pd

b = Path.cwd()
td = b / "outputs" / "tables"
rd = b / "outputs" / "results"
rd.mkdir(parents=True, exist_ok=True)

heights = ["rss2.0", "rss2.5", "rss3.0"]
rows = []

for h in heights:
    cf = td / h / "pfss_omni_ballistic_correlation_summary_compact.csv"
    nf = td / h / "pfss_omni_ballistic_shift_null_summary_compact.csv"

    if not cf.exists():
        print(f"Missing ballistic correlation table for {h}: {cf}")
        continue

    c = pd.read_csv(cf)
    c = c[(c["cr"].astype(str) == "all") & (c["r"].notna())].copy()
    c["rss"] = float(h.replace("rss", ""))
    c["abs_r"] = c["r"].abs()

    if nf.exists():
        n = pd.read_csv(nf)
        keep = ["x", "y", "method", "p_two_sided_abs", "abs_percentile", "null_abs_p95"]
        c = c.merge(n[keep], on=["x", "y", "method"], how="left")
    else:
        c["p_two_sided_abs"] = np.nan
        c["abs_percentile"] = np.nan
        c["null_abs_p95"] = np.nan

    rows.append(c)

if len(rows) == 0:
    raise RuntimeError("No height-specific ballistic tables were found.")

s = pd.concat(rows, ignore_index=True)

s = s.sort_values(["rss", "abs_r"], ascending=[True, False])
op = td / "pfss_omni_height_sensitivity_summary.csv"
s.to_csv(op, index=False)

target = s[
    (s["x"] == "pfss_equator_abs_mean") &
    (s["y"] == "bmag_mean") &
    (s["method"] == "spearman")
].copy()
target = target.sort_values("rss")

tp = td / "pfss_omni_height_sensitivity_eqabs_bmag_spearman.csv"
target.to_csv(tp, index=False)

best = s.sort_values("abs_r", ascending=False).head(15)

txt = []
txt.append("PFSS source-surface-height sensitivity summary")
txt.append("")
txt.append("Target poster relation:")
txt.append(target[["rss", "x", "y", "method", "r", "n", "p_two_sided_abs", "abs_percentile"]].to_string(index=False))
txt.append("")
txt.append("Best all-height ballistic correlations:")
txt.append(best[["rss", "x", "y", "method", "r", "n", "p_two_sided_abs", "abs_percentile"]].to_string(index=False))

rp = rd / "pfss_omni_height_sensitivity_summary.txt"
rp.write_text("\n".join(txt))

print("Height sensitivity target relation:")
print(target[["rss", "x", "y", "method", "r", "n", "p_two_sided_abs", "abs_percentile"]].to_string(index=False))
print()
print("Best all-height ballistic correlations:")
print(best[["rss", "x", "y", "method", "r", "n", "p_two_sided_abs", "abs_percentile"]].to_string(index=False))
print()
print(f"Saved table: {op}")
print(f"Saved target table: {tp}")
print(f"Saved text summary: {rp}")
print()
print("Status: pass")
'''

(scripts / "20b_make_height_sensitivity_summary.py").write_text(summary)

print("Created height-safe matching scripts for 17, 18, 19.")
print("Created scripts/20b_make_height_sensitivity_summary.py")
