from pathlib import Path
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
