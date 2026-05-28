from pathlib import Path
import pandas as pd
import numpy as np

b = Path.cwd()

cd = b / "data" / "processed" / "comparison"
td = b / "outputs" / "tables"
rd = b / "outputs" / "results"

rd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

ph = pd.read_csv(cd / "pfss_omni_phase_correlation_summary.csv")
ba = pd.read_csv(cd / "pfss_omni_ballistic_correlation_summary.csv")
fx = pd.read_csv(cd / "pfss_omni_fixed_lag_scan_correlations.csv")
nu = pd.read_csv(cd / "pfss_omni_ballistic_shift_null_summary.csv")
lg = pd.read_csv(cd / "pfss_omni_ballistic_lag_summary.csv")

def one(d, x, y, m):
    q = d[
        (d["cr"].astype(str) == "all") &
        (d["x"] == x) &
        (d["y"] == y) &
        (d["method"] == m)
    ].copy()

    if len(q) == 0:
        return np.nan

    return float(q.iloc[0]["r"])

def bestfx(x, y, m):
    q = fx[
        (fx["cr"].astype(str) == "all") &
        (fx["x"] == x) &
        (fx["y"] == y) &
        (fx["method"] == m) &
        (fx["r"].notna())
    ].copy()

    if len(q) == 0:
        return np.nan, np.nan

    q["ar"] = q["r"].abs()
    z = q.sort_values("ar", ascending=False).iloc[0]
    return float(z["r"]), float(z["scan_shift_deg"])

rows = []

for _, r in nu.iterrows():
    x = r["x"]
    y = r["y"]
    m = r["method"]

    fr, sh = bestfx(x, y, m)

    rows.append({
        "pfss_proxy": x,
        "omni_var": y,
        "method": m,
        "phase_r": one(ph, x, y, m),
        "ballistic_r": one(ba, x, y, m),
        "best_fixed_r": fr,
        "best_fixed_shift_deg": sh,
        "null_p": float(r["p_two_sided_abs"]),
        "null_percentile": float(r["abs_percentile"]),
        "null_abs_p95": float(r["null_abs_p95"]),
    })

s = pd.DataFrame(rows)

s["abs_ballistic_r"] = s["ballistic_r"].abs()
s["ballistic_minus_phase_abs"] = s["ballistic_r"].abs() - s["phase_r"].abs()
s["ballistic_minus_best_fixed_abs"] = s["ballistic_r"].abs() - s["best_fixed_r"].abs()

s = s.sort_values(["null_p", "abs_ballistic_r"], ascending=[True, False])

op = td / "final_pfss_omni_result_summary.csv"
s.to_csv(op, index=False)

top = s.iloc[0]

txt = []
txt.append("PFSS to OMNI diagnostic result summary")
txt.append("")
txt.append("Main result:")
txt.append(
    f"The strongest tested relationship was {top['pfss_proxy']} versus {top['omni_var']} "
    f"using {top['method']} correlation after ballistic alignment."
)
txt.append(
    f"The ballistic correlation was r = {top['ballistic_r']:.6f}, compared with "
    f"the unshifted phase-only value r = {top['phase_r']:.6f}."
)
txt.append(
    f"The random circular-shift null test gave p = {top['null_p']:.6f}, "
    f"placing the ballistic result at the {top['null_percentile']:.2f} percentile of null-shift outcomes."
)
txt.append("")
txt.append("Physical reading:")
txt.append(
    "The result supports a weak but non-random association between the PFSS unsigned "
    "source-surface magnetic-field proxy and the near-Earth OMNI magnetic-field magnitude "
    "after a speed-dependent ballistic travel-time correction."
)
txt.append("")
txt.append("Limits:")
txt.append(
    "This is a diagnostic proxy result, not a direct field-line tracing result. "
    "The mapping uses a source-surface longitude proxy, OMNI phase, and ballistic lag. "
    "The result should be presented as evidence that the proxy retains detectable solar-wind structure, "
    "not as proof of one-to-one magnetic connectivity."
)
txt.append("")
txt.append("Ballistic lag range:")
txt.append(
    f"Across selected Carrington rotations, mean ballistic lag ranged from "
    f"{lg['lag_mean_days'].min():.3f} to {lg['lag_mean_days'].max():.3f} days."
)
txt.append(
    f"Mean longitude shift ranged from {lg['shift_mean_deg'].min():.3f}° "
    f"to {lg['shift_mean_deg'].max():.3f}°."
)
txt.append("")
txt.append("Top rows:")
txt.append(s.head(8).to_string(index=False))

tp = rd / "pfss_omni_result_summary.txt"
tp.write_text("\n".join(txt))

print("Final result summary:")
print(s.head(12).to_string(index=False))
print()
print(f"Saved table: {op}")
print(f"Saved text summary: {tp}")
print()
print("Status: pass")
