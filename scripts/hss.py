from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

b = Path.cwd()

heights = [2.0, 2.5, 3.0]

outtab = b / "outputs" / "tables"
outfig = b / "outputs" / "figures" / "diag"
outres = b / "outputs" / "results"

outtab.mkdir(parents=True, exist_ok=True)
outfig.mkdir(parents=True, exist_ok=True)
outres.mkdir(parents=True, exist_ok=True)

def circular_signed_diff_deg(a, b):
    # signed difference a - b on [-180, 180)
    return ((a - b + 180.0) % 360.0) - 180.0

def find_local_peaks(q, phase_col, value_col, threshold):
    q = q.sort_values(phase_col).reset_index(drop=True).copy()
    vals = q[value_col].to_numpy(dtype=float)
    phases = q[phase_col].to_numpy(dtype=float)

    if len(q) == 0:
        return pd.DataFrame(columns=[phase_col, value_col])

    peaks = []
    n = len(q)

    for i in range(n):
        v = vals[i]
        if not np.isfinite(v):
            continue

        left = vals[(i - 1) % n]
        right = vals[(i + 1) % n]

        if v >= threshold and v >= left and v >= right:
            peaks.append(i)

    # Fallback: keep the strongest point if threshold/local-max rule finds nothing.
    if len(peaks) == 0:
        finite = np.where(np.isfinite(vals))[0]
        if len(finite) == 0:
            return pd.DataFrame(columns=[phase_col, value_col])
        peaks = [finite[np.nanargmax(vals[finite])]]

    p = q.loc[peaks, [phase_col, value_col]].copy()
    return p.sort_values(value_col, ascending=False).reset_index(drop=True)

detail_rows = []
summary_rows = []

for rss in heights:
    h = f"rss{rss:.1f}"

    gf = b / "data" / "processed" / "comparison" / h / "pfss_omni_ballistic_phase10_binned.csv"
    lf = b / "data" / "processed" / "comparison" / h / "pfss_omni_ballistic_lag_summary.csv"

    if not gf.exists():
        raise RuntimeError(f"Missing binned table for {h}: {gf}")

    g = pd.read_csv(gf)
    l = pd.read_csv(lf) if lf.exists() else pd.DataFrame()

    need = ["cr", "ballistic_phase10_deg", "speed_mean", "pfss_equator_abs_mean"]
    miss = [c for c in need if c not in g.columns]
    if miss:
        raise RuntimeError(f"{h} missing columns: {miss}")

    for c in need:
        g[c] = pd.to_numeric(g[c], errors="coerce")

    for cr in sorted(g["cr"].dropna().astype(int).unique()):
        q = g[g["cr"].astype(int) == cr].copy()
        q = q.dropna(subset=["ballistic_phase10_deg", "speed_mean", "pfss_equator_abs_mean"])

        if len(q) == 0:
            continue

        if len(l) > 0 and "period_days" in l.columns:
            lr = l[l["cr"].astype(int) == cr]
            period_days = float(lr["period_days"].iloc[0]) if len(lr) else 27.2753
        else:
            period_days = 27.2753

        speed_threshold = max(450.0, float(q["speed_mean"].quantile(0.75)))
        pfss_threshold = float(q["pfss_equator_abs_mean"].quantile(0.75))

        speed_peaks = find_local_peaks(q, "ballistic_phase10_deg", "speed_mean", speed_threshold)
        pfss_peaks = find_local_peaks(q, "ballistic_phase10_deg", "pfss_equator_abs_mean", pfss_threshold)

        for _, sp in speed_peaks.iterrows():
            speed_phase = float(sp["ballistic_phase10_deg"])

            if len(pfss_peaks) == 0:
                continue

            diffs = circular_signed_diff_deg(speed_phase, pfss_peaks["ballistic_phase10_deg"].to_numpy(dtype=float))
            k = int(np.argmin(np.abs(diffs)))

            pfss_phase = float(pfss_peaks.loc[k, "ballistic_phase10_deg"])
            signed_err_deg = float(diffs[k])
            abs_err_deg = abs(signed_err_deg)
            signed_err_days = signed_err_deg / 360.0 * period_days
            abs_err_days = abs_err_deg / 360.0 * period_days

            detail_rows.append({
                "rss": rss,
                "cr": cr,
                "period_days": period_days,
                "speed_peak_phase_deg": speed_phase,
                "speed_peak_km_s": float(sp["speed_mean"]),
                "nearest_pfss_peak_phase_deg": pfss_phase,
                "nearest_pfss_peak_value": float(pfss_peaks.loc[k, "pfss_equator_abs_mean"]),
                "signed_phase_error_deg": signed_err_deg,
                "abs_phase_error_deg": abs_err_deg,
                "signed_timing_error_days": signed_err_days,
                "abs_timing_error_days": abs_err_days,
                "speed_threshold_km_s": speed_threshold,
                "pfss_threshold": pfss_threshold,
                "speed_peak_count": len(speed_peaks),
                "pfss_peak_count": len(pfss_peaks),
            })

detail = pd.DataFrame(detail_rows)

if len(detail) == 0:
    raise RuntimeError("No high-speed timing peak matches were created.")

summary = detail.groupby("rss", as_index=False).agg(
    matched_stream_peaks=("abs_timing_error_days", "size"),
    median_abs_timing_error_days=("abs_timing_error_days", "median"),
    mean_abs_timing_error_days=("abs_timing_error_days", "mean"),
    median_abs_phase_error_deg=("abs_phase_error_deg", "median"),
    mean_abs_phase_error_deg=("abs_phase_error_deg", "mean"),
)

by_cr = detail.groupby(["rss", "cr"], as_index=False).agg(
    matched_stream_peaks=("abs_timing_error_days", "size"),
    median_abs_timing_error_days=("abs_timing_error_days", "median"),
    mean_abs_timing_error_days=("abs_timing_error_days", "mean"),
    median_abs_phase_error_deg=("abs_phase_error_deg", "median"),
)

detail_path = outtab / "pfss_omni_high_speed_stream_timing_detail.csv"
summary_path = outtab / "pfss_omni_high_speed_stream_timing_summary.csv"
by_cr_path = outtab / "pfss_omni_high_speed_stream_timing_by_cr.csv"

detail.to_csv(detail_path, index=False)
summary.to_csv(summary_path, index=False)
by_cr.to_csv(by_cr_path, index=False)

plt.figure(figsize=(7, 5))
plt.plot(summary["rss"], summary["median_abs_timing_error_days"], marker="o", linewidth=2)
plt.xlabel("PFSS source-surface height, Rs")
plt.ylabel("Median absolute timing error, days")
plt.title("High-speed stream timing diagnostic")
plt.xticks(summary["rss"])
plt.tight_layout()
fig_path = outfig / "fig_high_speed_stream_timing_by_height.png"
plt.savefig(fig_path, dpi=300)
plt.close()

txt = []
txt.append("High-speed stream timing diagnostic")
txt.append("")
txt.append("Definition:")
txt.append("For each CR and source-surface height, OMNI speed peaks were compared with the nearest PFSS equatorial unsigned Br peak after speed-dependent ballistic alignment.")
txt.append("This is a peak-alignment diagnostic, not a calibrated solar-wind speed prediction.")
txt.append("")
txt.append("Summary by source-surface height:")
txt.append(summary.to_string(index=False))
txt.append("")
txt.append("Summary by CR:")
txt.append(by_cr.to_string(index=False))

txt_path = outres / "pfss_omni_high_speed_stream_timing_summary.txt"
txt_path.write_text("\n".join(txt))

print("High-speed stream timing summary:")
print(summary.to_string(index=False))
print()
print("By CR:")
print(by_cr.to_string(index=False))
print()
print(f"Saved detail: {detail_path}")
print(f"Saved summary: {summary_path}")
print(f"Saved by-CR summary: {by_cr_path}")
print(f"Saved figure: {fig_path}")
print(f"Saved text summary: {txt_path}")
print()
print("Status: pass")
