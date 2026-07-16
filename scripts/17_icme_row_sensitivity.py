from pathlib import Path
import numpy as np
import pandas as pd

project_root = Path(__file__).resolve().parents[1]

final_crs = [2281, 2283, 2284, 2286, 2290]
debug_crs = [2287]
rss_values = ["rss2.0", "rss2.5", "rss3.0"]

input_template = project_root / "data/processed/comparison/{rss}/pfss_omni_ballistic_matched_rows_icme_flagged.csv"
outdir = project_root / "data/processed/comparison"
outdir.mkdir(parents=True, exist_ok=True)

final_out = outdir / "final_metrics_icme_filtered.csv"
by_cr_out = outdir / "final_metrics_by_cr_icme_filtered.csv"
sensitivity_out = outdir / "icme_filter_sensitivity.csv"

def clean_bool(series):
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1", "yes", "y"])

def sign_values(series):
    x = pd.to_numeric(series, errors="coerce")
    return np.where(x > 0, 1, np.where(x < 0, -1, np.nan))

def corr_metric(df, xcol, ycol, method):
    sub = df[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 3:
        return np.nan
    if sub[xcol].nunique() < 2 or sub[ycol].nunique() < 2:
        return np.nan
    return sub[xcol].corr(sub[ycol], method=method)

def polarity_accuracy(df, pred_col):
    pred = sign_values(df[pred_col])
    obs = -sign_values(df["bx_gse_nt"])
    mask = ~pd.isna(pred) & ~pd.isna(obs)
    if mask.sum() == 0:
        return np.nan
    return float((pred[mask] == obs[mask]).mean())

def polarity_n(df, pred_col):
    pred = sign_values(df[pred_col])
    obs = -sign_values(df["bx_gse_nt"])
    mask = ~pd.isna(pred) & ~pd.isna(obs)
    return int(mask.sum())

def metric_row(df, rss, sample_name):
    return {
        "rss": rss,
        "sample": sample_name,
        "crs": " ".join(str(x) for x in sorted(df["cr"].dropna().astype(int).unique())),
        "n_rows": len(df),
        "icme_rows": int(clean_bool(df["icme_flag"]).sum()) if "icme_flag" in df.columns else 0,
        "median_speed_km_s": pd.to_numeric(df["speed_km_s"], errors="coerce").median(),
        "median_bmag_nt": pd.to_numeric(df["bmag_nt"], errors="coerce").median(),
        "spearman_speed_equator_abs_br": corr_metric(df, "equator_abs_br", "speed_km_s", "spearman"),
        "spearman_speed_midlat_abs_br": corr_metric(df, "midlat_abs_br", "speed_km_s", "spearman"),
        "spearman_speed_global_abs_br": corr_metric(df, "global_abs_br", "speed_km_s", "spearman"),
        "pearson_speed_equator_abs_br": corr_metric(df, "equator_abs_br", "speed_km_s", "pearson"),
        "pearson_speed_midlat_abs_br": corr_metric(df, "midlat_abs_br", "speed_km_s", "pearson"),
        "pearson_speed_global_abs_br": corr_metric(df, "global_abs_br", "speed_km_s", "pearson"),
        "polarity_n_equator": polarity_n(df, "equator_polarity"),
        "polarity_accuracy_equator": polarity_accuracy(df, "equator_polarity"),
        "polarity_n_global": polarity_n(df, "global_polarity"),
        "polarity_accuracy_global": polarity_accuracy(df, "global_polarity"),
    }

final_rows = []
by_cr_rows = []
sensitivity_rows = []

for rss in rss_values:
    path = Path(str(input_template).format(rss=rss))
    df = pd.read_csv(path)
    df["cr"] = pd.to_numeric(df["cr"], errors="coerce").astype("Int64")
    df["icme_flag"] = clean_bool(df["icme_flag"])

    final_all = df[df["cr"].isin(final_crs)].copy()
    final_no_icme = final_all[~final_all["icme_flag"]].copy()
    debug = df[df["cr"].isin(debug_crs)].copy()

    final_rows.append(metric_row(final_no_icme, rss, "final_no_icme"))
    sensitivity_rows.append(metric_row(final_all, rss, "final_all_rows"))
    sensitivity_rows.append(metric_row(final_no_icme, rss, "final_no_icme"))
    sensitivity_rows.append(metric_row(debug, rss, "debug_cr_2287"))

    for cr in final_crs:
        g = final_no_icme[final_no_icme["cr"] == cr].copy()
        row = metric_row(g, rss, "final_no_icme_by_cr")
        row["cr"] = cr
        by_cr_rows.append(row)

final_metrics = pd.DataFrame(final_rows)
by_cr_metrics = pd.DataFrame(by_cr_rows)
sensitivity = pd.DataFrame(sensitivity_rows)

final_metrics.to_csv(final_out, index=False)
by_cr_metrics.to_csv(by_cr_out, index=False)
sensitivity.to_csv(sensitivity_out, index=False)

print("Wrote:", final_out)
print("Wrote:", by_cr_out)
print("Wrote:", sensitivity_out)
print()
print("Final ICME-filtered metrics:")
print(final_metrics.to_string(index=False))
print()
print("By-CR ICME-filtered metrics:")
print(by_cr_metrics[["rss", "cr", "n_rows", "icme_rows", "median_speed_km_s", "spearman_speed_global_abs_br", "polarity_accuracy_global"]].to_string(index=False))
print()
print("Sensitivity table samples:")
print(sensitivity[["rss", "sample", "n_rows", "icme_rows", "median_speed_km_s", "spearman_speed_global_abs_br", "polarity_accuracy_global"]].to_string(index=False))
