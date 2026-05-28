from pathlib import Path
import numpy as np
import pandas as pd
import math

b = Path.cwd()

heights = ["rss2.0", "rss2.5", "rss3.0"]
out = b / "outputs" / "tables"
figout = b / "outputs" / "poster_figures"
out.mkdir(parents=True, exist_ok=True)
figout.mkdir(parents=True, exist_ok=True)

rows = []
detail_rows = []

def sign_nonzero(x):
    x = np.asarray(x, dtype=float)
    y = np.sign(x)
    y[~np.isfinite(x)] = np.nan
    y[y == 0] = np.nan
    return y

def normal_binom_p(k, n):
    if n <= 0:
        return np.nan
    z = abs(k - 0.5 * n) / math.sqrt(0.25 * n)
    return math.erfc(z / math.sqrt(2.0))

for h in heights:
    f = b / "data" / "processed" / "comparison" / h / "pfss_omni_ballistic_matched_rows.csv"

    if not f.exists():
        print(f"Missing {h}: {f}")
        continue

    d = pd.read_csv(f)

    needed = [
        "cr",
        "time",
        "ballistic_phase10_deg",
        "equator_signed_br",
        "global_signed_br",
        "bx_gse_nt",
        "by_gse_nt",
        "bz_gse_nt",
        "bmag_nt",
    ]

    miss = [c for c in needed if c not in d.columns]
    if miss:
        raise RuntimeError(f"{h} missing columns: {miss}")

    for c in needed:
        if c != "time":
            d[c] = pd.to_numeric(d[c], errors="coerce")

    d["rss"] = float(h.replace("rss", ""))

    # Near-Earth GSE -> RTN:
    # GSE +X points from Earth toward Sun.
    # RTN +R points from Sun outward to Earth.
    # Therefore R = -X_GSE, T = -Y_GSE, N = +Z_GSE.
    d["br_rtn_nt"] = -d["bx_gse_nt"]
    d["bt_rtn_nt"] = -d["by_gse_nt"]
    d["bn_rtn_nt"] = d["bz_gse_nt"]
    d["btan_rtn_nt"] = np.sqrt(d["bt_rtn_nt"]**2 + d["bn_rtn_nt"]**2)
    d["bmag_rtn_check_nt"] = np.sqrt(d["br_rtn_nt"]**2 + d["bt_rtn_nt"]**2 + d["bn_rtn_nt"]**2)
    d["rtn_spiral_angle_deg"] = np.degrees(np.arctan2(d["bt_rtn_nt"], d["br_rtn_nt"]))

    # Use 10 degree binned values to match the other ballistic correlation products.
    g = d.groupby(["rss", "cr", "ballistic_phase10_deg"], as_index=False).agg(
        n=("cr", "size"),
        bmag_mean=("bmag_nt", "mean"),
        br_rtn_mean=("br_rtn_nt", "mean"),
        bt_rtn_mean=("bt_rtn_nt", "mean"),
        bn_rtn_mean=("bn_rtn_nt", "mean"),
        btan_rtn_mean=("btan_rtn_nt", "mean"),
        rtn_spiral_angle_mean_deg=("rtn_spiral_angle_deg", "mean"),
        pfss_equator_signed_mean=("equator_signed_br", "mean"),
        pfss_global_signed_mean=("global_signed_br", "mean"),
    )

    g["omni_rtn_radial_polarity"] = sign_nonzero(g["br_rtn_mean"])
    g["pfss_equator_polarity"] = sign_nonzero(g["pfss_equator_signed_mean"])
    g["pfss_global_polarity"] = sign_nonzero(g["pfss_global_signed_mean"])

    for proxy_col, label in [
        ("pfss_equator_polarity", "equator_signed_br"),
        ("pfss_global_polarity", "global_signed_br"),
    ]:
        q = g[[proxy_col, "omni_rtn_radial_polarity"]].dropna().copy()

        direct = q[proxy_col].to_numpy() == q["omni_rtn_radial_polarity"].to_numpy()
        inverted = -q[proxy_col].to_numpy() == q["omni_rtn_radial_polarity"].to_numpy()

        n = int(len(q))
        direct_n = int(direct.sum())
        inverted_n = int(inverted.sum())

        rows.append({
            "rss": float(h.replace("rss", "")),
            "pfss_polarity_proxy": label,
            "omni_polarity": "RTN Br sign",
            "bins": n,
            "direct_agreement_fraction": direct_n / n if n else np.nan,
            "direct_agreement_count": direct_n,
            "direct_binomial_p_approx": normal_binom_p(direct_n, n),
            "inverted_agreement_fraction": inverted_n / n if n else np.nan,
            "inverted_agreement_count": inverted_n,
            "inverted_binomial_p_approx": normal_binom_p(inverted_n, n),
            "best_agreement_fraction": max(direct_n, inverted_n) / n if n else np.nan,
            "best_orientation": "direct" if direct_n >= inverted_n else "inverted",
        })

    detail_rows.append(g)

summary = pd.DataFrame(rows).sort_values(["rss", "best_agreement_fraction"], ascending=[True, False])
detail = pd.concat(detail_rows, ignore_index=True) if detail_rows else pd.DataFrame()

summary_path = out / "pfss_omni_rtn_imf_polarity_agreement_summary.csv"
detail_path = out / "pfss_omni_rtn_imf_polarity_agreement_binned.csv"

summary.to_csv(summary_path, index=False)
detail.to_csv(detail_path, index=False)

print("RTN IMF polarity agreement summary:")
print(summary.to_string(index=False))
print()
print(f"Saved summary: {summary_path}")
print(f"Saved binned detail: {detail_path}")
print()
print("Status: pass")
