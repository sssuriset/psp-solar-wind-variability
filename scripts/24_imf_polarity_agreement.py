from pathlib import Path
import numpy as np
import pandas as pd
import math

b = Path.cwd()

heights = ["rss2.0", "rss2.5", "rss3.0"]
out = b / "outputs" / "tables"
out.mkdir(parents=True, exist_ok=True)

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
        "ballistic_phase10_deg",
        "equator_signed_br",
        "global_signed_br",
        "bx_gse_nt",
        "bmag_nt",
    ]

    miss = [c for c in needed if c not in d.columns]
    if miss:
        raise RuntimeError(f"{h} missing columns: {miss}")

    d["rss"] = float(h.replace("rss", ""))

    # GSE X points approximately Sunward, so -Bx_GSE is used as a near-Earth radial IMF polarity proxy.
    # This is a proxy, not a full RTN polarity transformation.
    d["omni_radial_polarity_proxy"] = sign_nonzero(-d["bx_gse_nt"])
    d["pfss_equator_polarity_proxy"] = sign_nonzero(d["equator_signed_br"])
    d["pfss_global_polarity_proxy"] = sign_nonzero(d["global_signed_br"])

    # Use 10 degree binned values so the sample matches the correlation products better than hourly rows.
    g = d.groupby(["rss", "cr", "ballistic_phase10_deg"], as_index=False).agg(
        n=("cr", "size"),
        bx_gse_mean=("bx_gse_nt", "mean"),
        bmag_mean=("bmag_nt", "mean"),
        pfss_equator_signed_mean=("equator_signed_br", "mean"),
        pfss_global_signed_mean=("global_signed_br", "mean"),
    )

    g["omni_radial_polarity_proxy"] = sign_nonzero(-g["bx_gse_mean"])
    g["pfss_equator_polarity_proxy"] = sign_nonzero(g["pfss_equator_signed_mean"])
    g["pfss_global_polarity_proxy"] = sign_nonzero(g["pfss_global_signed_mean"])

    for proxy_col, label in [
        ("pfss_equator_polarity_proxy", "equator_signed_br"),
        ("pfss_global_polarity_proxy", "global_signed_br"),
    ]:
        q = g[[proxy_col, "omni_radial_polarity_proxy"]].dropna().copy()

        direct = q[proxy_col].to_numpy() == q["omni_radial_polarity_proxy"].to_numpy()
        inverted = -q[proxy_col].to_numpy() == q["omni_radial_polarity_proxy"].to_numpy()

        n = int(len(q))
        direct_n = int(direct.sum())
        inverted_n = int(inverted.sum())

        rows.append({
            "rss": float(h.replace("rss", "")),
            "pfss_polarity_proxy": label,
            "omni_polarity_proxy": "-Bx_GSE sign",
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

summary_path = out / "pfss_omni_imf_polarity_agreement_summary.csv"
detail_path = out / "pfss_omni_imf_polarity_agreement_binned.csv"

summary.to_csv(summary_path, index=False)
detail.to_csv(detail_path, index=False)

print("IMF polarity agreement summary:")
print(summary.to_string(index=False))
print()
print(f"Saved summary: {summary_path}")
print(f"Saved binned detail: {detail_path}")
print()
print("Status: pass")
