from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ICME_PATH = ROOT / "data" / "external" / "wind_icme_intervals_2024.csv"

FINAL_CRS = {2281, 2283, 2284, 2286, 2290}
DEBUG_CRS = {2287}
TAGS = ["rss2.0", "rss2.5", "rss3.0"]

if not ICME_PATH.exists():
    raise FileNotFoundError(
        f"ICME interval catalog not found: {ICME_PATH}\n"
        "Expected columns: icme_start_utc, icme_end_utc, save_name, mo_type. "
        "See the README data provenance section."
    )

icme = pd.read_csv(ICME_PATH)
icme["icme_start_utc"] = pd.to_datetime(icme["icme_start_utc"], utc=True)
icme["icme_end_utc"] = pd.to_datetime(icme["icme_end_utc"], utc=True)

summary_rows = []

for tag in TAGS:
    rss_dir = ROOT / "data" / "processed" / "comparison" / tag
    in_path = rss_dir / "pfss_omni_ballistic_matched_rows.csv"
    if not in_path.exists():
        raise FileNotFoundError(f"Missing matched rows (run 07_match_ballistic.py --rss first): {in_path}")

    df = pd.read_csv(in_path)
    df["time"] = pd.to_datetime(df["time"], utc=True)

    df["icme_flag"] = False
    df["icme_save_name"] = ""
    df["icme_type"] = ""

    for _, r in icme.iterrows():
        mask = (df["time"] >= r["icme_start_utc"]) & (df["time"] <= r["icme_end_utc"])
        df.loc[mask, "icme_flag"] = True
        df.loc[mask, "icme_save_name"] = str(r["save_name"])
        df.loc[mask, "icme_type"] = str(r["mo_type"])

    out_path = rss_dir / "pfss_omni_ballistic_matched_rows_icme_flagged.csv"
    df.to_csv(out_path, index=False)

    for cr, g in df.groupby("cr"):
        cr_int = int(cr)
        role = "final" if cr_int in FINAL_CRS else "debug" if cr_int in DEBUG_CRS else "other"
        no_icme = g[~g["icme_flag"]]

        summary_rows.append({
            "rss": tag,
            "cr": cr_int,
            "role": role,
            "n_rows": len(g),
            "icme_rows": int(g["icme_flag"].sum()),
            "icme_fraction": round(float(g["icme_flag"].mean()), 4),
            "median_speed_all_km_s": g["speed_km_s"].median(),
            "median_speed_no_icme_km_s": no_icme["speed_km_s"].median() if len(no_icme) else "",
            "median_bmag_all_nt": g["bmag_nt"].median(),
            "median_bmag_no_icme_nt": no_icme["bmag_nt"].median() if len(no_icme) else "",
        })

    print("Wrote:", out_path)

summary = pd.DataFrame(summary_rows).sort_values(["rss", "cr"])
summary_path = ROOT / "data" / "processed" / "comparison" / "icme_overlap_summary.csv"
summary.to_csv(summary_path, index=False)

print()
print("Wrote:", summary_path)
print()
print(summary.to_string(index=False))
