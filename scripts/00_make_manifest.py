from pathlib import Path
import pandas as pd

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"

rows = [
    {
        "cr": 2281,
        "start_date": "2024-02-14",
        "end_date": "2024-03-12",
        "omni_start": "2024-02-12",
        "omni_end": "2024-03-19",
        "role": "ambient",
        "magnetogram_source": "GONG",
        "magnetogram_status": "",
        "omni_status": "",
        "pfss_status": "",
        "mapping_status": "",
        "notes": "Baseline early 2024 sample"
    },
    {
        "cr": 2283,
        "start_date": "2024-04-09",
        "end_date": "2024-05-06",
        "omni_start": "2024-04-07",
        "omni_end": "2024-05-13",
        "role": "ambient",
        "magnetogram_source": "GONG",
        "magnetogram_status": "",
        "omni_status": "",
        "pfss_status": "",
        "mapping_status": "",
        "notes": "Pre-storm or transition sample"
    },
    {
        "cr": 2284,
        "start_date": "2024-05-06",
        "end_date": "2024-06-02",
        "omni_start": "2024-05-04",
        "omni_end": "2024-06-09",
        "role": "stress_test",
        "magnetogram_source": "GONG",
        "magnetogram_status": "",
        "omni_status": "",
        "pfss_status": "",
        "mapping_status": "",
        "notes": "Disturbed May 2024 interval. Interpret separately from ambient sample."
    },
    {
        "cr": 2286,
        "start_date": "2024-06-29",
        "end_date": "2024-07-27",
        "omni_start": "2024-06-27",
        "omni_end": "2024-08-03",
        "role": "ambient",
        "magnetogram_source": "GONG",
        "magnetogram_status": "",
        "omni_status": "",
        "pfss_status": "",
        "mapping_status": "",
        "notes": "Post-active comparison sample"
    },
    {
        "cr": 2287,
        "start_date": "2024-07-27",
        "end_date": "2024-08-23",
        "omni_start": "2024-07-25",
        "omni_end": "2024-08-30",
        "role": "debug",
        "magnetogram_source": "GONG",
        "magnetogram_status": "",
        "omni_status": "",
        "pfss_status": "",
        "mapping_status": "",
        "notes": "Debug rotation only. Do not use for final claims unless promoted later."
    },
    {
        "cr": 2290,
        "start_date": "2024-10-16",
        "end_date": "2024-11-13",
        "omni_start": "2024-10-14",
        "omni_end": "2024-11-20",
        "role": "ambient",
        "magnetogram_source": "GONG",
        "magnetogram_status": "",
        "omni_status": "",
        "pfss_status": "",
        "mapping_status": "",
        "notes": "Near-maximum sample"
    },
]

df = pd.DataFrame(rows)

out = base / "metadata" / "rotation_manifest.csv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(f"Wrote manifest to: {out}")
print(df[["cr", "start_date", "end_date", "omni_start", "omni_end", "role"]])
