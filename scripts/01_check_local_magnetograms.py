from pathlib import Path
import pandas as pd

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"
manifest_path = base / "metadata" / "rotation_manifest.csv"

df = pd.read_csv(manifest_path)

gong_dir = base / "data" / "raw" / "magnetograms" / "gong"
hmi_dir = base / "data" / "raw" / "magnetograms" / "hmi"

def find_files_for_cr(folder, cr):
    patterns = [
        f"*{cr}*.fits",
        f"*{cr}*.fts",
        f"*cr{cr}*.fits",
        f"*CR{cr}*.fits",
    ]

    matches = []
    for pattern in patterns:
        matches.extend(folder.glob(pattern))

    return sorted(set(matches))

for i, row in df.iterrows():
    cr = int(row["cr"])

    gong_files = find_files_for_cr(gong_dir, cr)
    hmi_files = find_files_for_cr(hmi_dir, cr)

    if gong_files:
        df.loc[i, "magnetogram_status"] = "pass_gong_local"
        df.loc[i, "notes"] = str(row["notes"]) + f" | Found GONG file: {gong_files[0].name}"
    elif hmi_files:
        df.loc[i, "magnetogram_status"] = "backup_hmi_local"
        df.loc[i, "notes"] = str(row["notes"]) + f" | Found HMI file: {hmi_files[0].name}"
    else:
        df.loc[i, "magnetogram_status"] = "not_found_local"

df.to_csv(manifest_path, index=False)

print("Updated magnetogram local-file status.")
print(df[["cr", "role", "magnetogram_source", "magnetogram_status"]])
