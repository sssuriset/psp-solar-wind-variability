from pathlib import Path
import pandas as pd

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"
manifest_path = base / "metadata" / "rotation_manifest.csv"

df = pd.read_csv(manifest_path)

cols = [
    "cr",
    "role",
    "start_date",
    "end_date",
    "magnetogram_status",
    "omni_status",
    "pfss_status",
    "mapping_status",
    "notes",
]

print()
print("Phase 2 feasibility manifest")
print()

print(df[cols].to_string(index=False))

ambient = df[df["role"] == "ambient"]
ambient_mag_pass = ambient["magnetogram_status"].fillna("").str.contains("pass|backup", regex=True).sum()

print()
print(f"Ambient rotations with local magnetogram candidates: {ambient_mag_pass} / {len(ambient)}")

debug = df[df["role"] == "debug"]
if not debug.empty:
    print()
    print("Debug rotation:")
    print(debug[["cr", "magnetogram_status", "pfss_status", "mapping_status"]].to_string(index=False))
