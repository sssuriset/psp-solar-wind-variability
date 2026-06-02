from pathlib import Path
import pandas as pd

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"
manifest_path = base / "metadata" / "rotation_manifest.csv"

df = pd.read_csv(manifest_path)

# Empty status columns were read as float NaN.
# Force them to string/object columns before writing text statuses.
status_cols = [
    "magnetogram_status",
    "omni_status",
    "pfss_status",
    "mapping_status",
    "notes",
]

for col in status_cols:
    df[col] = df[col].fillna("").astype(str)

mask = df["cr"] == 2287

df.loc[mask, "pfss_status"] = "pass_sample_pfss"

old_note = df.loc[mask, "notes"].iloc[0]
extra_note = "PFSS sample sanity test passed at rss 2.0, 2.5, 3.0 Rs."

if extra_note not in old_note:
    df.loc[mask, "notes"] = old_note + " | " + extra_note

df.to_csv(manifest_path, index=False)

print("Updated CR 2287 PFSS debug status.")
print(df[df["cr"] == 2287][["cr", "role", "pfss_status", "notes"]].to_string(index=False))
