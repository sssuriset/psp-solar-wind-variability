from pathlib import Path
import pandas as pd
import numpy as np

b = Path.cwd()
od = b / "data" / "raw" / "omni"
pdn = b / "data" / "processed" / "omni"
man = b / "metadata" / "rotation_manifest.csv"

print(f"Project base: {b}")
print(f"Raw OMNI dir: {od}")
print()

if not od.exists():
    raise FileNotFoundError(f"Missing OMNI dir: {od}")

fs = sorted(od.glob("omni_cr*.csv"))

if len(fs) == 0:
    raise RuntimeError("No omni_cr*.csv files found.")

print("OMNI files:")
for f in fs:
    print(f"- {f.name}")
print()

if man.exists():
    m = pd.read_csv(man)
    print("Manifest OMNI status:")
    cs = [c for c in ["cr", "omni_status", "pfss_status", "proxy_status"] if c in m.columns]
    print(m[cs].to_string(index=False))
    print()

for f in fs:
    print("=" * 80)
    print(f.name)

    d = pd.read_csv(f)

    print(f"Rows: {len(d)}")
    print(f"Columns: {len(d.columns)}")
    print()

    print("Column names:")
    for c in d.columns:
        print(f"- {c}")

    print()
    print("First 3 rows:")
    print(d.head(3).to_string(index=False))

    print()
    print("Numeric column quick check:")
    ns = d.select_dtypes(include=[np.number]).columns.tolist()

    if len(ns) == 0:
        print("No numeric columns detected.")
    else:
        q = []
        for c in ns:
            x = pd.to_numeric(d[c], errors="coerce")
            q.append({
                "col": c,
                "n": int(x.notna().sum()),
                "nan": int(x.isna().sum()),
                "min": float(x.min()) if x.notna().any() else np.nan,
                "max": float(x.max()) if x.notna().any() else np.nan,
                "mean": float(x.mean()) if x.notna().any() else np.nan,
            })

        s = pd.DataFrame(q)
        print(s.to_string(index=False))

    print()
    print("Status: inspected")
    print()

print("=" * 80)
print("Status: pass")
