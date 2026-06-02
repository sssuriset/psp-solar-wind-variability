from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

b = Path.cwd()

rd = b / "data" / "raw" / "omni"
od = b / "data" / "processed" / "omni"
fd = b / "outputs" / "figures" / "omni"
td = b / "outputs" / "tables"
mf = b / "metadata" / "rotation_manifest.csv"

od.mkdir(parents=True, exist_ok=True)
fd.mkdir(parents=True, exist_ok=True)
td.mkdir(parents=True, exist_ok=True)

m = pd.read_csv(mf)

if "omni_feature_status" not in m.columns:
    m["omni_feature_status"] = ""

for c in ["omni_feature_status", "notes"]:
    if c not in m.columns:
        m[c] = ""
    m[c] = m[c].fillna("").astype(str)

def gc(x, ks):
    for k in ks:
        if k in x.columns:
            return k
    return None

def st(x, nm):
    z = pd.to_numeric(x, errors="coerce")
    n = int(z.notna().sum())
    q = {
        f"{nm}_n": n,
        f"{nm}_nan": int(z.isna().sum()),
        f"{nm}_mean": float(z.mean()) if n else np.nan,
        f"{nm}_median": float(z.median()) if n else np.nan,
        f"{nm}_std": float(z.std()) if n else np.nan,
        f"{nm}_min": float(z.min()) if n else np.nan,
        f"{nm}_max": float(z.max()) if n else np.nan,
        f"{nm}_p10": float(z.quantile(0.10)) if n else np.nan,
        f"{nm}_p90": float(z.quantile(0.90)) if n else np.nan,
    }
    return q

ss = [
    "start_time", "start", "start_utc", "start_date", "start_datetime",
    "cr_start", "cr_start_time", "cr_start_utc",
    "rotation_start", "rotation_start_time", "t0", "begin", "begin_time"
]

ee = [
    "end_time", "end", "end_utc", "end_date", "end_datetime",
    "cr_end", "cr_end_time", "cr_end_utc",
    "rotation_end", "rotation_end_time", "t1", "finish", "finish_time"
]

sc = gc(m, ss)
ec = gc(m, ee)

fs = sorted(rd.glob("omni_cr*.csv"))

if len(fs) == 0:
    raise RuntimeError("No OMNI CSV files found.")

print(f"Project base: {b}")
print(f"OMNI files: {len(fs)}")
print(f"Manifest window columns: start={sc}, end={ec}")
print()

rows = []
tabs = []

for f in fs:
    mt = re.search(r"cr(\d+)", f.name)

    if mt is None:
        print(f"Skipping file with no CR number: {f.name}")
        continue

    cr = int(mt.group(1))

    print("=" * 80)
    print(f"CR {cr}")
    print(f"File: {f.name}")

    d = pd.read_csv(f)

    need = ["time", "speed_km_s", "bx_gse_nt", "by_gse_nt", "bz_gse_nt", "bmag_nt", "density_cm3"]
    miss = [c for c in need if c not in d.columns]

    if len(miss) > 0:
        print(f"Missing columns: {miss}")
        m.loc[m["cr"].astype(int) == cr, "omni_feature_status"] = "fail_missing_columns"
        continue

    d["time"] = pd.to_datetime(d["time"], utc=True, errors="coerce")
    d = d.dropna(subset=["time"]).sort_values("time").reset_index(drop=True)

    for c in need[1:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    wk = "full_file"
    a = d["time"].min()
    z = d["time"].max()

    q = m[m["cr"].astype(int) == cr]

    if len(q) == 1 and sc is not None and ec is not None:
        aa = pd.to_datetime(q.iloc[0][sc], utc=True, errors="coerce")
        zz = pd.to_datetime(q.iloc[0][ec], utc=True, errors="coerce")

        if pd.notna(aa) and pd.notna(zz):
            wk = "manifest"
            a = aa
            z = zz
            d = d[(d["time"] >= a) & (d["time"] < z)].copy()

    if len(d) == 0:
        print("No rows after time window crop.")
        m.loc[m["cr"].astype(int) == cr, "omni_feature_status"] = "fail_empty_after_crop"
        continue

    d = d.reset_index(drop=True)

    sp = (z - a).total_seconds()

    if sp <= 0:
        a = d["time"].min()
        z = d["time"].max()
        sp = max((z - a).total_seconds(), 1)

    d["cr"] = cr
    d["tday"] = (d["time"] - a).dt.total_seconds() / 86400.0
    d["rot_frac"] = (d["time"] - a).dt.total_seconds() / sp
    d["lon_deg_proxy"] = 360.0 * d["rot_frac"]

    goodn = d["density_cm3"] > 0
    d["pdyn_npa"] = 1.6726e-6 * d["density_cm3"] * d["speed_km_s"] ** 2
    d["va_km_s"] = np.nan
    d.loc[goodn, "va_km_s"] = 21.812 * d.loc[goodn, "bmag_nt"] / np.sqrt(d.loc[goodn, "density_cm3"])
    d["ma"] = d["speed_km_s"] / d["va_km_s"]

    op = od / f"omni_cr{cr}_clean_with_phase.csv"
    d.to_csv(op, index=False)

    row = {
        "cr": cr,
        "rows": int(len(d)),
        "window_source": wk,
        "time_min": d["time"].min().isoformat(),
        "time_max": d["time"].max().isoformat(),
        "duration_days": float((d["time"].max() - d["time"].min()).total_seconds() / 86400.0),
        "clean_path": str(op),
    }

    for c, nm in [
        ("speed_km_s", "speed"),
        ("bmag_nt", "bmag"),
        ("density_cm3", "density"),
        ("bx_gse_nt", "bx"),
        ("by_gse_nt", "by"),
        ("bz_gse_nt", "bz"),
        ("pdyn_npa", "pdyn"),
        ("va_km_s", "va"),
        ("ma", "ma"),
    ]:
        row.update(st(d[c], nm))

    rows.append(row)
    tabs.append(d)

    msk = m["cr"].astype(int) == cr
    m.loc[msk, "omni_feature_status"] = "pass_omni_features"

    nt = "OMNI feature and phase table created."
    old = m.loc[msk, "notes"].iloc[0]

    if nt not in old:
        m.loc[msk, "notes"] = old + " | " + nt

    print(f"Rows used: {len(d)}")
    print(f"Window source: {wk}")
    print(f"Saved clean table: {op}")
    print(f"Mean speed: {d['speed_km_s'].mean():.3f} km/s")
    print(f"Mean |B|: {d['bmag_nt'].mean():.3f} nT")
    print(f"Mean density: {d['density_cm3'].mean():.3f} cm^-3")

if len(rows) == 0:
    raise RuntimeError("No OMNI feature tables were created.")

s = pd.DataFrame(rows).sort_values("cr")
ap = pd.concat(tabs, ignore_index=True).sort_values(["cr", "time"])

sp = od / "omni_selected_crs_feature_summary.csv"
apath = od / "omni_selected_crs_clean_with_phase.csv"
cp = td / "omni_feature_summary_compact.csv"

s.to_csv(sp, index=False)
s.to_csv(cp, index=False)
ap.to_csv(apath, index=False)
m.to_csv(mf, index=False)

for y, lab, fn in [
    ("speed_km_s", "OMNI solar wind speed, km/s", "omni_speed_by_rotation_phase.png"),
    ("bmag_nt", "OMNI magnetic field magnitude, nT", "omni_bmag_by_rotation_phase.png"),
    ("density_cm3", "OMNI proton density, cm^-3", "omni_density_by_rotation_phase.png"),
    ("pdyn_npa", "OMNI dynamic pressure, nPa", "omni_pdyn_by_rotation_phase.png"),
]:
    plt.figure(figsize=(10, 5))

    for cr in sorted(ap["cr"].unique()):
        q = ap[ap["cr"] == cr]
        plt.plot(q["lon_deg_proxy"], q[y], label=f"CR {cr}", linewidth=1)

    plt.xlabel("Rotation phase proxy, degrees")
    plt.ylabel(lab)
    plt.title(lab)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(fd / fn, dpi=200)
    plt.close()

print()
print("OMNI feature summary:")
print(s[[
    "cr",
    "rows",
    "window_source",
    "duration_days",
    "speed_mean",
    "speed_p10",
    "speed_p90",
    "bmag_mean",
    "density_mean",
    "pdyn_mean",
    "ma_mean",
]].to_string(index=False))

print()
print(f"Saved summary: {sp}")
print(f"Saved compact table: {cp}")
print(f"Saved all clean phase data: {apath}")
print(f"Saved figures: {fd}")

print()
print("Manifest OMNI feature status:")
print(m[["cr", "omni_status", "pfss_status", "proxy_status", "omni_feature_status"]].to_string(index=False))

print()
print("Status: pass")
