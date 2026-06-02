from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunkit_magex import pfss

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"
outdir = base / "outputs" / "phase2"
outdir.mkdir(parents=True, exist_ok=True)

print("PFSS sample sanity test")
print("This uses sample GONG data, not the selected SPD 57 science rotations.")
print()

print("Loading sample GONG magnetogram...")
gong_file = pfss.sample_data.get_gong_map()
print(f"Sample file: {gong_file}")

gong_map = sunpy.map.Map(gong_file)

print("Map summary:")
print(f"  Shape: {gong_map.data.shape}")
print(f"  Unit: {gong_map.unit}")
print(f"  Date: {gong_map.date}")

finite = np.isfinite(gong_map.data)
finite_fraction = finite.sum() / gong_map.data.size
print(f"  Finite fraction: {finite_fraction:.3f}")

if finite_fraction < 0.80:
    raise RuntimeError("Map has too many missing values.")

nrho = 40
rss_values = [2.0, 2.5, 3.0]

results = []

for rss in rss_values:
    print()
    print(f"Running PFSS with source surface height rss = {rss} Rs...")

    pfss_input = pfss.Input(gong_map, nrho, rss)
    pfss_output = pfss.pfss(pfss_input)

    nsteps = 30

    lon_1d = np.linspace(0, 2 * np.pi, nsteps * 2 + 1)
    lat_1d = np.arcsin(np.linspace(-1, 1, nsteps + 1))

    lon, lat = np.meshgrid(lon_1d, lat_1d, indexing="ij")

    seeds = SkyCoord(
        lon.ravel() * u.rad,
        lat.ravel() * u.rad,
        const.R_sun,
        frame=pfss_output.coordinate_frame,
    )

    tracer = pfss.tracing.PerformanceTracer()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        field_lines = tracer.trace(seeds, pfss_output)

    polarities = field_lines.polarities

    open_count = int(np.count_nonzero(polarities != 0))
    closed_count = int(np.count_nonzero(polarities == 0))
    total_count = len(polarities)

    print(f"  Total traced lines: {total_count}")
    print(f"  Open lines: {open_count}")
    print(f"  Closed lines: {closed_count}")

    if open_count == 0 or closed_count == 0:
        status = "warning_all_open_or_all_closed"
    else:
        status = "pass"

    results.append({
        "rss": rss,
        "total": total_count,
        "open": open_count,
        "closed": closed_count,
        "status": status,
    })

figure_path = outdir / "pfss_sample_result.txt"

with open(figure_path, "w") as f:
    f.write("PFSS sample sanity test results\n")
    for row in results:
        f.write(
            f"rss={row['rss']}, "
            f"total={row['total']}, "
            f"open={row['open']}, "
            f"closed={row['closed']}, "
            f"status={row['status']}\n"
        )

print()
print("PFSS sample results:")
for row in results:
    print(
        f"  rss={row['rss']}: "
        f"open={row['open']}, closed={row['closed']}, status={row['status']}"
    )

print()
print(f"Saved result file: {figure_path}")

if not any(row["status"] == "pass" for row in results):
    raise RuntimeError("PFSS ran, but all source-surface heights gave suspicious open/closed results.")

print()
print("Status: pass")
