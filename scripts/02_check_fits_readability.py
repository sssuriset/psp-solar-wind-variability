from pathlib import Path
import numpy as np
from astropy.io import fits

base = Path.home() / "Desktop" / "astro_project" / "pfss_omni_spd57"
mag_dirs = [
    base / "data" / "raw" / "magnetograms" / "gong",
    base / "data" / "raw" / "magnetograms" / "hmi",
]

fits_files = []
for folder in mag_dirs:
    fits_files.extend(folder.glob("*.fits"))
    fits_files.extend(folder.glob("*.fts"))

if not fits_files:
    print("No FITS files found yet.")
    print("This is okay before the magnetogram download step.")
    raise SystemExit

for path in sorted(fits_files):
    print()
    print(f"Checking: {path}")

    try:
        with fits.open(path) as hdul:
            print(f"  HDUs: {len(hdul)}")

            found_data = False

            for hdu_index, hdu in enumerate(hdul):
                data = hdu.data

                if data is None:
                    continue

                arr = np.asarray(data)

                if arr.ndim < 2:
                    continue

                finite = np.isfinite(arr)
                finite_fraction = finite.sum() / arr.size

                print(f"  First 2D data HDU: {hdu_index}")
                print(f"  Shape: {arr.shape}")
                print(f"  Finite fraction: {finite_fraction:.3f}")

                if finite_fraction > 0.80:
                    print("  Status: pass")
                else:
                    print("  Status: warning_low_finite_fraction")

                found_data = True
                break

            if not found_data:
                print("  Status: fail_no_2d_data")

    except Exception as exc:
        print(f"  Status: fail_read_error")
        print(f"  Error: {exc}")
