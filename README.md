# PFSS Source-Surface Height Sensitivity

How sensitive are coronal-to-1 AU magnetic mapping predictions to the PFSS source-surface height? This repository builds potential-field source-surface (PFSS) models from HMI synoptic magnetograms at three source-surface heights (2.0, 2.5, 3.0 solar radii), ballistically maps the sub-Earth open-field structure to 1 AU, and tests the mapped longitude proxies against OMNI hourly solar-wind observations across five Carrington rotations.

Presented at the AAS Solar Physics Division meeting SPD 57 (Baltimore, August 2026).

## Headline result

The equatorial open-flux proxy correlates with near-Earth field strength (Spearman r about 0.30 for the 180-point pooled sample), and the correlation survives a phase-shift null test (two-sided p between 0.002 and 0.007 depending on height and sample). Height sensitivity is modest but consistent, with rss = 2.5 giving the strongest correlation in the pooled sample and in 4 of 5 leave-one-rotation-out subsets. Full tables are in `outputs/tables/final/` and the plain-text summary in `outputs/results/final/summary.txt`.

## Sample and caveats

The final science sample is CRs 2281, 2283, 2284, 2286, 2290 (CR 2287 was used for debugging only and is excluded). CR 2284 is a deliberate stress-test rotation and results are reported both with and without it. Magnetograms are HMI-only (`hmi.synoptic_mr_polfil_720s`). IMF sector polarity uses the approximation Br_RTN = -Bx_GSE. With n = 5 rotations, pooled correlations are supported by the shift-null p-values rather than raw r alone.

## Data provenance

- HMI synoptic maps from JSOC, series `hmi.synoptic_mr_polfil_720s`, fetched via sunpy Fido. Requires a JSOC-registered email in the `JSOC_EMAIL` environment variable (register at http://jsoc.stanford.edu/ajax/register_email.html).
- OMNI hourly merged data (`OMNI2_H0_MRG1HR`) fetched from CDAWeb via cdasws.
- ICME intervals in `data/external/wind_icme_intervals_2024.csv` with columns `icme_start_utc`, `icme_end_utc`, `save_name`, `mo_type`, derived from [insert ICME catalog citation and version here]. This small reference file is committed so the ICME robustness check reproduces from a clean checkout.

No raw or intermediate data is committed. Only the final tables (`outputs/tables/final/`) and result summaries (`outputs/results/final/`) are tracked.

## Setup

```
pip install -r requirements.txt
export JSOC_EMAIL='your registered email'
```

PFSS extrapolation uses sunkit-magex (the maintained successor to pfsspy).

## Running the pipeline

```
./run_pipeline.sh
```

runs every stage in order. The stages, individually:

| Stage | Script | Output |
|---|---|---|
| 1 | `01_build_manifest.py` | `metadata/rotation_manifest.csv`, the CR sample definition |
| 2 | `02_download_hmi.py` | HMI synoptic FITS in `data/raw/` |
| 3 | `03_download_omni.py` | OMNI hourly CSVs in `data/raw/omni/` |
| 4 | `04_clean_omni.py` | cleaned, phase-tagged OMNI table in `data/processed/omni/` |
| 5 | `05_run_pfss.py --rss 2.0 2.5 3.0` | PFSS solutions and metrics in `data/processed/pfss/hmi/` |
| 6 | `06_build_proxies.py --rss 2.0 2.5 3.0` | longitude proxy profiles in `data/processed/pfss/hmi_longitude_proxies/` |
| 7 | `07_match_ballistic.py`, `08_phase_correlations.py`, `09_scan_lags.py`, `10_null_tests.py`, each with `--rss X` | ballistic matching, phase correlations, fixed-lag scans, and shift-null tests per height in `data/processed/comparison/rssX/` |
| 8 | `12_earth_track.py`, `13_hss_timing.py`, `14_polarity_agreement.py`, `15_expansion_proxy.py` | cross-height secondary analyses (Earth-connected track coverage, high-speed-stream timing, IMF polarity agreement, radial expansion proxy) |
| 9 | `11_flag_icme.py` | flags matched hours falling inside catalogued ICME intervals; writes `_icme_flagged` matched rows per height |
| 10 | `16_final_results.py` | locked science sample, pooled and per-CR correlations, leave-one-out table, null tests in `outputs/tables/final/` and `outputs/results/final/` |
| 11 | `16_final_results.py --icme-filter` | same analysis with in-ICME hours removed, in `outputs/tables/final_no_icme/` and `outputs/results/final_no_icme/` |
| 12 | `17_icme_row_sensitivity.py` | row-level correlation and polarity sensitivity table comparing all hours, ICME-removed hours, and the debug rotation |

`check_local_magnetograms.py`, `check_fits_readability.py`, `inspect_omni_dataset.py`, and `inspect_omni_table.py` are optional inspection utilities and are not part of the pipeline.

## License

MIT. See `LICENSE`.
