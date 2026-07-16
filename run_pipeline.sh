#!/usr/bin/env bash
# Full pipeline for the PFSS source-surface height sensitivity study.
# Run from anywhere; scripts anchor paths to the repo root.
#
# Requirements before running:
#   pip install -r requirements.txt
#   export JSOC_EMAIL='your JSOC-registered email'   (http://jsoc.stanford.edu/ajax/register_email.html)
#
# Stages 1-3 download data (JSOC HMI synoptic maps, CDAWeb OMNI hourly data)
# and are slow on first run. Later stages are pure local computation.

set -euo pipefail
cd "$(dirname "$0")/scripts"

HEIGHTS="2.0 2.5 3.0"

# 1. Rotation manifest (defines the CR sample)
python 01_build_manifest.py

# 2. HMI synoptic magnetograms from JSOC (requires JSOC_EMAIL)
python 02_download_hmi.py

# 3. OMNI hourly data from CDAWeb -> data/raw/omni/omni_cr*.csv
python 03_download_omni.py

# 4. Clean and phase-tag OMNI -> data/processed/omni/
python 04_clean_omni.py

# 5. PFSS extrapolations, all heights -> data/processed/pfss/hmi/
python 05_run_pfss.py --rss $HEIGHTS

# 6. Longitude proxies, all heights -> data/processed/pfss/hmi_longitude_proxies/
python 06_build_proxies.py --rss $HEIGHTS

# 7. Per-height ballistic matching, phase correlations, lag scans, null tests
for RSS in $HEIGHTS; do
    python 07_match_ballistic.py    --rss "$RSS"
    python 08_phase_correlations.py --rss "$RSS"
    python 09_scan_lags.py          --rss "$RSS"
    python 10_null_tests.py         --rss "$RSS"
done

# 8. Cross-height secondary analyses
python 12_earth_track.py
python 13_hss_timing.py
python 14_polarity_agreement.py
python 15_expansion_proxy.py

# 9. ICME robustness check (requires the ICME interval catalog; see README)
ICME_CATALOG="../data/external/wind_icme_intervals_2024.csv"
if [ -f "$ICME_CATALOG" ]; then
    python 11_flag_icme.py
    RUN_ICME=1
else
    echo "ICME catalog not found at data/external/; skipping ICME-filtered results."
    RUN_ICME=0
fi

# 10. Locked final sample and headline tables -> outputs/tables/final/, outputs/results/final/
python 16_final_results.py

# 11. ICME-filtered final results and row-level sensitivity
if [ "$RUN_ICME" = "1" ]; then
    python 16_final_results.py --icme-filter
    python 17_icme_row_sensitivity.py
fi

echo "Pipeline complete. See outputs/results/final/summary.txt"
