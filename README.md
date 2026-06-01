# PFSS Source-Surface Sensitivity and Near-Earth Solar-Wind Structure

This project tests how PFSS source-surface height and ballistic longitude mapping assumptions affect agreement between photospheric magnetic-field structure and near-Earth solar-wind measurements.

The project is not a solar-wind prediction model. It is a model-assumption sensitivity study. The main question is whether changing the PFSS source-surface height changes the agreement between source-surface magnetic-field proxies and OMNI near-Earth solar-wind quantities.

## Scientific question

How sensitive are PFSS and ballistic-mapping comparisons to the assumed source-surface height?

Tested source-surface heights:

| Height label | Meaning |
|---|---|
| `rss = 2.0` | source surface at 2.0 solar radii |
| `rss = 2.5` | source surface at 2.5 solar radii |
| `rss = 3.0` | source surface at 3.0 solar radii |

## Data sources

| Data product | Role |
|---|---|
| HMI synoptic magnetograms | photospheric magnetic-field input for PFSS |
| PFSS source-surface maps | magnetic-field structure at the model source surface |
| OMNI near-Earth solar-wind data | comparison data for solar-wind speed, magnetic-field magnitude, density, dynamic pressure, and Alfvénic quantities |

## Final science sample

Final science sample:

`CR 2281, 2283, 2284, 2286, 2290`

CR 2287 was used as a debug rotation and is excluded from final science metrics.

| Carrington rotation | Final role |
|---|---|
| 2281 | final science sample |
| 2283 | final science sample |
| 2284 | final science sample and disturbed/stress-test rotation |
| 2286 | final science sample |
| 2287 | debug-only, excluded from final science metrics |
| 2290 | final science sample |

## Pipeline summary

The pipeline does the following:

1. Load selected HMI synoptic magnetograms
2. Run PFSS at multiple source-surface heights
3. Extract longitude-dependent source-surface magnetic-field proxies
4. Apply ballistic longitude mapping to OMNI near-Earth measurements
5. Bin mapped quantities into 10-degree longitude bins
6. Compare PFSS magnetic-field proxies against OMNI solar-wind quantities
7. Run height sensitivity, ambient/disturbed comparison, and leave-one-rotation-out validation

Final target relation:

`pfss_equator_abs_mean` vs `bmag_mean`

The comparison uses Spearman correlation across 10-degree ballistic longitude bins.

## Final result

In the locked five-rotation sample, `rss = 2.5` gives the strongest full-sample association between equatorial unsigned PFSS magnetic-field strength and OMNI magnetic-field magnitude.

| Sample | rss | Spearman r | bins n | shift-null p |
|---|---:|---:|---:|---:|
| final sample with CR 2284 | 2.0 | 0.280349 | 180 | 0.005497 |
| final sample with CR 2284 | 2.5 | 0.298558 | 180 | 0.001999 |
| final sample with CR 2284 | 3.0 | 0.252821 | 180 | 0.012994 |
| ambient comparison without CR 2284 | 2.0 | 0.310749 | 144 | 0.004498 |
| ambient comparison without CR 2284 | 2.5 | 0.313508 | 144 | 0.006997 |
| ambient comparison without CR 2284 | 3.0 | 0.251493 | 144 | 0.032984 |

Removing CR 2284 does not remove the association. The association is slightly stronger in the ambient-only comparison for `rss = 2.5`.

Leave-one-rotation-out validation shows that `rss = 2.5` is strongest in four of five dropout tests. `rss = 2.0` is strongest only when CR 2290 is removed.

## Output files supporting the final result

| File | Purpose |
|---|---|
| `outputs/tables/final_sample/final_abstract_support_summary.csv` | single table for abstract and README support |
| `outputs/tables/final_sample/final_sample_height_sensitivity_eqabs_bmag_spearman.csv` | final height-sensitivity result for the target relation |
| `outputs/tables/final_sample/leave_one_rotation_out_eqabs_bmag_spearman.csv` | leave-one-rotation-out validation |
| `outputs/tables/final_sample/final_sample_all_correlations_with_null.csv` | broader correlation table with circular-shift null tests |
| `outputs/tables/final_sample/final_sample_ballistic_lag_summary.csv` | ballistic lag and longitude-shift summary |
| `outputs/tables/final_sample/final_sample_rtn_imf_polarity_agreement.csv` | IMF polarity agreement check |
| `outputs/results/final_sample/final_abstract_support_summary.txt` | text summary of final abstract-support values |
| `outputs/results/final_sample/final_sample_summary.txt` | final metric summary |

## Reproducing the final metric tables

From the repository root:

    conda activate pfss-spd57
    python scripts/28_final_sample_metrics.py

This assumes that the height-specific matched comparison tables already exist at:

    data/processed/comparison/rss2.0/pfss_omni_ballistic_matched_rows.csv
    data/processed/comparison/rss2.5/pfss_omni_ballistic_matched_rows.csv
    data/processed/comparison/rss3.0/pfss_omni_ballistic_matched_rows.csv

If those files need to be regenerated, run:

    python scripts/17_build_pfss_omni_ballistic_match_rss2p0.py
    python scripts/17_build_pfss_omni_ballistic_match_rss2p5.py
    python scripts/17_build_pfss_omni_ballistic_match_rss3p0.py
    python scripts/28_final_sample_metrics.py

## Requirements

Install dependencies with:

    python -m pip install -r requirements.txt

The core packages are listed in `requirements.txt`.

## Limitations

This project should be interpreted as a sensitivity study, not as a predictive solar-wind model.

| Limitation | Effect |
|---|---|
| PFSS assumes a current-free corona | coronal magnetic structure is simplified |
| Source-surface height is imposed | results depend on the chosen `rss` value |
| Ballistic mapping assumes simple longitudinal transport | transient structure and acceleration are not fully modeled |
| OMNI is near-Earth data | comparison is not a direct Parker Solar Probe in situ match |
| CR 2284 may include disturbed solar-wind structure | ambient and disturbed intervals can influence the correlations |
| IMF polarity uses a near-Earth radial-field proxy | polarity agreement should be treated as approximate |

## Final headline claim

For the locked final science sample, `rss = 2.5` gives the strongest association between equatorial unsigned PFSS magnetic-field strength and OMNI magnetic-field magnitude. The result persists when CR 2284 is removed, and `rss = 2.5` remains strongest in four of five leave-one-rotation-out tests.
