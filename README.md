> [!IMPORTANT]  
> **Before you begin — upstream package & prerequisites**
>
> This repository is **not** a replacement for the NOS skill-assessment workflow. It is an
> **equivalence checker** that consumes outputs from the official packages and verifies that the
> modern Python workflow reproduces the legacy Fortran results.
>
> - Start with the NOS package under active development:
>   **Next-Gen NOS OFS Skill Assessment** — https://github.com/NOAA-CO-OPS/Next-Gen-NOS-OFS-Skill-Assessment
>   
>   This package contains important background for understanding what this script does. Please default to that repository for any and all questions regarding the Python skill assessment software.  
>   -Run that repository end-to-end for your target OFS and date window.
>
> - (If you have access) also run the **legacy Fortran skill-assessment package** to produce its
>   corresponding outputs. *Note: the Fortran workflow is access-restricted to NOAA/NOS users.*
>   
>   Repository for the legacy/operational package: https://github.com/NOAA-CO-OPS/NOS-OFS-Skill-Assessment-Code
>
> - This script consumes:
>   - The **Python package outputs** (both `control_files/` and `data/` for nowcast/forecast)
>   - The **Fortran outputs** (WL/CU/Salt/Temp tables and raw time series)
>   
>   It then performs datum alignment, table normalization, statistical equivalence testing, and
>   agreement checks between the two workflows.
>
> - Folder names and file expectations are documented below in **Folder Layout & Required Inputs**.
>
> If you have not yet run the NOS package linked above, please do that first. Those outputs are the
> required inputs for the comparisons performed here.
> 
> **No access to both outputs?** You can still explore results. This repository includes example
> comparison workbooks and time series plots generated from both workflows so you can review the methodology
> and expected outputs without running the restricted Fortran package. Included results detail
> comparisons for cbofs from Jan-Jun 2024 and dbofs in Jan 2024.
>
> Please see the abstract and presentation files for a summary of this project.




# General Information / Capabilities

The AnalysisScript.py file is a script for verifying that a new, Python-based NOS skill-assessment package produces results statistically indistinguishable from the legacy Fortran package that is currently operational.

Given raw outputs from both packages (time series control files and skill assessment score tables) the script can:

- Convert Fortran water-level (WL) files from their native MSL reference to any tidal datum the Python run uses (MLLW, NAVD88, MLW, etc.).
- Transform Fortran skill-assessment tables into CSVs that match the Python column format exactly.
- Run a suite of equivalence tests (paired t, Wilcoxon signed-rank, two-sample KS, Cohen’s d, two-one-sided TOST) plus tolerance checks, at every station and across the OFS
- Produce a consolidated Excel workbook summarising pass/fail and confidence metrics for nowcast & forecast, each physical variable, and the total system.
- Interactively visualise model vs. observation behaviour in two ways:
  1. Time-series overlays (model vs obs, Python vs Fortran)
  2. Frequency domain Bode plots of magnitude & phase error
- Download the correct <ofs>_vdatums.nc grid automatically from NOAA’s public S3 bucket, so users never have to manage datum grids manually
- Work with any NOS OFS. Supply the OFS tag (e.g. cbofs, dbofs, sscofs) and the script infers all folder names and file prefixes.


## Pre-processing Stage

1. Interactive prompts ask the user for
   - OFS tag (cbofs, dbofs, …)
   - Start & end UTC (YYYY-MM-DD HH:MM)
   - Water-level datum used by the Python run  
     (MHHW MHW MTL MSL DTL MLW MLLW NAVD STND)
   - Whether to run the optional plotting and bode modules

2. `configure()` sets all global paths, downloads the `vdatum` grid if it is not already present,
   and creates a datum-specific folder like `fortran run/wl_mllw_shifted/`.

3. Fortran WL shift (optional)
   - KD-tree finds the nearest offset in `<ofs>_vdatums.nc` for each station.
   - Offsets (MSL→ user datum) are written to `station_offsets.csv`.
   - If the user agrees, every raw WL `.dat` file is rewritten with the offset applied and saved as
     `*_mllw.dat`, `*_navd.dat`, etc.

If the Python run already uses MSL, both the offset build and shift are skipped automatically.


## Comparison Stage

**1. Table Conversion**

`convert_var()` walks through every variable folder (WL, CU, Salt, Temp) and

- Reads each Fortran fixed-width skill-assessment table
- Extracts phase bias, observation mean, RMSE / bias / σ rows
- Writes a uniform CSV so the package loaders don’t need special parsing rules

**2. Mean-Accuracy Check (pass-ratio test)**

For each metric we first compare the mean values at each station, and across the OFS. If these values sit inside a relative tolerance band at each station, we pass a “Yes” equivalence flag.

Tolerance Bands Used in the Analysis:

| Metric category               | Tolerance band                 |
| ---                           | ---                            |
| Water temperature bias        | ± 1.00 °C                      |
| Water temperature RMSE        | ≤ 1.25 °C                      |
| Water level bias              | ± 0.12 m                       |
| Water level RMSE              | ≤ 0.18 m                       |
| RMSE, bias, SD                | ± 15 % of mean                 |
| Bias %                        | ± 15 percentage points         |
| Central-frequency             | model CF ≥ 85 % (± 5 pp window)|
| Positive / negative outlier freq | ≤ 1.5 × NOS 1 % limit (i.e. ≤ 1.5 %) |

If ≥ 50 % of stations pass then mean accuracy = YES, else NO.

**3. Statistical Equivalence (composite test)**

For every metric (RMSE, bias, bias %, bias SD, r, central freq, outlier freq) we run five statistical tests (α = 0.05):

| Test                         | Null hypothesis                                         | Pass criteria |
| ---                          | ---                                                     | ---           |
| Two-one-sided (TOST)         | Model and obs differ by more than a defined relative tolerance | p < α         |
| Distribution test (paired t and Wilcoxon) | Mean difference ≠ 0                           | p > α         |
| Shape test (two-sample KS)   | Model and obs drawn from different distributions        | p > α         |
| Size test (Cohens D)         | Size of the difference between the two sets is significant | p < 0.3     |

A weighted vote (TOST×3 + 4 other tests) produces a final equivalent = YES/NO flag for each statistical variable.

**4. Pass/Fail Agreement Verification**

After computing equivalence separately for Python and Fortran, we verify that the two packages agree by running McNemar’s test on the 2 × 2 contingency table of (pass, fail) outcomes.

| Outcome | Interpretation                                   |
| ---     | ---                                              |
| p > α   | Python & Fortran are statistically consistent    |
| p ≤ α   | Packages disagree more often than random chance  |

Agreement status is reported in the “agreement” column, and the p-value is stored in `p_mcnemar`. This generates a third set of equivalent = YES/NO flags.

All three test results (mean accuracy, statistical equivalence, binary verification) are consolidated into `compare_<ofs>_<year>.xlsx` with dedicated sheets for nowcast, forecast, and a summary of the overall comparison.

## Visual Diagnostics

Interactive Time-Series Plots (make_plots)

Raw time series overlays for every station and every variable:

- Pyt-model vs. Pyt-obs
- Frt-model vs. Frt-obs (after datum shift)
- Pyt-model vs. Frt-model

Saved under plots

Bode Plots (make_bode_all)

- Log-binned magnitude & phase error in the frequency domain.
- Plotted for each station plus one “ALL” curve that aggregates the entire OFS.
  - Inaccuracies may occur here due to stations being physically distant
 


## Folder Layout & Required Inputs

This script was designed to be run with data in the following folders. Either use this convention, or change the directory definitions within the script to match your configuration. Within your working directory create folders “fortran run” and “python run”. In the python folder please create “nowcast” and “forecast” folders, then copy and paste the respective python package outputs (both the “control_files” and “data” folders) into these folders. Within the fortran folder please create 4 folders, one for each variable, named “ofs” “year” CU, “ofs” “year” WL, “ofs” “year” Salt, and “ofs” “year” Temp; where ofs is the tag for the OFS run at (like cbofs, dbofs, sscofs, etc) and year is the year the package ran in (2024, 2022, etc). Example shown below. Then, for each of your runs on the server (for each variable) please copy the entirety of your “work” folder into its respective folder in your working directory.
In the “python run” folder.

<img width="580" height="115" alt="image" src="https://github.com/user-attachments/assets/6a546428-3d09-4024-ad9f-b467169efa06" />
<img width="395" height="93" alt="image" src="https://github.com/user-attachments/assets/ff2aaff2-0300-4e59-ba80-e4e3e3f89636" />

In the “fortran run” folder.

<img width="443" height="165" alt="image" src="https://github.com/user-attachments/assets/3c7c2237-8b29-4a6c-99f7-3062f51c3786" />
<img width="225" height="245" alt="image" src="https://github.com/user-attachments/assets/38f75964-1140-45de-a7b8-534c107eddbe" />



Please make sure the provided script is copied to your working directory.
```text
working_dir/
|
├─ fortran run/
|  ├─ <ofs> <year> WL/            ← raw WL *.dat files (full “work” folder from fortran)
|  ├─ <ofs> <year> CU/
|  ├─ <ofs> <year> Salt/
|  └─ <ofs> <year> Temp/
|
├─ python run/
|  ├─ nowcast/
|  |  ├─ control_files/
|  |  └─ data
|  └─ forecast/
|     ├─ control_files/
|     └─ data
|
└─ AnalysisScript.py
```

## Running the Script

Run the below commands to download all required python packages and run the script.

```bash
pip install numpy pandas scipy statsmodels matplotlib pyproj xarray netCDF4 openpyxl
cd /workingdirectory
python AnalysisScript.py
```
Interactive prompts will guide you through OFS, date window, datum, and which optional modules to run.


## Example Results (browse/download)
You will need to download the raw html file to view any time series plots, as they are interactively generated with plotly. 

- [`results/cbofs/`](results/cbofs/) — workbook + plots for CBOFS (2024)
- [`results/dbofs/`](results/dbofs/) — workbook + plots for DBOFS (2024)


## Output Format

|                     |                                           |
| ---                 | ---                                       |
| fortran_csv/        | converted Fortran skill tables            |
| wl_<datum>_shifted/ | datum-converted WL *.dat files            |
| station_coords.csv  | lat/lon lookup table                      |
| station_offsets.csv | MSL→ user datum offsets (m)               |
| vdatums.nc          | downloaded datum grid for this OFS        |
| compare_<ofs>_<year>.xlsx | ← MAIN summary workbook            |
| plots/              | Interactive time-series plots             |
| bode/               | Bode plots                                |

## Module Reference

| Group        | Direct Use        | Internal helpers |
| ---          | ---               | --- |
| **Shift**    | `shift_fortran_wl` | `build_station_coords`, `make_station_offsets`, `load_station_offsets` |
| **Convert Tables** | `convert_var` | `parse_fortran_table`, `parse_phase_bias`, `parse_obs_mean` |
| **Comparison** | `comparison_main` | `analyse`, `extract_tost_pvals`, `build_summary`, `add_accuracy`, `csv_path`, `norm_id` |
| **Plotters** | `make_plots`, `make_bode_all` | `build_file_index`, `_plot_triplet`, `_station_bode_plot`, `_log_bin`, `parse_fixed_width` |
| **Driver**   | `run_all`         | — |

## Troubleshooting Pointers

| Symptom                         | Likely cause |
| ---                             | --- |
| `wl_<datum>_shifted` is empty   | No raw fortran WL folder found. Verify `FR_WL_RAW` path |
| Plot folders created but are empty | The `build_file_index()` regex didn’t match your file naming scheme. Confirm filenames contain `_{OFS_ABBR}_nowcast.dat` / `forecast.dat` or tweak the pattern. |
| Workbook shows “equivalent = NO” everywhere | Python and Fortran date windows don’t overlap. Check `start` / `end` prompts. |
| `ValueError: datum not supported` | User entered datum not in `DATUM_VAR`. The accepted list is printed in the prompt. |




