# dosefitr -- NanoBRET & Cell Viability Dose-Response Analysis

`dosefitr` is an R package for end-to-end dose-response analysis of **NanoBRET kinase binding assays** and **cell viability experiments**. Starting from raw BMG PHERAstar plate-reader exports, it computes BRET ratios, removes outliers (ROUT method), and fits three- or four-parameter logistic models to derive IC50 and pIC50 values with confidence intervals. It also provides automatic quality control metrics (Z-prime, assay window, CV%), batch processing across multiple plates, multi-compound overlay and plate-comparison plots, Excel export, and Scarab-formatted tables for NanoBRET kinase profiling data submission.

---

## Installation

```r
remotes::install_github("AKKDataAnalysis/dosefitr")
```

**Dependencies:** `openxlsx`, `dplyr`, `ggplot2`, `scales`

---

## Overview

The package follows a linear pipeline. Each function feeds into the next:

```
Raw Excel files
      |
      v
batch_ratio_analysis()        <- BRET: reads plates, computes BRET ratios
batch_viability_analysis()    <- Viability: reads plates, normalises signal
      |
      v
rout_outliers_batch()         <- detects and removes outliers (optional)
      |
      v
batch_drc_analysis()          <- fits dose-response curves (3PL or 4PL)
      |
      |---> batch_save_all_drc_plots()    <- saves one plot per compound
      |---> plot_multiple_compounds()     <- overlays selected compounds
      |---> compare_plates_drc()          <- compares same compound across plates
      |---> scarab_table()               <- generates Scarab-format export (NanoBRET)
      `---> scarab_viability()           <- generates Scarab-format export (viability)
```

---

## Input Files

### Directory layout

`batch_ratio_analysis()` expects all files for an experiment to live in the same directory (the working directory by default):

```
experiment_folder/
|-- info_tables.xlsx        <- metadata file (one sheet per plate)
|-- raw_data_1.xlsx         <- plate-reader export for plate 1
|-- raw_data_2.xlsx         <- plate-reader export for plate 2
`-- raw_data_3.xlsx         <- plate-reader export for plate 3
```

The function matches each sheet in `info_tables.xlsx` to a raw data file by the trailing number in the filename (e.g. sheet `Sheet1` or `Plate_1` -> `raw_data_1.xlsx`).

---

### `info_tables.xlsx` -- the metadata file

This is the only file you need to create manually. It is an Excel workbook with **one sheet per plate**. Each sheet name must end with a number (e.g. `Sheet1`, `Plate_1`, `2`) that matches the trailing number in the corresponding raw data filename.

#### Required columns

Each sheet must have **at least 4 columns** in this exact order (column names can be anything, but position matters):

| Column | Position | Content | Example |
|---|---|---|---|
| log(inhibitor) | 1 | Log10 inhibitor concentrations, one per row | `-9`, `-8.5`, `-8`, ... |
| Plate_Row | 2 | Plate row letter(s) assigned to this compound | `A`, `B`, `C`, ... |
| Construct | 3 | Protein/kinase name | `KinaseA`, `KinaseB` |
| Compound | 4 | Compound/inhibitor name | `CpdA`, `CpdB` |

#### How the columns are used

- **Column 1 (log concentrations):** Becomes the first column of the final ratio table (`log(inhibitor).[M]`). Each row corresponds to one concentration point.
- **Column 2 (Plate_Row):** Maps each row of the info table to a row of the plate (A, B, C, ...). This is how the function knows which plate rows belong to which compound.
- **Columns 3 & 4 (Construct + Compound):** Combined as `Construct:Compound` to label the output columns (e.g. `KinaseA:CpdA`). This label is used throughout the pipeline -- in plots, summary tables, and the Scarab export.

#### Example sheet layout

A full 384-well plate (rows A-P) with one compound per row and three constructs.
The log(inhibitor) column lists the concentrations used in the experiment -- one
per row, in any order. These values are not tied to a specific plate row; the
function reads them as a sequence and assigns them to the data columns in order.

| log(inhibitor) | Plate_Row | Construct | Compound |
|---|---|---|---|
|  | A | KinaseA | CpdA |
| -4.61 | B | KinaseA | CpdB |
| -5.13 | C | KinaseA | CpdC |
| -5.61 | D | KinaseA | CpdD |
| -6.13 | E | KinaseA | CpdE |
| -6.61 | F | KinaseA | CpdF |
| -7.13 | G | KinaseA | CpdG |
| -7.61 | H | KinaseB | CpdH |
| -8.00 | I | KinaseB | CpdI |
| -8.43 | J | KinaseB | CpdJ |
| -8.91 | K | KinaseB | CpdK |
| -9.38 | L | KinaseC | CpdL |
|  | M | KinaseC | CpdM |
|  | N | KinaseC | CpdN |
|  | O | KinaseC | CpdO |
|  | P | KinaseC | CpdP |

> Each plate row letter (Plate_Row column) maps to a physical row on the plate.
> The log(inhibitor) column is a list of concentrations used in the experiment
> and is not tied to a specific plate row -- it is read as a sequence.

#### Duplicate constructs (biological replicates)

If the same `Construct:Compound` combination appears more than once in the sheet (e.g. two independent transfections of the same kinase), the function automatically appends a suffix to distinguish them:

| Construct | Compound | -> Output column label |
|---|---|---|
| KinaseA | CpdA | `KinaseA:CpdA` |
| KinaseA | CpdA | `KinaseA_2:CpdA` |

#### Multiple plates

Add one sheet per plate to the same `info_tables.xlsx` file. Sheet names must end with the plate number:

```
info_tables.xlsx
|-- Sheet1   <- matches raw_data_1.xlsx
|-- Sheet2   <- matches raw_data_2.xlsx
`-- Sheet3   <- matches raw_data_3.xlsx
```

Alternatively, use any naming convention as long as the trailing number matches:

```
info_tables.xlsx
|-- Plate_1   <- matches experiment_plate_1.xlsx
`-- Plate_2   <- matches experiment_plate_2.xlsx
```

#### Changing the info file name or location

By default the function looks for `info_tables.xlsx` in the working directory. Both can be changed:

```r
batch_ratio_analysis(
  directory = "path/to/experiment",
  info_file = "my_plate_layout.xlsx"
)
```

### Raw plate-reader Excel files

Each raw data file is a standard Excel export from the plate reader. The function reads only the first sheet.

#### BMG PHERAstar format (primary)

The PHERAstar exports two emission-channel measurement tables in a single sheet, separated by a few header rows. The function locates them by searching column 2 for rows matching the pattern `"N. Raw Data ..."` (e.g. `"1. Raw Data (450-80 B)"`).

The expected layout after the instrument header block is:

```
Row  1-N  |  Instrument metadata, plate ID, protocol info  (ignored)
          |
          |  -- First emission channel ------------------------------------------
Row  ?    |  "1. Raw Data (450-80 B)"   <- title row detected in col 2
Row  ?+1  |   [blank]  1   2   3  ...  24   <- column-number header
Row  ?+2  |   A        [well values]         <- data rows (rows A-H or A-P)
  ...     |   B, C, D, ...
          |
          |  -- Second emission channel ------------------------------------------
Row  ?    |  "2. Raw Data (610-20 B)"   <- second title row
Row  ?+1  |   [blank]  1   2   3  ...  24   <- column-number header
Row  ?+2  |   A        [well values]
  ...     |   B, C, D, ...
```

The function automatically finds the column-number header row (the row immediately after each `"N. Raw Data"` title) and slices the data from there, so the exact row numbers in the file do not matter.

After the column-number header is located, the function internally re-indexes the data as follows -- these are the row positions the ratio functions operate on:

| Rows (internal) | Content |
|---|---|
| Row 9 | Column names (well numbers `1, 2, ..., 24`) |
| Rows 10-25 | **Subtable 1** -- first emission channel (e.g. donor / luciferase, 16 rows = rows A-P of a 384-well plate) |
| Rows 26-27 | Separator / second channel header (ignored) |
| Rows 28-43 | **Subtable 2** -- second emission channel (e.g. acceptor / NanoBRET dye, 16 rows) |

> For a **96-well plate** (rows A-H, 8 rows), the same structure applies but only 8 data rows are present per subtable. The function auto-detects plate format from the number of column-header integers (<=12 columns -> 96-well; >12 columns -> 384-well). Use `plate_format = "96"` or `plate_format = "384"` to override if auto-detection fails.

The BRET ratio is computed as:

```
BRET ratio = (Subtable 2 / Subtable 1) x 1000
```

Wells where the donor-channel value (Subtable 1) falls below `low_value_threshold` (default `1000`) are set to `NA` before the ratio is calculated, to exclude failed or empty wells.

#### Generic / other plate readers (fallback)

If no `"Raw Data"` title row is found in column 2, the function falls back to a generic detection strategy: it scans the first 30 rows for a row whose integer values form a consecutive sequence starting at 1 (e.g. `1, 2, ..., 12` or `1, 2, ..., 24`). This makes the function compatible with exports from other instruments (e.g. Tecan Spark, BioTek Synergy) as long as they export a standard plate layout with column numbers as a header row.

#### What the function validates

Before passing data to the ratio functions, the parser checks that:

- Both emission-channel data blocks are present
- Row labels in each block match the expected plate row letters (A-H for 96-well, A-P for 384-well)

If validation fails, a descriptive error message is shown indicating which block and which row labels are unexpected.

#### Common issues

| Problem | Likely cause | Fix |
|---|---|---|
| `"Could not locate the plate column-number header"` | File is not a standard plate-reader export, or the sheet is empty | Check that the correct sheet is being exported and the file is not corrupted |
| `"Row label validation failed"` | Row letters in the data block are missing or out of order | Check that the plate-reader export includes row labels (A, B, C, ...) in the first column of each data block |
| `"Data must have at least 43 rows"` | File has too few rows after the header | Ensure both emission channels are exported; partial exports (one channel only) are not supported |
| Auto-detection picks wrong plate format | 384-well plate exported with only 12 columns selected | Pass `plate_format = "384"` explicitly |

---

## Workflow 1 -- NanoBRET Kinase Binding Assay

### Step 1 -- Process raw plates

```r
library(dosefitr)

results <- batch_ratio_analysis(
  control_0perc   = 24,    # column position of the 0% control (e.g. DMSO)
  control_100perc = 12,    # column position(s) of the 100% control (e.g. staurosporine)
  output_dir      = "./qc_output",
  verbose         = TRUE,
  low_value_threshold = 1000   # donor-channel values below this are set to NA
)
```

`batch_ratio_analysis()` reads all raw BMG PHERAstar Excel files in the working directory, pairs them with the `info_tables.xlsx` metadata file, computes BRET ratios, and returns a named list -- one entry per plate (e.g. `results$plate_01`, `results$plate_02`).

**Key parameters:**

| Parameter | Description |
|---|---|
| `control_0perc` | Column index or name for the 0% inhibition control |
| `control_100perc` | Column index(es) or name for the 100% inhibition control |
| `split_replicates` | Split technical duplicates into separate `.2` columns (default `TRUE`) |
| `low_value_threshold` | Donor-channel threshold below which values are treated as NA (default `1000`) |
| `function_version` | `"v1"` (default) or `"v2"` -- see below |
| `generate_reports` | Save per-plate and consolidated Excel reports (default `TRUE`) |

**Verbose output (`verbose = TRUE`):**

When `verbose = TRUE`, the function prints a concise progress log. The settings block and redundant path messages have been removed -- only actionable information is shown:

```
Found 4 plate sheet(s): plate_01, plate_02, plate_03, plate_04

Processing plate_01 -- 260323_LRRK2_plate1_Raw_01.xlsx
  [read] found 2 emission table(s) in '260323_LRRK2_plate1_Raw_01.xlsx'
  Saved: results_01.xlsx
v plate_01 done

Processing plate_02 -- 260323_LRRK2_plate1_2_Raw_02.xlsx
  [read] found 2 emission table(s) in '260323_LRRK2_plate1_2_Raw_02.xlsx'
  Saved: results_02.xlsx
v plate_02 done

...

=== BATCH COMPLETE ===
Plates processed: 4
Report saved:     ./drc_quality/batch_analysis_report.xlsx
======================
```

The `[read fallback]` line only appears when the file is not a standard PHERAstar export and the generic column-header detection path is used. Row-label validation is silent on success and only speaks up on failure.

---

#### `function_version`: v1 vs v2

`batch_ratio_analysis()` supports two internal ratio-calculation engines, selected via `function_version`. **v1 is the default** and is suitable for most experiments. v2 adds support for fixed-value 0% controls and multiple 100% control columns.

| | v1 (default) | v2 |
|---|---|---|
| `control_0perc` | Column **position** or column **name** (e.g. `24` or `"DMSO"`) | Column position, column name, **or** a fixed numeric value applied to all wells |
| `control_100perc` | **Single** column position or name (e.g. `12` or `"Staurosporine"`) | **One or more** column positions or names (e.g. `c(12, 24)`) |
| Multiple 100% controls | Not supported -- only one column | Supported -- averages them into a single `Mean_100perc` column |
| Fixed 0% control | Not supported | Supported -- creates a `Fixed_0perc` column with a constant value |
| Default `low_value_threshold` | `1000` | `3000` |
| Use when | One 0% control column and one 100% control column | Duplicate 100% control columns that need averaging, or a fixed background value for 0% |

**v1 -- single-control mode (default)**

Both controls are a single column, identified by position (most common) or by name. Column positions are 1-indexed from the first data column (i.e. column 1 in the plate reader export = position 1).

```r
# Column index (most common -- use the column number in the plate export)
results <- batch_ratio_analysis(
  control_0perc   = 24,   # column 24 is the 0% control (e.g. DMSO)
  control_100perc = 12,   # column 12 is the 100% control (e.g. staurosporine)
  verbose         = TRUE
  # function_version = "v1" is the default -- no need to specify
)

# Column name (must match the header row exactly as it appears in the file)
results <- batch_ratio_analysis(
  control_0perc   = "DMSO",
  control_100perc = "Staurosporine",
  verbose         = TRUE
)
```

**v2 -- fixed-value / multi-column mode**

Use v2 when you have **duplicate 100% control columns** that should be averaged before normalisation, or when your 0% control is a fixed instrument background value rather than a measured well. When multiple positions are given for `control_100perc`, their values are averaged row-wise into a single `Mean_100perc` column; the original control columns are then removed from the data.

```r
results <- batch_ratio_analysis(
  control_0perc       = 16,        # fixed value used as 0% control for all wells
  control_100perc     = c(12, 24), # columns 12 and 24 are duplicate 100% controls -- averaged
  output_dir          = "./drc_quality",
  verbose             = TRUE,
  low_value_threshold = 3000,
  function_version    = "v2"
)
```

Both versions produce the same output structure and feed identically into all downstream functions (`rout_outliers_batch()`, `batch_drc_analysis()`, etc.).

#### `label_sep`: configurable construct/compound separator

By default, construct and compound names are joined with `":"` (e.g. `EPHA1:KK135`). If your construct or compound names already contain colons, choose a different separator:

```r
results <- batch_ratio_analysis(
  control_0perc   = 24,
  control_100perc = 12,
  label_sep       = "/"   # columns will be named e.g. "EPHA1/KK135"
)
```

The separator is stored as an attribute on the returned list and propagated automatically through the entire pipeline — `rout_outliers_batch()`, `batch_drc_analysis()`, `plot_multiple_compounds()`, `scarab_table()`, etc. You only need to set it once. The same parameter is available in `batch_viability_analysis()`.

| Value | Column label example | When to use |
|---|---|---|
| `":"` (default) | `EPHA1:KK135` | Standard — use unless names contain colons |
| `"/"` | `EPHA1/KK135` | Construct or compound names contain `:` |
| `"\|"` | `EPHA1\|KK135` | Either name contains `/` |

---

### Step 2 -- Remove outliers (optional)

```r
results_clean <- rout_outliers_batch(results, Q = 0.01)

# Inspect what was flagged
outliers  <- results_clean$outlier_summary
rescued   <- results_clean$rescued_summary

# Access the cleaned and original tables
original_table <- results_clean$plate_01$result$modified_ratio_table_original
clean_table    <- results_clean$plate_01$result$modified_ratio_table

# Visualise outlier curves
plot_outliers_batch_curves(results_clean)
```

`rout_outliers_batch()` applies ROUT outlier detection to each plate.

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `Q` | `0.01` | ROUT false-discovery rate threshold. Lower = more conservative (fewer outliers flagged). Typical range: `0.001` (strict) to `0.05` (lenient). |
| `n_param` | `4L` | Hill-model parameters: `4` (free Hill slope, default) or `3` (Hill slope fixed at -1). Use `3` for sparse data. |
| `direction` | `"inhibition"` | Curve direction: `"inhibition"` (response decreases with concentration) or `"activation"` / `"agonist"` (response increases). |
| `ntry_retry` | `3L` | Number of random-restart attempts when the initial fit does not converge. Increase to `10` for difficult curves. |
| `min_dynamic_range` | `20` | Warn when the estimated dynamic range (%) is below this value. Flat curves with small dynamic range are harder to fit reliably. |
| `keep_cytotoxic` | `FALSE` | When `TRUE`, retains data points that resemble outliers but are part of a genuine cytotoxic response (sharp drop at high concentrations). Recommended for viability assays. |
| `seed` | `42L` | Integer random seed set before the compound loop, ensuring reproducible outlier calls on repeated runs. Set to `NULL` to disable. |
| `verbose` | `TRUE` | Print per-plate progress and a final summary. |

The ROUT method is described in: Motulsky HJ & Brown RE, BMC Bioinformatics 2006, 7:123 (doi:10.1186/1471-2105-7-123).

---

### Step 3 -- Fit dose-response curves

Pass whichever variable holds your processed data to `batch_results`. If you ran outlier removal in Step 2, pass `results_clean` -- **not** the original `results`. `rout_outliers_batch()` overwrites `modified_ratio_table` in-place with the cleaned data (the original is preserved as `modified_ratio_table_original`), so `batch_drc_analysis()` will automatically use the outlier-removed values.

```r
# If you ran rout_outliers_batch() -- use results_clean
drc_results <- batch_drc_analysis(
  batch_results = results_clean,   # <- cleaned data, outliers removed
  normalize     = TRUE,
  output_dir    = "./drc_results",
  verbose       = TRUE
)

# If you skipped outlier removal -- use the original results
drc_results <- batch_drc_analysis(
  batch_results = results,         # <- original data, no outlier removal
  normalize     = TRUE,
  output_dir    = "./drc_results",
  verbose       = TRUE
)
```

**Key parameters:**

| Parameter | Description |
|---|---|
| `normalize` | Normalise responses to 0-100% before fitting |
| `model` | `"3pl"` (default, Hill slope fixed at +/-1) or `"4pl"` (Hill slope freely estimated) |
| `r_sqr_threshold` | Minimum R2 to accept a fit (default `0.8`) |
| `enforce_bottom_threshold` | Reject fits where the bottom plateau exceeds a threshold |
| `nd_if_activation` | Set IC50 and pIC50 to `"N/D"` for activation curves (default `FALSE`) -- see below |
| `generate_reports` | Save a consolidated `batch_drc_analysis_report.xlsx` |

**Choosing between 3PL and 4PL:**
- **3PL** (`model = "3pl"`): Hill slope is fixed at -1. Recommended when you expect a standard sigmoidal response and want fewer free parameters.
- **4PL** (`model = "4pl"`): Hill slope is freely estimated. Use when the curve shape is expected to deviate from a standard Hill slope (e.g. cooperative binding, mixed mechanisms).

**IC50 above the tested concentration range:**

When the fitted IC50 exceeds the highest concentration tested, the `IC50 (uM)` and `IC50 (nM)` columns in the `Pharmacology_Summary` table display a `>X` value (e.g. `>25`) rather than the extrapolated number. The compound is also flagged in the `Exclusion` column as `"IC50 above tested range (>X uM)"`. The pIC50 is still reported as a numeric value in this case.

**N/D for non-inhibitory curves (`nd_if_activation`):**

Compounds whose fitted curve is classified as **flat** (no meaningful response across the concentration range) always have their `IC50 (uM)`, `IC50 (nM)`, and `pIC50` set to `"N/D"` (not determined) in the `Pharmacology_Summary` table. This is the default behaviour and cannot be disabled.

Compounds whose curve goes **up** (activation -- response increases with concentration) are treated the same way when `nd_if_activation = TRUE`. The default is `FALSE`, which reports the fitted IC50 for activation curves as usual.

| `curve_type` | `nd_if_activation = FALSE` (default) | `nd_if_activation = TRUE` |
|---|---|---|
| `"flat"` | N/D | N/D |
| `"activation"` | IC50/pIC50 reported | N/D |
| `"inhibition"` | IC50/pIC50 reported | IC50/pIC50 reported |

```r
# Default -- only flat curves get N/D
drc_results <- batch_drc_analysis(batch_results = results, normalize = TRUE)

# Also suppress IC50 for activation curves
drc_results <- batch_drc_analysis(batch_results = results, normalize = TRUE,
                                  nd_if_activation = TRUE)
```

---

### Step 4 -- Generate Scarab export table

```r
table <- scarab_table(
  results_list      = results,
  drc_results_list  = drc_results,
  plate_name        = "plate_01",
  date              = "260323",
  experimenter_abbrev = "TL",
  nLuc_orientation  = "C",
  tracer            = "Tracer K10",
  tracer_kd_app     = -7.5,
  tracer_concentration_used = -7.5,
  decimal_separator = ",",   # use "." for English format
  eubscarab_ready   = "No"   # or "Yes", or a named/unnamed vector per compound
)
```

`scarab_table()` generates a 76-row Scarab-format data frame combining experimental metadata, curve parameters, QC metrics, and raw BRET values. It is saved as an `.xlsx` file with two sheets: **`Scarab_Table`** (fields as rows, compounds as columns -- the standard Scarab layout) and **`Scarab_Table_Transposed`** (compounds as rows, fields as columns -- easier for filtering and sorting in Excel).

**Per-kinase parameters -- three ways to specify:**

`nLuc_orientation`, `tracer`, `tracer_kd_app`, and `tracer_concentration_used` all support the same flexible input format:

```r
# Option 1: Single value -- applies to all kinases
scarab_table(..., nLuc_orientation = "C", tracer_kd_app = -7.5)

# Option 2: Named vector -- specify per kinase by name
scarab_table(...,
  nLuc_orientation  = c(KinaseA = "C", KinaseB = "N"),
  tracer_kd_app     = c(KinaseA = -7.5, KinaseB = -8.0),
  tracer            = c(KinaseA = "Tracer K10", KinaseB = "Tracer 236")
)

# Option 3: Unnamed vector -- applied in order of unique kinases in the plate
# (if unique kinases are KinaseA, KinaseB in that order)
scarab_table(...,
  nLuc_orientation  = c("C", "N"),
  tracer_kd_app     = c(-7.5, -8.0)
)
```

The Construct ID is generated automatically:
- `nLuc_orientation = "N"` -> `{Kinase}A-nb001`
- `nLuc_orientation = "C"` -> `{Kinase}A-nb002`

**`eubscarab_ready` -- EUbScarab readiness flag:**

Controls the value of the **"Is EUbScarab Ready?"** row in the output table. Accepts the same three input formats as the per-kinase parameters above:

```r
# Single value -- applies to all compounds (default: "No")
scarab_table(..., eubscarab_ready = "No")

# Named vector -- specify per compound (Kinase:Compound label)
scarab_table(..., eubscarab_ready = c("KinaseA:CpdA" = "Yes", "KinaseB:CpdB" = "No"))

# Unnamed vector -- applied in order of compounds in the plate
scarab_table(..., eubscarab_ready = c("Yes", "No", "Yes"))
```

If a compound is not found in a named vector, or the unnamed vector is shorter than the number of compounds, the value defaults to `"No"`.

---

## Workflow 2 -- Cell Viability Assay

The viability workflow mirrors the NanoBRET pipeline exactly. Replace `batch_ratio_analysis()` with `batch_viability_analysis()` in Step 1. `batch_drc_analysis()` automatically detects the assay type from its input -- no extra parameter needed. All downstream functions (plotting, table export) work identically.

### Step 1 -- Process raw plates

```r
via_results <- batch_viability_analysis(
  directory       = "data/viability/",
  control_0perc   = 1,
  control_100perc = 24,
  output_dir      = "./qc_output",
  verbose         = TRUE
)
```

`batch_viability_analysis()` reads all raw viability Excel files in the target directory, pairs them with the `info_tables.xlsx` metadata file, and returns a named list -- one entry per plate -- in the same format as `batch_ratio_analysis()`, so all downstream functions accept it without modification.

**Key parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `directory` | character | `getwd()` | Directory containing raw data files and the info file |
| `control_0perc` | integer (1-24) | `NULL` | Plate column index for the 0% viability control (background / vehicle) |
| `control_100perc` | integer (1-24) | `NULL` | Plate column index for the 100% viability control (untreated / positive control) |
| `split_replicates` | logical | `TRUE` | Split technical replicates into separate `.2` columns |
| `low_value_threshold` | numeric | `0` | Values below this are replaced with `NA` before processing (default `0` = no filtering) |
| `apply_control_means` | logical | `TRUE` | Replace individual control values with construct-specific means |
| `auto_detect` | logical | `TRUE` | Auto-detect the table layout inside each raw file |
| `info_file` | character | `"info_tables.xlsx"` | Name of the metadata Excel file (must be in `directory`) |
| `data_pattern` | character | `"_\d+\.xlsx$"` | Regex to identify raw data files; files starting with `~$` are always excluded |
| `output_dir` | character | `directory` | Directory for output files; created automatically if it does not exist |
| `generate_reports` | logical | `TRUE` | Save per-plate and consolidated Excel reports |
| `selected_columns` | integer vector | `NULL` | 1-based column indices to include (e.g. `2:23`); `NULL` uses all 24 columns |
| `verbose` | logical | `TRUE` | Print progress messages |

---

#### `control_0perc` and `control_100perc`

Both parameters accept a **single integer** between 1 and 24 (the 1-based column index in the plate layout). They are applied identically to every plate in the batch.

- **`control_0perc`** -- the column used as the 0% viability reference (background signal, e.g. a vehicle-only or cell-free well). This anchors the lower end of the normalisation.
- **`control_100perc`** -- the column used as the 100% viability reference (untreated cells, maximum signal). This anchors the upper end of the normalisation.

Both can be `NULL` if you do not want normalisation controls (quality metrics will not be computed in that case).

```r
# Typical setup: column 1 = background, column 24 = untreated cells
via_results <- batch_viability_analysis(
  control_0perc   = 1,
  control_100perc = 24
)

# Exclude outer columns and use inner controls
via_results <- batch_viability_analysis(
  control_0perc    = 13,
  control_100perc  = 12,
  selected_columns = 2:23   # exclude columns 1 and 24 from compound data
)
```

> **Note:** Unlike `batch_ratio_analysis()` v2, `batch_viability_analysis()` does not support multiple 100% control columns or fixed numeric 0% values. Each control must be a single column index.

---

#### `selected_columns`

Restricts which plate columns are included in the output table. Useful when the outer columns of the plate are reserved for controls and should not appear as compound data.

```r
# Include only columns 2-23 (exclude columns 1 and 24, which are controls)
via_results <- batch_viability_analysis(
  control_0perc    = 1,
  control_100perc  = 24,
  selected_columns = 2:23
)
```

`NULL` (default) includes all 24 columns.

---

#### `low_value_threshold`

Values in the raw data below this threshold are replaced with `NA` before processing. The default is `0`, meaning no filtering is applied. Increase this to exclude wells with near-zero signal that likely represent failed or empty wells.

```r
via_results <- batch_viability_analysis(
  control_0perc       = 1,
  control_100perc     = 24,
  low_value_threshold = 500   # exclude wells with signal < 500
)
```

> This parameter serves the same purpose as `low_value_threshold` in `batch_ratio_analysis()`, but the appropriate threshold depends on your instrument and assay signal levels. For viability assays the default of `0` is usually appropriate; for NanoBRET the default is `1000`.

---

#### `apply_control_means`

When `TRUE` (default), the raw values in the control columns are replaced with the **construct-specific mean** of that control across all replicates of the same construct. This reduces the effect of single-well noise on normalisation.

Set to `FALSE` to use the raw per-well control values directly.

---

#### `split_replicates`

When `TRUE` (default), technical replicates (multiple plate rows assigned to the same `Construct:Compound`) are split into separate columns in the output table, labelled with a `.2` suffix for the second replicate. This matches the behaviour of `batch_ratio_analysis()`.

---

#### `generate_reports` and output files

When `generate_reports = TRUE`, the function writes two types of Excel files into a `drc_quality/` sub-folder inside `output_dir`:

| File | Contents |
|---|---|
| `drc_quality/viability_results_<N>.xlsx` | Per-plate workbook with four sheets: `Quality_Metrics`, `Modified_Table`, `Original_Table`, `Processing_Info` |
| `drc_quality/batch_viability_report.xlsx` | Consolidated workbook with a `Summary` sheet (one row per plate) |

The `Quality_Metrics` sheet reports per-construct quality metrics based on the CV% of each control:

| Metric | Description |
|---|---|
| `Mean_Background` | Mean of 0% control replicates for this construct |
| `SD_Background` | Standard deviation of 0% control replicates |
| `CV_Background_pct` | CV% of 0% control (SD / Mean x 100) |
| `Mean_Positive_Ctrl` | Mean of 100% control replicates |
| `SD_Positive_Ctrl` | Standard deviation of 100% control replicates |
| `CV_Positive_Ctrl_pct` | CV% of 100% control |
| `CV_Background_Comment` | Quality label: `high (<=10%)`, `medium (10-20%)`, `low (>20%)` |
| `CV_PosCtrl_Comment` | Quality label for 100% control CV% |
| `Overall_Quality` | Lowest quality level across both controls |
| `Signal_to_Background` | Mean_Positive / Mean_Background (descriptive only -- not used in Overall_Quality) |
| `Rows` | Construct name and plate row range used |
| `Rows_Count` | Number of replicate rows |

> `Overall_Quality` is driven exclusively by CV% of the controls, which are scale-independent. `Signal_to_Background` is reported for information only and does not affect the quality assessment, because its absolute value depends on instrument signal levels and is not interpretable without context before normalisation.

---

#### Return value

A named list with one element per successfully processed plate. The element name is the info-sheet name (e.g. `"Sheet1"`). Each element contains:

| Field | Description |
|---|---|
| `data_file` | Filename of the raw data file used |
| `info_sheet` | Name of the info sheet used |
| `sheet_number` | Plate number extracted from the sheet name |
| `control_0perc` | Value passed as `control_0perc` |
| `control_100perc` | Value passed as `control_100perc` |
| `selected_columns` | Value passed as `selected_columns` |
| `result` | Output of `process_viability_data()`, with `$modified_table` renamed to `$modified_ratio_table` for downstream compatibility |

Plates that fail processing are omitted and a warning is issued. If no plates succeed, an empty list is returned with a warning.

---

### Step 2 -- Remove outliers (optional)

```r
results_clean <- rout_outliers_batch(via_results, Q = 0.01, keep_cytotoxic = TRUE)
```

Use `keep_cytotoxic = TRUE` to retain points that look like outliers but are part of a genuine cytotoxic response (sharp drop at high concentrations). See Workflow 1 Step 2 for full details.

---

### Step 3 -- Fit dose-response curves

```r
# If you ran rout_outliers_batch() -- use results_clean
drc_results <- batch_drc_analysis(
  batch_results = results_clean,   # <- cleaned data, outliers removed
  normalize     = TRUE,
  output_dir    = "./drc_results",
  verbose       = TRUE
)

# If you skipped outlier removal -- use the original results
drc_results <- batch_drc_analysis(
  batch_results = via_results,     # <- original data, no outlier removal
  normalize     = TRUE,
  output_dir    = "./drc_results",
  verbose       = TRUE
)
```

`batch_drc_analysis()` automatically applies the appropriate parameter limits based on the input source. See Workflow 1 Step 3 for full `batch_drc_analysis()` parameter details.

---

### Step 4 -- Generate Scarab export table (viability)

```r
scarab_viability(
  results_list        = via_results,
  drc_results_list    = drc_results,
  plate_name          = "plate_01",
  date                = "260323",
  experimenter_abbrev = "TL",
  cell_line           = "HeLa",
  cell_type           = "cervical adenocarcinoma",
  treatment_time      = "72h",
  measurement_method  = "CellTiter-Glo",
  assay_volume        = "40",
  plate_format        = "384-well",
  decimal_separator   = ","   # use "." for English format
)
```

`scarab_viability()` generates a Scarab-format data frame for cell viability experiments. It is the viability counterpart of `scarab_table()` and produces an `.xlsx` file with two sheets: **`Scarab_Table`** (fields as rows, compounds as columns) and **`Scarab_Table_Transposed`** (compounds as rows -- easier for filtering in Excel).

**Key parameters (viability-specific):**

| Parameter | Default | Description |
|---|---|---|
| `cell_line` | `"HeLa"` | Cell line used in the experiment. |
| `cell_type` | `"cervical adenocarcinoma"` | Cell type description. |
| `treatment_time` | `"72h"` | Duration of compound treatment (e.g. `"24h"`, `"48h"`, `"72h"`). |
| `measurement_method` | `"CellTiter-Glo"` | Viability measurement method (e.g. `"CellTiter-Glo"`, `"resazurin"`, `"MTT"`). |
| `assay_volume` | `"40"` | Assay volume in µL. |
| `plate_format` | `"384-well"` | Plate format string. |
| `plate_manufacturer` | `"3570"` | Plate manufacturer or catalog number. |
| `plate_material` | `"PS"` | Plate material. |
| `sgc_compound_id` | `NA` | SGC global compound identifier / batch number. |
| `eln_id` | `NA` | Electronic Lab Notebook ID. |
| `comments` | `NA` | Free-text comments. |
| `decimal_separator` | `"."` | Decimal separator for numeric values in the output: `"."` (English) or `","` (European). |
| `save` | `TRUE` | Save the output as an `.xlsx` file. |
| `file_name` | `NULL` | Custom filename; `NULL` auto-generates a name from the plate name and date. |

---

## Merging Replicate Plates (optional)

When the same experiment is run across multiple plates (biological replicates), `merge_plate_replicates()` combines their `modified_ratio_table`s into a single merged entry that feeds directly into all downstream functions.

### How it works

Each plate from `batch_ratio_analysis()` already contains two replicate columns per compound (e.g. `KinaseA:CpdA` and `KinaseA:CpdA.2`). The merge function collects all replicate columns across plates and renumbers them sequentially:

| Source | Columns in plate | Columns after merge |
|---|---|---|
| Plate 1 | `KinaseA:CpdA`, `KinaseA:CpdA.2` | `KinaseA:CpdA`, `KinaseA:CpdA.2` |
| Plate 2 | `KinaseA:CpdA`, `KinaseA:CpdA.2` | `KinaseA:CpdA.3`, `KinaseA:CpdA.4` |
| Plate 3 | `KinaseA:CpdA`, `KinaseA:CpdA.2` | `KinaseA:CpdA.5`, `KinaseA:CpdA.6` |

The individual plate entries are removed from the list and replaced by a single `"merged"` entry. The merged result is a drop-in replacement -- all downstream functions (`rout_outliers_batch()`, `batch_drc_analysis()`, `plot_multiple_compounds()`, etc.) accept it without modification.

### Basic usage

```r
# Run batch_ratio_analysis as normal
results <- batch_ratio_analysis(
  control_0perc   = 24,
  control_100perc = 12
)

# Merge all plates (default -- merges everything in results)
merged <- merge_plate_replicates(results)

# Feed directly into the DRC step
drc_results <- batch_drc_analysis(
  batch_results = merged,
  normalize     = TRUE,
  output_dir    = "./drc_results"
)
```

### Merge a subset of plates

```r
# Keep Sheet3 separate; merge only Sheet1 and Sheet2
merged <- merge_plate_replicates(results,
  plates      = c("Sheet1", "Sheet2"),
  merged_name = "plates_1_2"
)
# merged$plates_1_2  <- the combined entry
# merged$Sheet3      <- kept as-is
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `results` | list | -- | Named list from `batch_ratio_analysis()` or `batch_viability_analysis()` |
| `plates` | character vector | `NULL` | Plate names to merge; `NULL` merges all plates |
| `merged_name` | character | `"merged"` | Name of the new merged entry in the returned list |
| `check_concentrations` | logical | `TRUE` | Stop with an error if the `log(inhibitor)` columns differ between plates |
| `output_dir` | character | `NULL` | Directory for the Excel report; `NULL` writes to `drc_quality/` in the working directory |
| `generate_reports` | logical | `TRUE` | Save `merged_results_<merged_name>.xlsx` into `drc_quality/`; set `FALSE` to skip |
| `verbose` | logical | `TRUE` | Print a summary of compounds and replicate counts |

### Asymmetric plates

Compounds that appear in only some of the merged plates are included with however many replicates they have -- no NA-padding is added:

```r
# Sheet1 has KinaseA:CpdA + KinaseA:CpdB
# Sheet2 has KinaseA:CpdA only
merged <- merge_plate_replicates(results, plates = c("Sheet1", "Sheet2"))
# Result: KinaseA:CpdA gets 4 replicates, KinaseA:CpdB gets 2
```

### Concentration check

All plates being merged must share an identical `log(inhibitor)` column. If they differ, the function stops with a descriptive error:

```
Error: Concentration mismatch between plate 'Sheet1' and plate 'Sheet2'.
  Plate 'Sheet1' has 9 concentration points: -9, -8.5, -8, -7.5, -7 ...
  Plate 'Sheet2' has 9 concentration points: -9, -8.5, -8, -7.5, -6 ...
```

Set `check_concentrations = FALSE` to bypass this check (not recommended).

### Provenance

The merged entry records which plates it was built from in a `$merged_from` field:

```r
merged$merged$merged_from
#> [1] "Sheet1" "Sheet2" "Sheet3"
```

### Output file

When `generate_reports = TRUE` (default), a `drc_quality/` folder is created in the working directory (matching the convention of `batch_ratio_analysis()`) and the Excel workbook is saved inside it:

```
drc_quality/
`-- merged_results_<merged_name>.xlsx
    |-- Merged_Table        <- combined modified table with all replicates renumbered
    |-- Original_Sheet1     <- modified_ratio_table from Sheet1 before merging
    |-- Original_Sheet2     <- modified_ratio_table from Sheet2 before merging
    `-- Provenance          <- one row per plate: data file, n columns, n compounds; plus merge date
```

```r
# Default -- saves to drc_quality/ in the working directory
merged <- merge_plate_replicates(results)

# Save report to a custom directory
merged <- merge_plate_replicates(results, output_dir = "./my_reports")

# Skip the report entirely
merged <- merge_plate_replicates(results, generate_reports = FALSE)
```

---
## Plotting

### Save all individual compound plots

```r
batch_save_all_drc_plots(
  batch_drc_results    = drc_results,
  output_dir           = "DRC_Plots",
  verbose              = TRUE,
  y_axis_title         = "% Cell Survival",   # optional override; auto-detects "Cell Viability (%)" by default
  y_limits             = NULL,                # NULL auto-scales each plot; use c(0, 100) for a fixed scale
  organize_by          = "plate",             # "plate" (default) or "compound"
  format               = "png",              # "png", "pdf", "svg", "tiff", "eps"
  width                = 10,
  height               = 10,
  dpi                  = 600,
  compounds_to_plot    = NULL,               # character vector to restrict to specific compounds
  plates_to_plot       = NULL,               # character vector to restrict to specific plates
  save_panel           = TRUE,               # also save a combined panel PNG per plate
  panel_ncol           = 4L,                 # columns in the panel grid
  panel_width_per_col  = 6,                  # inches per panel column
  panel_height_per_row = 6,                  # inches per panel row
  panel_spacing        = 1,                  # spacing between sub-plots (in cm)
  subplot_title        = "auto",             # "auto", "compound", "construct", or a fixed string
  show_ic50_line       = FALSE,
  point_color          = "black",
  point_size           = 2
)
```

Saves one PNG per compound across all plates into `output_dir`. When `save_panel = TRUE` (default), also saves a combined panel image per plate with all compounds arranged in a grid.

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `output_dir` | `"DRC_Plots"` | Directory for saved files; created automatically if it does not exist. |
| `organize_by` | `"plate"` | Sub-folder structure: `"plate"` creates one sub-folder per plate; `"compound"` creates one sub-folder per compound. |
| `format` | `"png"` | File format for individual compound plots. Any format accepted by `ggplot2::ggsave()`. |
| `compounds_to_plot` | `NULL` | Character vector of compound names to include; `NULL` plots all compounds. |
| `plates_to_plot` | `NULL` | Character vector of plate names to include; `NULL` plots all plates. |
| `save_panel` | `TRUE` | Save a combined panel image per plate in addition to individual plots. |
| `panel_ncol` | `4L` | Number of columns in the panel grid. |
| `panel_width_per_col` | `6` | Width in inches allocated to each panel column. |
| `panel_height_per_row` | `6` | Height in inches allocated to each panel row. |
| `panel_spacing` | `1` | Spacing between sub-plots in the panel (in cm). |
| `subplot_title` | `"auto"` | Title shown on each sub-plot: `"auto"` uses the compound label; `"compound"` shows only the compound name; `"construct"` shows only the construct name; any other string is used as a fixed title. |
| `y_axis_title` | `NULL` | Y-axis label. `NULL` auto-detects `"Cell Viability (%)"` for viability assays and `"Normalized BRET ratio [%]"` for NanoBRET. |
| `y_limits` | `NULL` | Y-axis limits as `c(min, max)`. `NULL` auto-scales each plot independently. |
| `show_ic50_line` | `FALSE` | Draw a vertical dashed line at the fitted IC50. |

---

### Overlay multiple compounds on one plot

```r
# Pass the full batch result -- first plate is used automatically
plot_multiple_compounds(drc_results, compound_indices = 1:5)

# Select a specific plate
plot_multiple_compounds(drc_results, plate = "plate_02", compound_indices = 9:11)

# Select compounds by name or partial match
plot_multiple_compounds(drc_results, target_compound = "KinaseA")
plot_multiple_compounds(drc_results, target_compound = "KinaseA:CpdA")

# Select a single compound by position
plot_multiple_compounds(drc_results, position = 3)
```

**Saving plots:**

```r
# Auto-named PNG
plot_multiple_compounds(drc_results, compound_indices = 1:5, save_plot = TRUE)

# Custom filename (extension determines format: .png, .pdf, .svg, .tiff, .eps)
plot_multiple_compounds(drc_results, compound_indices = 1:5,
                        save_plot = "my_compounds.png")
```

**Customising appearance:**

```r
plot_multiple_compounds(
  drc_results,
  plate            = "plate_01",
  compound_indices = c(1, 3, 5),
  color_palette    = "colorblind",   # see palette list below
  y_limits         = c(0, 100),      # or NULL to auto-scale to the data
  y_axis_title     = "% Normalized Luminescence",
  plot_title       = "GSK Compounds",
  axis_text_size   = 18,
  axis_title_size  = 20,
  legend_text_size = 14,
  legend_position  = "right",
  legend_ncol      = 1,
  show_error_bars  = TRUE,
  show_grid        = FALSE,
  shape_by_compound = FALSE,
  plot_width       = 12,
  plot_height      = 8,
  plot_dpi         = 600,
  save_plot        = TRUE
)
```

**Additional appearance controls** (new parameters):

```r
plot_multiple_compounds(
  drc_results,
  compound_indices       = 1:5,
  # Axis range
  x_limits               = c(-9, -5),   # log10 molar; or NULL to auto-scale
  x_limits_scale         = "log10",      # "log10" (default), "molar", "uM", "nM"
  x_axis_title           = NULL,         # NULL uses default Log10 Concentration [M]
  # Curve appearance
  curve_linewidth        = 1,            # fitted curve line width
  curve_alpha            = 0.7,          # fitted curve opacity (0-1)
  # IC50 reference lines
  show_ic50_lines        = TRUE,         # dashed vertical line at each IC50
  # Title and axis styling
  plot_title_size        = 16,           # NULL inherits axis_title_size + 2
  axis_line_color        = "black",      # colour of axis lines and ticks
  show_border            = FALSE,        # rectangular panel border
  # Background
  transparent_background = FALSE         # TRUE for transparent PNG export
)
```

**Available color palettes** (pass to `color_palette`):

| Category | Options |
|---|---|
| Default | `"hue"`, `"ggplot2"` |
| ColorBrewer qualitative | `"set1"`, `"set2"`, `"set3"`, `"dark2"`, `"paired"`, `"accent"` |
| ColorBrewer sequential | `"blues"`, `"reds"`, `"greens"`, `"purples"` |
| Colorblind-friendly | `"okabe_ito"`, `"colorblind"`, `"cud"`, `"tol"` |
| Scientific journals | `"nature"`, `"science"`, `"cell"`, `"plos"`, `"elife"` |
| Publishers | `"bmc"`, `"frontiers"`, `"wiley"`, `"elsevier"`, `"springer"` |
| Viridis | `"viridis"`, `"magma"`, `"inferno"`, `"plasma"` |
| Gradient | `"blue_red"`, `"green_red"`, `"cool_warm"` |
| Classic R | `"rainbow"`, `"heat"`, `"terrain"` |

You can also supply a custom color vector directly:

```r
plot_multiple_compounds(drc_results,
  colors = c("#E41A1C", "#377EB8", "#4DAF4A"))
```

---

### Compare the same compound across plates

```r
compare_plates_drc(
  drc_results,
  compare_by = "compound"   # title shows compound name
  # compare_by = "construct" -- title shows construct/kinase name instead
)
```

Generates one PNG per unique `Construct:Compound` pair, overlaying curves from all plates that contain it. Useful for assessing inter-plate reproducibility. Each plate is shown in a distinct colour; the legend is labelled with plate names.

**Key parameters:**

```r
compare_plates_drc(
  drc_results,
  compare_by       = "compound",
  output_dir       = "plate_comparison_plots",
  color_palette    = "set1",
  y_limits         = c(0, 100),   # or NULL to auto-scale to the data
  y_axis_title     = "Normalized BRET ratio [%]",
  show_error_bars  = TRUE,
  min_plates       = 2,           # skip entities found in fewer than 2 plates
  selected_entities = c("CpdA", "CpdB"),  # restrict to specific compounds
  plot_width       = 10,
  plot_height      = 8,
  plot_dpi         = 600
)
```

**Additional appearance controls** (mirrors `plot_multiple_compounds()`):

```r
compare_plates_drc(
  drc_results,
  compare_by             = "compound",
  # Axis range
  x_limits               = c(-9, -5),   # log10 molar; or NULL to auto-scale
  x_limits_scale         = "log10",      # "log10" (default), "molar", "uM", "nM"
  x_axis_title           = NULL,         # NULL uses default Log10 Concentration [M]
  # Curve appearance
  curve_linewidth        = 1,
  curve_alpha            = 0.7,
  # IC50 reference lines
  show_ic50_lines        = TRUE,
  # Title and axis styling
  plot_title_size        = 16,
  axis_line_color        = "black",
  show_border            = FALSE,
  # Background
  transparent_background = FALSE
)
```

---

## Accessing Results Directly

The output of `batch_drc_analysis()` is a nested list. Common access patterns:

```r
# Full DRC result object for a specific plate
plate1_drc <- drc_results$drc_results$plate_01$drc_result

# Final summary table (IC50, R2, Hill slope, etc.)
summary <- plate1_drc$final_summary_table

# Detailed per-compound results
compound_1 <- plate1_drc$detailed_results[[1]]

# Quality metrics
quality <- plate1_drc$curve_quality_table
```

---

## Complete Example -- NanoBRET

```r
library(dosefitr)

setwd("path/to/experiment")

# 1. Process raw plates
results <- batch_ratio_analysis(
  control_0perc   = 24,
  control_100perc = 12,
  output_dir      = "./qc_output",
  verbose         = TRUE,
  low_value_threshold = 1000
)

# 2. Remove outliers
results_clean <- rout_outliers_batch(results, Q = 0.01)
plot_outliers_batch_curves(results_clean)

# 3. Fit dose-response curves
drc_results <- batch_drc_analysis(
  batch_results = results,
  normalize     = TRUE,
  output_dir    = "./drc_results",
  verbose       = TRUE
)

# 4. Save all individual plots
batch_save_all_drc_plots(drc_results, verbose = TRUE)

# 5. Overlay selected compounds
plot_multiple_compounds(drc_results,
  plate            = "plate_02",
  compound_indices = 9:11,
  color_palette    = "colorblind",
  y_limits         = c(0, 100),
  axis_text_size   = 18,
  axis_title_size  = 20,
  legend_text_size = 20,
  save_plot        = TRUE,
  plot_width       = 12
)

# 6. Compare plates
compare_plates_drc(drc_results, compare_by = "compound")

# 7. Export Scarab table
scarab_table(
  results_list      = results,
  drc_results_list  = drc_results,
  plate_name        = "plate_02",
  experimenter_abbrev = "TL",
  nLuc_orientation  = "C",
  tracer            = "Tracer K10",
  tracer_kd_app     = -7.5,
  tracer_concentration_used = -7.5,
  decimal_separator = ",",
  eubscarab_ready   = "No"
)
```

---

## Complete Example -- Cell Viability

```r
library(dosefitr)

setwd("path/to/experiment")

# 1. Process raw plates
via_results <- batch_viability_analysis(
  control_0perc    = 13,
  control_100perc  = 12,
  selected_columns = c(2:23)
)

# 2. Remove outliers (keep cytotoxic points)
results_clean <- rout_outliers_batch(via_results, Q = 0.01, keep_cytotoxic = TRUE)
plot_outliers_batch_curves(results_clean)

# 3. Fit dose-response curves
drc_results <- batch_drc_analysis(
  batch_results = via_results,
  normalize     = TRUE,
  output_dir    = "./drc_results",
  verbose       = TRUE
)

# 4. Save all individual plots
batch_save_all_drc_plots(drc_results,
  verbose      = TRUE,
  y_axis_title = "% Cell Survival"  # optional override; auto-detects "Cell Viability (%)" by default
)

# 5. Overlay selected compounds
plot_multiple_compounds(drc_results,
  plate            = "plate_01",
  compound_indices = 9:11,
  color_palette    = "colorblind",
  y_limits         = c(0, 100),
  y_axis_title     = "% Cell Survival",  # optional override; auto-detects "Cell Viability (%)" by default
  plot_title       = "",
  axis_text_size   = 18,
  axis_title_size  = 20,
  legend_text_size = 20,
  save_plot        = TRUE,
  plot_width       = 12
)

# 6. Export Scarab viability table
scarab_viability(
  results_list        = via_results,
  drc_results_list    = drc_results,
  plate_name          = "plate_01",
  date                = "260323",
  experimenter_abbrev = "TL",
  cell_line           = "HeLa",
  treatment_time      = "72h",
  measurement_method  = "CellTiter-Glo",
  decimal_separator   = ","
)
```

---

## Function Reference

### Batch pipeline functions

| Function | Description |
|---|---|
| `batch_ratio_analysis()` | Read raw NanoBRET plates, compute BRET ratios |
| `batch_viability_analysis()` | Read raw viability plates, normalise signal |
| `merge_plate_replicates()` | Merge replicate plates into a single combined result |
| `rout_outliers_batch()` | ROUT outlier detection across all plates |
| `plot_outliers_batch_curves()` | Visualise outlier-flagged curves (batch) |
| `batch_drc_analysis()` | Fit 3PL or 4PL dose-response curves across all plates |
| `batch_save_all_drc_plots()` | Save one plot per compound across all plates |
| `plot_multiple_compounds()` | Overlay selected compounds on one plot |
| `compare_plates_drc()` | Compare the same compound across plates |
| `scarab_table()` | Generate Scarab-format export table (NanoBRET) |
| `scarab_viability()` | Generate Scarab-format export table (cell viability) |

### Single-plate / lower-level functions

These functions are called internally by the batch pipeline but can also be
used directly for interactive exploration of a single plate.

| Function | Description |
|---|---|
| `ratio_dose_response()` | Compute BRET ratios for a single plate (v1 engine) |
| `ratio_dose_response_v2()` | Compute BRET ratios for a single plate (v2 engine — supports multiple 100% control columns and fixed 0% values) |
| `process_viability_data()` | Normalise viability signal for a single plate |
| `fit_drc_3pl()` | Fit a single plate with a 3-parameter logistic model |
| `fit_drc_4pl()` | Fit a single plate with a 4-parameter logistic model |
| `plot_dose_response()` | Plot a single dose-response curve |
| `rout_outliers()` | ROUT outlier detection for a single plate |
| `plot_outliers_curves()` | Visualise outlier-flagged curves for a single plate |
