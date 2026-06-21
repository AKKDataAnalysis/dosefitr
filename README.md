# dosefitr

> Dose-Response Curve Fitting for NanoBRET and Cell Viability Assays

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/dosefitr)](https://CRAN.R-project.org/package=dosefitr)
[![R-CMD-check](https://github.com/AKKDataAnalysis/dosefitr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AKKDataAnalysis/dosefitr/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/mit/)
<!-- badges: end -->

End-to-end dose-response analysis for **NanoBRET kinase binding assays**
and **cell viability experiments**. **dosefitr** reads raw BMG
PHERAstar plate-reader exports, computes BRET ratios, detects
outliers (ROUT method), fits 3PL or 4PL logistic models, reports
IC50 / pIC50 with confidence intervals, and exports
publication-ready plots, multi-sheet Excel workbooks, and
Scarab-formatted submission tables.

## Installation

You can install the stable release from CRAN with

``` r
install.packages("dosefitr")
```

To install the development version from GitHub, use

``` r
# install.packages("remotes")
remotes::install_github("AKKDataAnalysis/dosefitr")
```

`dosefitr` requires R >= 4.1.0 and one non-CRAN dependency,
`OptimModel` (used internally by the ROUT outlier test). On most
systems `install.packages("dosefitr")` will pull `OptimModel` from
CRAN automatically.

## Quick example

The package ships with two small example plates and an info table in
`inst/extdata/` so you can copy the example below straight into your
own R session:

``` r
library(dosefitr)

# Find and stage the bundled plates
extdata_dir <- system.file("extdata", package = "dosefitr")
work_dir    <- file.path(tempdir(), "dosefitr_demo")
dir.create(work_dir, showWarnings = FALSE)
file.copy(
  list.files(extdata_dir, pattern = "^nanobret_", full.names = TRUE),
  work_dir
)

# 1. Compute BRET ratios across both plates
ratio_res <- batch_ratio_analysis(
  directory        = work_dir,
  info_file        = "nanobret_info.xlsx",
  data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
  control_0perc    = "1",
  control_100perc  = "24",
  selected_columns = 1:24,
  generate_reports = FALSE,
  output_dir       = file.path(tempdir(), "dosefitr_ratio")
)

# 2. Fit 4PL dose-response curves
drc_res <- batch_drc_analysis(
  batch_results    = ratio_res,
  model            = "4pl",
  generate_reports = FALSE,
  output_dir       = file.path(tempdir(), "dosefitr_drc")
)

# 3. Inspect the per-compound summary
drc_res$drc_results$plate_01$drc_result$summary_table[, c(
  "Compound", "LogIC50", "HillSlope", "R_squared", "Curve_Quality"
)]

# 4. Plot one compound
plot_drc_batch(
  batch_drc_results  = drc_res,
  construct_compound = "KinaseA:Cpd1"
)
```

For the full pipeline (outlier detection, multi-compound overlay,
cross-plate comparison, export), see the pipeline vignettes.

## Pipeline overview

```
Raw Excel files
      |
      v
batch_ratio_analysis()        <- NanoBRET:   reads plates, computes BRET ratios
batch_viability_analysis()    <- Viability:  reads plates, normalises signal
      |
      v
rout_outliers_batch()         <- detects and removes outliers (optional)
      |
      v
batch_drc_analysis()          <- fits dose-response curves (3PL or 4PL)
      |
      |--> plot_drc_batch()             <- one or several compounds, faceted
      |--> plot_multiple_compounds()    <- overlays selected compounds
      |--> compare_plates_drc()         <- compares same compound across plates
      |--> scarab_table()               <- Scarab-format export (NanoBRET)
      `--> scarab_viability()           <- Scarab-format export (viability)
```

## Function reference

Each function has its own help page; below is the grouped index.

### Reading and normalising raw data

| Function | Purpose |
|----------|---------|
| `batch_ratio_analysis()`     | Read NanoBRET plates and compute BRET ratios |
| `batch_viability_analysis()` | Read viability plates and normalise signal (`version = "v1"` or `"v2"`) |
| `process_viability_data()`   | Single-plate v1 viability processor |
| `process_viability_data_v2()` | Single-plate v2 viability processor |
| `ratio_dose_response()`      | Single-plate v1 NanoBRET ratio processor |
| `ratio_dose_response_v2()`   | Single-plate v2 NanoBRET ratio processor |
| `merge_plate_replicates()`   | Merge technical replicates across plates |

### Quality control

| Function | Purpose |
|----------|---------|
| `rout_outliers()`       | ROUT outlier test on a single fitted curve |
| `rout_outliers_batch()` | ROUT outlier test on all curves in a batch result |

### Dose-response fitting

| Function | Purpose |
|----------|---------|
| `batch_drc_analysis()` | Fit 3PL or 4PL models on all compounds in a batch |
| `fit_drc_4pl()`        | Single-compound 4PL fit |
| `fit_drc_3pl()`        | Single-compound 3PL fit |
| `reshape_dr_table()`   | Reshape DRC results between wide and long forms |

### Plotting

| Function | Purpose |
|----------|---------|
| `plot_drc_batch()`            | Fitted curves for one or several construct-compound pairs |
| `plot_multiple_compounds()`   | Overlay several compounds in a single plot |
| `plot_dose_response()`        | Single-compound plot |
| `plot_all_dose_responses()`   | One plot per compound (returns a list) |
| `compare_plates_drc()`        | Compare the same compound across plates |
| `plot_outliers_curves()`      | Highlight ROUT-flagged wells on one plate |
| `plot_outliers_batch_curves()`| Highlight ROUT-flagged wells across plates |
| `batch_save_all_drc_plots()`  | Save one plot per compound to a directory |

### Export

| Function | Purpose |
|----------|---------|
| `scarab_table()`        | Scarab-format submission table (NanoBRET) |
| `scarab_viability()`    | Scarab-format submission table (viability) |
| `save_multiple_sheets()`| Write any number of data frames to one Excel file |

## Vignettes

Two pipeline vignettes walk through the full workflow on the bundled
example data; both are fully evaluated at build time.

-   `vignette("dosefitr-nanobret", package = "dosefitr")` --
    BRET-ratio computation, outlier detection, 4PL fitting, plotting,
    and Scarab export for NanoBRET kinase binding assays.
-   `vignette("dosefitr-viability", package = "dosefitr")` --
    Single-channel viability normalisation, 4PL fitting, plotting,
    and Scarab export for cell viability experiments.

## Citation

When using **dosefitr** in published work, please cite the package:

``` r
citation("dosefitr")
```

The ROUT outlier-detection algorithm used by `rout_outliers()` was
developed by Motulsky & Brown (2006); please cite the original paper:

> Motulsky, H. J., & Brown, R. E. (2006). Detecting outliers when
> fitting data with nonlinear regression -- a new method based on
> robust nonlinear regression and the false discovery rate. *BMC
> Bioinformatics*, 7, 123.
> [doi:10.1186/1471-2105-7-123](https://doi.org/10.1186/1471-2105-7-123)

## Getting help

-   Browse the function index in the package help: `?dosefitr`
-   Report bugs or request features at
    <https://github.com/AKKDataAnalysis/dosefitr/issues>

## License

MIT (c) 2025 Thiago Loreto Matos. See `LICENSE` for details.
