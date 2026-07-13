# dosefitr

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/mit)
<!-- badges: end -->

End-to-end dose-response analysis for **NanoBRET kinase binding assays**
and **cell viability experiments**. `dosefitr` reads raw BMG PHERAstar
plate-reader exports, computes BRET ratios (or normalises luminescence),
detects outliers with the ROUT method, fits 3PL / 4PL logistic models,
reports IC50 / pIC50 with confidence intervals, and exports
publication-ready plots, multi-sheet Excel workbooks, and
Scarab-formatted submission tables.

## Hero example

    library(dosefitr)

    # Stage the two bundled NanoBRET plates into a working directory
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

    # 2. Fit 4PL curves and 3. save panel plots
    drc_res <- batch_drc_analysis(ratio_res, model = "4pl",
                                  generate_reports = FALSE,
                                  output_dir = file.path(tempdir(), "dosefitr_drc"))
    batch_save_all_drc_plots(drc_res,
                             output_dir = file.path(tempdir(), "dosefitr_plots"))

## Installation

`dosefitr` is planned for submission to Bioconductor once the
accompanying manuscript is published. Once accepted, install the release
version with:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("dosefitr")

The in-development version is on GitHub:

    # install.packages("remotes")
    remotes::install_github("AKKDataAnalysis/dosefitr")

System requirements: R (&gt;= 4.4.0). All package dependencies
(including `OptimModel`, which provides the ROUT outlier test) are on
CRAN and `remotes::install_github()` pulls them automatically.

## What it does

`dosefitr` turns raw 384-well plate-reader exports into a full
dose-response analysis package in one R session, without stitching
together plate-parsing, curve-fitting, and export tools:

- **NanoBRET pipeline** (dual-channel donor/acceptor):
  `batch_ratio_analysis` reads BMG PHERAstar exports, extracts donor and
  acceptor matrices, computes BRET ratios per well, and normalises each
  row by its DMSO / saturated controls.
- **Viability pipeline** (single-channel luminescence,
  e.g. CellTiter-Glo): `batch_viability_analysis` reads single-channel
  exports, normalises to DMSO, and clamps the 0 % floor at the
  compound-saturated wells. Two info-table layouts are supported via the
  `version = "v1"` / `"v2"` switch.
- **Shared downstream stack**: `rout_outliers_batch` flags outlier wells
  with the ROUT method (Motulsky & Brown 2006); `batch_drc_analysis`
  fits 3PL or 4PL models with assay-aware plausibility limits on Bottom,
  Top, and Hill slope; `batch_save_all_drc_plots` writes one plot per
  compound plus a patchwork-assembled per-plate panel;
  `plot_multiple_compounds` overlays compounds in one figure;
  `compare_plates_drc` compares the same compound across plates.
- **Automatic QC**: Z-prime, assay window, signal-to-background, and
  per-plate CV%, plus a `Curve_Quality` label per compound
  (`Good curve`, `Wide logIC50 CI range`, `Top too low`, …) exposed on
  the per-plate summary table.
- **Export**: `scarab_table` (NanoBRET) and `scarab_viability` produce
  Scarab-format submission tables for SGC-style kinase profiling
  submission; `save_multiple_sheets` writes any number of data frames to
  one Excel file.

## What it doesn’t do

`dosefitr` does not read raw fluorescence images, does not perform
competition-binding thermodynamics beyond IC50 / pIC50, and is not
designed for combination-drug (Bliss / Loewe / HSA) analysis. For
large-scale pharmacogenomic meta-analyses across cell-line panels and
public screens, see `PharmacoGx` on Bioconductor.

## Pipeline overview

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
          |--> batch_save_all_drc_plots() <- one plot per compound + panel per plate
          |--> plot_multiple_compounds()  <- overlays selected compounds
          |--> compare_plates_drc()       <- same compound across plates
          |--> scarab_table()             <- Scarab-format export (NanoBRET)
          `--> scarab_viability()         <- Scarab-format export (viability)

## Function map

<table>
<colgroup>
<col style="width: 17%" />
<col style="width: 82%" />
</colgroup>
<thead>
<tr>
<th>Group</th>
<th>Functions</th>
</tr>
</thead>
<tbody>
<tr>
<td>Reading / normalising</td>
<td><code>batch_ratio_analysis</code>
<code>batch_viability_analysis</code>
<code>process_viability_data</code>
<code>process_viability_data_v2</code> <code>ratio_dose_response</code>
<code>ratio_dose_response_v2</code>
<code>merge_plate_replicates</code></td>
</tr>
<tr>
<td>Import (pre-computed)</td>
<td><code>batch_read_tables</code></td>
</tr>
<tr>
<td>Quality control</td>
<td><code>rout_outliers</code> <code>rout_outliers_batch</code></td>
</tr>
<tr>
<td>Fitting</td>
<td><code>batch_drc_analysis</code> <code>fit_drc_4pl</code>
<code>fit_drc_3pl</code> <code>reshape_dr_table</code></td>
</tr>
<tr>
<td>Plotting</td>
<td><code>batch_save_all_drc_plots</code>
<code>plot_multiple_compounds</code> <code>plot_dose_response</code>
<code>plot_all_dose_responses</code> <code>compare_plates_drc</code>
<code>plot_outliers_curves</code>
<code>plot_outliers_batch_curves</code> <code>plot_drc_batch</code></td>
</tr>
<tr>
<td>Export</td>
<td><code>scarab_table</code> <code>scarab_viability</code>
<code>save_multiple_sheets</code></td>
</tr>
</tbody>
</table>

## Documentation

Four vignettes ship inside the package. After installing, browse the
full list in your R session with:

    browseVignettes("dosefitr")

Or open each one directly:

- **`vignette("dosefitr")`** — introduction and package tour:
  motivation, comparison with related dose-response packages, and
  pointers into the pipeline vignettes.
- **`vignette("dosefitr-nanobret")`** — end-to-end NanoBRET walkthrough
  on the bundled example plates (BRET ratios, ROUT, 3PL / 4PL fits,
  plotting, Scarab export).
- **`vignette("dosefitr-viability")`** — end-to-end cell-viability
  walkthrough on the bundled example plates (single-channel
  normalisation, 3PL / 4PL fits, Scarab export).
- **`vignette("dosefitr-protocol")`** — the day-to-day analysis script,
  section by section, mirroring how you would run it against your own
  plates (ratios, ROUT, replicate merging, DRC fits, plotting, Scarab
  metadata).

Each exported function also has its own help page; open the grouped
index with `help(package = "dosefitr")`.

## Citation

Manuscript in preparation. Run `citation("dosefitr")` for the up-to-date
record. The ROUT outlier-detection algorithm used by `rout_outliers()`
was developed by Motulsky and Brown (2006); please cite the original
paper:

> Motulsky, H. J., & Brown, R. E. (2006). Detecting outliers when
> fitting data with nonlinear regression – a new method based on robust
> nonlinear regression and the false discovery rate. *BMC
> Bioinformatics*, 7, 123.
> [doi:10.1186/1471-2105-7-123](https://doi.org/10.1186/1471-2105-7-123)

## Getting help

- Function-level docs: `help(package = "dosefitr")`
- Bug reports and feature requests:
  <https://github.com/AKKDataAnalysis/dosefitr/issues>

## License

MIT (c) 2025 Thiago Loreto Matos. See `LICENSE` for details.
