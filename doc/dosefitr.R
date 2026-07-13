## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse  = TRUE,
    comment   = "#>",
    fig.align = "center"
)

## ----install-bioc, eval = FALSE-----------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("dosefitr")

## ----install-github, eval = FALSE---------------------------------------------
# # install.packages("remotes")
# remotes::install_github("AKKDataAnalysis/dosefitr")

## ----setup--------------------------------------------------------------------
library(dosefitr)

## ----stage--------------------------------------------------------------------
extdata_dir <- system.file("extdata", package = "dosefitr")
work_dir    <- file.path(tempdir(), "dosefitr_intro")
dir.create(work_dir, showWarnings = FALSE)
invisible(file.copy(
    from      = list.files(extdata_dir, pattern = "^nanobret_",
                           full.names = TRUE),
    to        = work_dir,
    overwrite = TRUE
))
list.files(work_dir)

## ----hero---------------------------------------------------------------------
# 1. BRET ratios
ratio_res <- batch_ratio_analysis(
    directory        = work_dir,
    info_file        = "nanobret_info.xlsx",
    data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc    = "1",
    control_100perc  = "24",
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = file.path(tempdir(), "dosefitr_intro_ratio"),
    verbose          = FALSE
)

# 2. 4PL dose-response fits
drc_res <- batch_drc_analysis(
    batch_results    = ratio_res,
    model            = "4pl",
    normalize        = FALSE,
    output_dir       = file.path(tempdir(), "dosefitr_intro_drc"),
    generate_reports = FALSE,
    verbose          = FALSE
)

# 3. Per-plate summary table
sum_01 <- drc_res$drc_results$plate_01$drc_result$summary_table
sum_01[, c("Compound", "LogIC50", "HillSlope",
           "R_squared", "Curve_Quality")]

## ----quality------------------------------------------------------------------
table(sum_01$Curve_Quality)

## ----import-precomputed, eval = FALSE-----------------------------------------
# # Minimal: directory defaults to getwd(), assay_source defaults to "nanobret",
# # and the first column is taken as the log-concentration column regardless
# # of its exact name.
# imported <- batch_read_tables()
# 
# # Or, spelling every argument out (useful for viability, for a custom
# # directory, or when the log-concentration column is not the first one).
# imported <- batch_read_tables(
#     directory         = "path/to/precomputed_tables",
#     assay_source      = "viability",            # or "nanobret" (default)
#     log_conc_col_name = "log(inhibitor).[M]"    # look up by name; if it is
#                                                 # not already column 1 it
#                                                 # will be moved there
# )
# 
# drc_res <- batch_drc_analysis(
#     batch_results = imported,
#     model         = "3pl",                      # or "4pl"
#     normalize     = TRUE                        # set FALSE if your table is
#                                                 # already control-normalised
# )

## ----import-file-map, eval = FALSE--------------------------------------------
# # Explicit file list. `file_pattern` is ignored when `file_map` is set.
# # Plates are named plate_01, plate_02, ... in the order given.
# imported <- batch_read_tables(
#     directory = "path/to/precomputed_tables",
#     file_map  = c("run_A.xlsx", "run_B.xlsx", "run_C.xlsx")
# )

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

