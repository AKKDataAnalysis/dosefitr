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

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

