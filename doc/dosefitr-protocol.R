## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  eval       = TRUE,
  message    = FALSE,
  warning    = FALSE
)

## ----install-bioc, eval = FALSE-----------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("dosefitr")

## ----install-github, eval = FALSE---------------------------------------------
# # install.packages("remotes")
# remotes::install_github("AKKDataAnalysis/dosefitr")

## ----libs---------------------------------------------------------------------
library(dosefitr)

## ----stage--------------------------------------------------------------------
extdata_dir <- system.file("extdata", package = "dosefitr")
work_dir    <- file.path(tempdir(), "dosefitr_protocol")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

invisible(file.copy(
    from      = list.files(extdata_dir, pattern = "^nanobret_",
                           full.names = TRUE),
    to        = work_dir,
    overwrite = TRUE
))

list.files(work_dir)

## ----ratio--------------------------------------------------------------------
ratio_out <- file.path(tempdir(), "dosefitr_protocol_ratio")
dir.create(ratio_out, showWarnings = FALSE)

ratio_res <- batch_ratio_analysis(
    directory           = work_dir,
    info_file           = "nanobret_info.xlsx",
    data_pattern        = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc       = "1",
    control_100perc     = "24",
    selected_columns    = 1:24,
    low_value_threshold = 3000,
    generate_reports    = FALSE,
    output_dir          = ratio_out,
    verbose             = FALSE
)

names(ratio_res)

mrt <- ratio_res$plate_01$result$modified_ratio_table
dim(mrt)
head(mrt[, 1:5], 6)

## ----own-project-ratio, eval = FALSE------------------------------------------
# ratio_res <- batch_ratio_analysis(
#     directory           = "C:/path/to/your/dataset",
#     control_0perc       = 24,
#     control_100perc     = 12,
#     low_value_threshold = 3000,
#     output_dir          = "./drc_quality",
#     verbose             = TRUE
# )

## ----ratio-v2, eval = FALSE---------------------------------------------------
# ratio_res_v2 <- batch_ratio_analysis(
#     directory           = "C:/path/to/your/dataset",
#     control_0perc       = 16,
#     control_100perc     = c(12, 24),
#     low_value_threshold = 3000,
#     function_version    = "v2",
#     output_dir          = "./drc_quality",
#     verbose             = TRUE
# )

## ----rout---------------------------------------------------------------------
rout_res <- rout_outliers_batch(ratio_res, Q = 0.01, verbose = FALSE)

nrow(rout_res$outlier_summary)
head(rout_res$outlier_summary, 5)

## ----rout-plot, fig.width = 8, fig.height = 6---------------------------------
pobc_out <- file.path(tempdir(), "dosefitr_protocol_pobc")
dir.create(pobc_out, showWarnings = FALSE)

pobc_res <- plot_outliers_batch_curves(
    rout_res,
    plates     = "plate_01",
    output_dir = pobc_out,
    verbose    = FALSE
)

## ----merge--------------------------------------------------------------------
merge_out <- file.path(tempdir(), "dosefitr_protocol_merge")
dir.create(merge_out, showWarnings = FALSE)

merged_res <- merge_plate_replicates(
    ratio_res,
    output_dir       = merge_out,
    generate_reports = FALSE
)

names(merged_res)

merged_mrt <- merged_res[[1]]$result$modified_ratio_table
dim(merged_mrt)
head(colnames(merged_mrt), 8)

## ----drc----------------------------------------------------------------------
drc_out <- file.path(tempdir(), "dosefitr_protocol_drc")
dir.create(drc_out, showWarnings = FALSE)

drc_res <- batch_drc_analysis(
    batch_results    = rout_res,
    model            = "3pl",
    normalize        = TRUE,
    output_dir       = drc_out,
    generate_reports = FALSE,
    verbose          = FALSE
)

sum_01 <- drc_res$drc_results$plate_01$drc_result$summary_table
sum_01[, c("Compound", "LogIC50", "R_squared", "Curve_Quality")]

## ----save-plots---------------------------------------------------------------
plots_out <- file.path(tempdir(), "dosefitr_protocol_plots")
dir.create(plots_out, showWarnings = FALSE)

save_res <- batch_save_all_drc_plots(
    batch_drc_results = drc_res,
    output_dir        = plots_out,
    dpi               = 72,
    save_panel        = FALSE,
    verbose           = FALSE
)

length(list.files(plots_out, recursive = TRUE, pattern = "\\.png$"))

## ----overlay, fig.width = 7, fig.height = 5-----------------------------------
overlay <- plot_multiple_compounds(
    drc_res,
    compound_indices  = 1:6,
    plate             = "plate_01",
    color_palette     = "colorblind",
    plot_title        = "",
    legend_text_size  = 12,
    axis_text_size    = 11,
    axis_title_size   = 12,
    save_plot         = NULL,
    verbose           = FALSE
)

overlay

## ----compare------------------------------------------------------------------
cmp_out <- file.path(tempdir(), "dosefitr_protocol_cmp")
dir.create(cmp_out, showWarnings = FALSE)

cmp_res <- compare_plates_drc(
    drc_res,
    compare_by        = "compound",
    selected_entities = c("KinaseA:Cpd1", "KinaseB:Cpd9"),
    min_plates        = 2,
    output_dir        = cmp_out,
    plot_dpi          = 72,
    verbose           = FALSE
)

length(cmp_res)
names(cmp_res)

## ----scarab-------------------------------------------------------------------
scarab <- scarab_table(
    ratio_res, drc_res,
    plate_name                 = "plate_01",
    date                       = "240101",
    experimenter_abbrev        = "ME",
    nLuc_orientation           = "C",
    tracer_kd_app              = -7.5,
    tracer_concentration_used  = -7.5,
    tracer                     = "Tracer K10",
    decimal_separator          = ","
)

dim(scarab)
head(scarab[, 1:5], 6)

## ----viability-sidebar, eval = FALSE------------------------------------------
# via_res <- batch_viability_analysis(
#     directory        = "C:/path/to/your/dataset",
#     control_0perc    = 13,
#     control_100perc  = 12,
#     selected_columns = c(2:23),
#     output_dir       = "./drc_quality",
#     verbose          = TRUE
# )
# 
# drc_res_via <- batch_drc_analysis(
#     batch_results = via_res,
#     model         = "3pl",
#     normalize     = TRUE,
#     output_dir    = "./drc_results",
#     verbose       = TRUE
# )

## ----session-info, echo = FALSE-----------------------------------------------
sessionInfo()

