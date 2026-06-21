## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  eval      = TRUE,
  message   = FALSE,
  warning   = FALSE
)

## ----libs---------------------------------------------------------------------
library(dosefitr)

## ----locate-data--------------------------------------------------------------
extdata_dir <- system.file("extdata", package = "dosefitr")

list.files(extdata_dir, pattern = "^nanobret_")

## ----stage--------------------------------------------------------------------
work_dir <- file.path(tempdir(), "dosefitr_nanobret")
dir.create(work_dir, showWarnings = FALSE)

invisible(file.copy(
  from      = list.files(extdata_dir, pattern = "^nanobret_",
                         full.names = TRUE),
  to        = work_dir,
  overwrite = TRUE
))

list.files(work_dir)

## ----show-info----------------------------------------------------------------
info_path <- file.path(work_dir, "nanobret_info.xlsx")

info_plate_01 <- openxlsx::read.xlsx(info_path, sheet = "plate_01")

head(info_plate_01, 8)

## ----ratio--------------------------------------------------------------------
ratio_out <- file.path(tempdir(), "dosefitr_nanobret_ratio")
dir.create(ratio_out, showWarnings = FALSE)

ratio_res <- batch_ratio_analysis(
  directory        = work_dir,
  info_file        = "nanobret_info.xlsx",
  data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
  control_0perc    = "1",
  control_100perc  = "24",
  selected_columns = 1:24,
  generate_reports = FALSE,
  output_dir       = ratio_out,
  verbose          = FALSE
)

names(ratio_res)

## ----ratio-inspect------------------------------------------------------------
mrt <- ratio_res$plate_01$result$modified_ratio_table

dim(mrt)

head(mrt[, 1:6], 8)

## ----ratio-direction----------------------------------------------------------
data.frame(
  log_conc = mrt[, 1],
  Cpd1_rep1 = round(mrt[, "KinaseA:Cpd1"],   1),
  Cpd1_rep2 = round(mrt[, "KinaseA:Cpd1.2"], 1)
)

## ----rout---------------------------------------------------------------------
rout_res <- rout_outliers_batch(
  batch_results = ratio_res,
  Q             = 0.01,
  n_param       = 4L,
  direction     = "inhibition",
  verbose       = FALSE
)

# Summary of detection
nrow(rout_res$outlier_summary)

table(rout_res$outlier_summary$plate)

## ----rout-summary-------------------------------------------------------------
head(rout_res$outlier_summary, 6)

## ----rout-plot, fig.height = 4------------------------------------------------
plot_outliers_curves(
  rout_output   = rout_res$plate_01,
  subplot_title = "compound",
  ncol          = 4L
)

## ----drc----------------------------------------------------------------------
drc_out <- file.path(tempdir(), "dosefitr_nanobret_drc")
dir.create(drc_out, showWarnings = FALSE)

drc_res <- batch_drc_analysis(
  batch_results    = ratio_res,
  model            = "4pl",
  normalize        = FALSE,
  output_dir       = drc_out,
  generate_reports = FALSE,
  verbose          = FALSE
)

drc_res$drc_results$plate_01$drc_result$successful_fits

## ----drc-summary--------------------------------------------------------------
sum_01 <- drc_res$drc_results$plate_01$drc_result$summary_table

sum_01[, c("Compound", "LogIC50", "HillSlope",
           "R_squared", "Curve_Quality")]

## ----drc-quality--------------------------------------------------------------
table(sum_01$Curve_Quality)

## ----plot-drc, fig.height = 4-------------------------------------------------
plot_drc_batch(
  batch_drc_results  = drc_res,
  construct_compound = "KinaseA:Cpd1",
  show_error_bars    = TRUE,
  verbose            = FALSE
)

## ----plot-overlay-------------------------------------------------------------
plot_multiple_compounds(
  results          = drc_res,
  compound_indices = 1:4,
  plate            = "plate_01",
  y_limits         = NULL,
  show_error_bars  = TRUE,
  verbose          = FALSE
)

## ----compare-plates, fig.show = "hide"----------------------------------------
compare_out <- file.path(tempdir(), "dosefitr_nanobret_compare")
dir.create(compare_out, showWarnings = FALSE)

cmp_res <- compare_plates_drc(
  batch_drc_result = drc_res,
  compare_by       = "compound",
  output_dir       = compare_out,
  min_plates       = 2,
  verbose          = FALSE
)

list.files(compare_out, pattern = "\\.png$")

## ----scarab-------------------------------------------------------------------
sc_01 <- scarab_table(
  results_list      = ratio_res,
  drc_results_list  = drc_res,
  plate_name        = "plate_01",
  date              = format(Sys.Date(), "%y%m%d"),
  experimenter_abbrev = "DEMO",
  save              = FALSE
)

dim(sc_01)

head(sc_01[, 1:5], 6)

## ----save-multi---------------------------------------------------------------
xlsx_out <- file.path(tempdir(), "dosefitr_nanobret_export.xlsx")

save_multiple_sheets(
  file_name     = xlsx_out,
  Ratios        = mrt,
  Fit_Summary   = sum_01,
  decimal_comma = FALSE
)

file.exists(xlsx_out)

