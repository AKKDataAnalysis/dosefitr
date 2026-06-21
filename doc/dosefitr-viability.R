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

list.files(extdata_dir, pattern = "^viability_")

## ----stage--------------------------------------------------------------------
work_dir <- file.path(tempdir(), "dosefitr_viability")
dir.create(work_dir, showWarnings = FALSE)

invisible(file.copy(
  from      = list.files(extdata_dir, pattern = "^viability_",
                         full.names = TRUE),
  to        = work_dir,
  overwrite = TRUE
))

list.files(work_dir)

## ----show-info----------------------------------------------------------------
info_path <- file.path(work_dir, "viability_info.xlsx")

info_plate_01 <- openxlsx::read.xlsx(info_path, sheet = "plate_01")

head(info_plate_01, 8)

## ----viability----------------------------------------------------------------
via_out <- file.path(tempdir(), "dosefitr_viability_proc")
dir.create(via_out, showWarnings = FALSE)

via_res <- batch_viability_analysis(
  directory        = work_dir,
  info_file        = "viability_info.xlsx",
  data_pattern     = "viability_plate_\\d+\\.xlsx$",
  control_0perc    = 1,
  control_100perc  = 24,
  selected_columns = 1:24,
  generate_reports = FALSE,
  output_dir       = via_out,
  verbose          = FALSE,
  version          = "v1"
)

names(via_res)

## ----viability-inspect--------------------------------------------------------
mrt <- via_res$plate_01$result$modified_ratio_table

dim(mrt)

head(mrt[, 1:6], 8)

## ----viability-direction------------------------------------------------------
data.frame(
  log_conc  = mrt[, 1],
  Cpd1_rep1 = round(mrt[, "KinaseA:Cpd1"],   0),
  Cpd1_rep2 = round(mrt[, "KinaseA:Cpd1.2"], 0)
)

## ----drc----------------------------------------------------------------------
drc_out <- file.path(tempdir(), "dosefitr_viability_drc")
dir.create(drc_out, showWarnings = FALSE)

drc_res <- batch_drc_analysis(
  batch_results    = via_res,
  model            = "4pl",
  normalize        = TRUE,
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

## ----batch-save-plots---------------------------------------------------------
plots_dir <- file.path(tempdir(), "dosefitr_viability_plots")
dir.create(plots_dir, showWarnings = FALSE)

results2 <- batch_save_all_drc_plots(
  batch_drc_results = drc_res,
  output_dir        = plots_dir,
  verbose           = FALSE,
  dpi               = 120,            # lighter for the vignette
  panel_width_per_col  = 4,
  panel_height_per_row = 4,
  panel_spacing        = 1.2
)

results2$successes
basename(results2$panel_files)

## ----show-panel, echo = FALSE, out.width = "100%"-----------------------------
knitr::include_graphics(results2$panel_files[1])

## ----plot-overlay-------------------------------------------------------------
plot_multiple_compounds(
  results          = drc_res,
  compound_indices = 1:4,
  plate            = "plate_01",
  show_error_bars  = TRUE,
  verbose          = FALSE
)

## ----compare-plates, fig.show = "hide"----------------------------------------
compare_out <- file.path(tempdir(), "dosefitr_viability_compare")
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
sc_01 <- scarab_viability(
  results_list      = via_res,
  drc_results_list  = drc_res,
  plate_name        = "plate_01",
  date              = format(Sys.Date(), "%y%m%d"),
  experimenter_abbrev = "DEMO",
  save              = FALSE
)

dim(sc_01)

head(sc_01[, 1:5], 6)

## ----save-multi---------------------------------------------------------------
xlsx_out <- file.path(tempdir(), "dosefitr_viability_export.xlsx")

save_multiple_sheets(
  file_name     = xlsx_out,
  Viability     = mrt,
  Fit_Summary   = sum_01,
  decimal_comma = FALSE
)

file.exists(xlsx_out)

## ----v2-note, eval = FALSE----------------------------------------------------
# via_res_v2 <- batch_viability_analysis(
#   directory        = work_dir,
#   info_file        = "viability_info.xlsx",
#   data_pattern     = "viability_plate_\\d+\\.xlsx$",
#   control_0perc    = 1,
#   control_100perc  = 24,
#   selected_columns = 1:24,
#   generate_reports = FALSE,
#   output_dir       = via_out,
#   version          = "v2"
# )

