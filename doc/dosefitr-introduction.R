## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = FALSE
)

## ----install------------------------------------------------------------------
# # Install from GitHub
# remotes::install_github("AKKDataAnalysis/dosefitr")

## ----load---------------------------------------------------------------------
# library(dosefitr)

## ----nanobret-step1-----------------------------------------------------------
# results <- batch_ratio_analysis(
#   control_0perc       = 24,    # plate column for 0% control (e.g. DMSO)
#   control_100perc     = 12,    # plate column for 100% control (e.g. staurosporine)
#   output_dir          = "./qc_output",
#   low_value_threshold = 1000,  # donor-channel values below this → NA
#   verbose             = TRUE
# )

## ----nanobret-step1-v2--------------------------------------------------------
# results <- batch_ratio_analysis(
#   control_0perc       = 16,         # fixed background value
#   control_100perc     = c(12, 24),  # two columns averaged into one
#   low_value_threshold = 3000,
#   function_version    = "v2",
#   verbose             = TRUE
# )

## ----nanobret-step2-----------------------------------------------------------
# results_clean <- rout_outliers_batch(results, Q = 0.01)
# 
# # Inspect flagged points
# head(results_clean$outlier_summary)
# 
# # Visualise outlier-flagged curves
# plot_outliers_batch_curves(results_clean)

## ----nanobret-step3-----------------------------------------------------------
# drc_results <- batch_drc_analysis(
#   batch_results = results_clean,  # use results if you skipped Step 2
#   normalize     = TRUE,
#   assay_type    = "nanobret",
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )

## ----nanobret-step3-nd--------------------------------------------------------
# drc_results <- batch_drc_analysis(
#   batch_results    = results_clean,
#   normalize        = TRUE,
#   nd_if_activation = TRUE   # activation curves → N/D
# )

## ----nanobret-step3-access----------------------------------------------------
# # Final summary table for plate 1 (IC50, pIC50, R², Hill slope, ...)
# summary_table <- drc_results$drc_results$plate_01$drc_result$final_summary_table
# 
# # Quality metrics table
# quality_table <- drc_results$drc_results$plate_01$drc_result$curve_quality_table

## ----nanobret-step4-all-------------------------------------------------------
# batch_save_all_drc_plots(drc_results, verbose = TRUE)

## ----nanobret-step4-overlay---------------------------------------------------
# # By index
# plot_multiple_compounds(drc_results,
#   plate            = "plate_01",
#   compound_indices = 1:6,
#   color_palette    = "colorblind",
#   y_limits         = c(0, 100),   # or NULL to auto-scale to the data
#   y_axis_title     = "Normalized BRET ratio [%]",
#   show_error_bars  = TRUE,
#   save_plot        = TRUE
# )
# 
# # By name or partial match
# plot_multiple_compounds(drc_results, target_compound = "LRRK2")
# 
# # Additional appearance controls
# plot_multiple_compounds(drc_results,
#   compound_indices       = 1:6,
#   x_limits               = c(-9, -5),   # log10 molar; NULL auto-scales
#   x_limits_scale         = "log10",      # "log10", "molar", "uM", or "nM"
#   x_axis_title           = NULL,         # NULL uses default Log10 Concentration [M]
#   curve_linewidth        = 1,
#   curve_alpha            = 0.7,
#   show_ic50_lines        = TRUE,         # dashed vertical line at each IC50
#   plot_title_size        = 16,
#   axis_line_color        = "black",
#   show_border            = FALSE,
#   transparent_background = FALSE
# )

## ----nanobret-step4-compare---------------------------------------------------
# compare_plates_drc(
#   drc_results,
#   compare_by    = "compound",
#   color_palette = "set1",
#   y_limits      = c(0, 100),   # or NULL to auto-scale to the data
#   y_axis_title  = "Normalized BRET ratio [%]",
#   min_plates    = 2,
#   # Additional appearance controls (same as plot_multiple_compounds)
#   x_limits               = NULL,
#   show_ic50_lines        = FALSE,
#   axis_line_color        = "black",
#   show_border            = FALSE,
#   transparent_background = FALSE
# )

## ----nanobret-step5-reshape---------------------------------------------------
# final_table <- drc_results$drc_results$plate_01$drc_result$final_summary_table
# 
# excel_table <- reshape_dr_table(
#   final_table,
#   decimal_comma = FALSE,
#   output_file   = "results_plate01.xlsx"
# )

## ----nanobret-step5-multi-----------------------------------------------------
# original_table <- results_clean$plate_01$result$modified_ratio_table_original
# 
# save_multiple_sheets(
#   "final_data.xlsx",
#   original_table,
#   excel_table,
#   decimal_comma  = TRUE,
#   decimal_places = 3,
#   round_sheets   = 1:2
# )

## ----nanobret-step5-scarab----------------------------------------------------
# scarab_table(
#   results_list              = results,
#   drc_results_list          = drc_results,
#   plate_name                = "plate_01",
#   date                      = "260323",
#   experimenter_abbrev       = "TL",
#   nLuc_orientation          = "C",       # "N" or "C" terminus
#   tracer                    = "Tracer K10",
#   tracer_kd_app             = -7.5,
#   tracer_concentration_used = -7.5,
#   decimal_separator         = ",",       # use "." for English format
#   eubscarab_ready           = "No"
# )

## ----nanobret-step5-scarab-multi----------------------------------------------
# scarab_table(
#   results_list     = results,
#   drc_results_list = drc_results,
#   plate_name       = "plate_01",
#   date             = "260323",
#   experimenter_abbrev = "TL",
#   nLuc_orientation = c(LRRK2 = "C", AAK1 = "N"),
#   tracer           = c(LRRK2 = "Tracer K10", AAK1 = "Tracer 236"),
#   tracer_kd_app    = c(LRRK2 = -7.5, AAK1 = -8.0),
#   tracer_concentration_used = c(LRRK2 = -7.5, AAK1 = -8.0),
#   decimal_separator = ","
# )

## ----viability-step1----------------------------------------------------------
# via_results <- batch_viability_analysis(
#   directory       = "data/viability/",
#   control_0perc   = 1,    # column for background / vehicle control
#   control_100perc = 24,   # column for untreated cells (maximum signal)
#   output_dir      = "./qc_output",
#   verbose         = TRUE
# )

## ----viability-step1-cols-----------------------------------------------------
# via_results <- batch_viability_analysis(
#   control_0perc    = 13,
#   control_100perc  = 12,
#   selected_columns = 2:23   # exclude columns 1 and 24
# )

## ----viability-step2----------------------------------------------------------
# results_clean <- rout_outliers_batch(via_results,
#   Q              = 0.01,
#   keep_cytotoxic = TRUE
# )

## ----viability-step3----------------------------------------------------------
# drc_results <- batch_drc_analysis(
#   batch_results = results_clean,  # use via_results if you skipped Step 2
#   normalize     = TRUE,
#   assay_type    = "viability",    # applies viability-appropriate parameter limits
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )

## ----viability-step4----------------------------------------------------------
# batch_save_all_drc_plots(drc_results,
#   verbose      = TRUE,
#   y_axis_title = "Cell Viability (%)",
#   y_limits     = NULL   # NULL auto-scales each plot; use c(0, 100) for a fixed scale
# )
# 
# plot_multiple_compounds(drc_results,
#   compound_indices = 1:6,
#   color_palette    = "colorblind",
#   y_limits         = c(0, 100),   # or NULL to auto-scale to the data
#   y_axis_title     = "Cell Viability (%)",
#   save_plot        = TRUE
# )
# 
# # Additional appearance controls (same as NanoBRET workflow)
# plot_multiple_compounds(drc_results,
#   compound_indices       = 1:6,
#   x_limits               = c(-9, -5),   # log10 molar; NULL auto-scales
#   x_limits_scale         = "log10",      # "log10", "molar", "uM", or "nM"
#   x_axis_title           = NULL,
#   curve_linewidth        = 1,
#   curve_alpha            = 0.7,
#   show_ic50_lines        = TRUE,
#   plot_title_size        = 16,
#   axis_line_color        = "black",
#   show_border            = FALSE,
#   transparent_background = FALSE
# )

## ----merge--------------------------------------------------------------------
# # Merge all plates
# merged <- merge_plate_replicates(results)
# 
# # Merge a subset, keeping others separate
# merged <- merge_plate_replicates(results,
#   plates      = c("plate_01", "plate_02"),
#   merged_name = "rep1_rep2"
# )
# 
# # Save the Excel report to a custom directory (default: drc_quality/ in working dir)
# merged <- merge_plate_replicates(results, output_dir = "./my_reports")
# 
# # Feed into the DRC step as usual
# drc_results <- batch_drc_analysis(
#   batch_results = merged,
#   normalize     = TRUE,
#   output_dir    = "./drc_results"
# )

## ----full-nanobret------------------------------------------------------------
# library(dosefitr)
# 
# setwd("path/to/experiment")
# 
# # 1. Process raw plates
# results <- batch_ratio_analysis(
#   control_0perc       = 24,
#   control_100perc     = 12,
#   output_dir          = "./qc_output",
#   low_value_threshold = 1000,
#   verbose             = TRUE
# )
# 
# # 2. Remove outliers
# results_clean <- rout_outliers_batch(results, Q = 0.01)
# plot_outliers_batch_curves(results_clean)
# 
# # 3. Fit dose-response curves
# drc_results <- batch_drc_analysis(
#   batch_results = results_clean,
#   normalize     = TRUE,
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )
# 
# # 4. Save all individual plots
# batch_save_all_drc_plots(drc_results, verbose = TRUE)
# 
# # 5. Overlay selected compounds
# plot_multiple_compounds(drc_results,
#   plate            = "plate_01",
#   compound_indices = 1:6,
#   color_palette    = "colorblind",
#   y_limits         = c(0, 100),   # or NULL to auto-scale
#   y_axis_title     = "Normalized BRET ratio [%]",
#   save_plot        = TRUE
# )
# 
# # 6. Compare the same compound across plates
# compare_plates_drc(drc_results, compare_by = "compound")
# 
# # 7. Export Scarab table
# scarab_table(
#   results_list              = results,
#   drc_results_list          = drc_results,
#   plate_name                = "plate_01",
#   date                      = "260323",
#   experimenter_abbrev       = "TL",
#   nLuc_orientation          = "C",
#   tracer                    = "Tracer K10",
#   tracer_kd_app             = -7.5,
#   tracer_concentration_used = -7.5,
#   decimal_separator         = ",",
#   eubscarab_ready           = "No"
# )

## ----full-viability-----------------------------------------------------------
# library(dosefitr)
# 
# setwd("path/to/experiment")
# 
# # 1. Process raw plates
# via_results <- batch_viability_analysis(
#   control_0perc    = 13,
#   control_100perc  = 12,
#   selected_columns = 2:23,
#   verbose          = TRUE
# )
# 
# # 2. Remove outliers (keep cytotoxic points)
# results_clean <- rout_outliers_batch(via_results,
#   Q              = 0.01,
#   keep_cytotoxic = TRUE
# )
# plot_outliers_batch_curves(results_clean)
# 
# # 3. Fit dose-response curves
# drc_results <- batch_drc_analysis(
#   batch_results = results_clean,
#   normalize     = TRUE,
#   assay_type    = "viability",
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )
# 
# # 4. Save all individual plots
# batch_save_all_drc_plots(drc_results,
#   verbose      = TRUE,
#   y_axis_title = "Cell Viability (%)"
# )
# 
# # 5. Overlay selected compounds
# plot_multiple_compounds(drc_results,
#   compound_indices = 1:6,
#   color_palette    = "colorblind",
#   y_limits         = c(0, 100),   # or NULL to auto-scale
#   y_axis_title     = "Cell Viability (%)",
#   save_plot        = TRUE
# )
# 
# # 6. Export summary table
# final_table <- drc_results$drc_results$plate_01$drc_result$final_summary_table
# excel_table <- reshape_dr_table(final_table,
#   decimal_comma = FALSE,
#   output_file   = "results.xlsx"
# )
# 
# # 7. Save combined workbook
# original_table <- results_clean$plate_01$result$modified_ratio_table_original
# save_multiple_sheets("final_data.xlsx",
#   original_table, excel_table,
#   decimal_comma  = TRUE,
#   decimal_places = 3,
#   round_sheets   = 1:2
# )

