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
#   low_value_threshold = 1000,  # donor-channel values below this -> NA
#   verbose             = TRUE
# )

## ----nanobret-step1-colname---------------------------------------------------
# results <- batch_ratio_analysis(
#   control_0perc   = "DMSO",
#   control_100perc = "Staurosporine",
#   verbose         = TRUE
# )

## ----nanobret-step1-v2--------------------------------------------------------
# results <- batch_ratio_analysis(
#   control_0perc       = 16,         # fixed background value
#   control_100perc     = c(12, 24),  # two columns averaged into one
#   low_value_threshold = 3000,
#   function_version    = "v2",
#   verbose             = TRUE
# )

## ----nanobret-step1-sep, eval=FALSE-------------------------------------------
# # Show "EPHA1/KK135" instead of "EPHA1:KK135" in legends and titles
# plot_multiple_compounds(drc_results, label_sep = "/")

## ----legend-width, eval=FALSE-------------------------------------------------
# plot_multiple_compounds(drc_results, legend_width = 4)

## ----legend-width-auto, eval=FALSE--------------------------------------------
# p1 <- plot_multiple_compounds(drc_results, compound_indices = 1:5,
#                               legend_width = "auto")
# max_w <- attr(p1, "metadata")$legend_width_cm
# 
# # Then pad all plots to the maximum width
# plot_multiple_compounds(drc_results, compound_indices = 1:5,
#                         legend_width = max_w)

## ----nanobret-step2-----------------------------------------------------------
# results_clean <- rout_outliers_batch(
#   batch_results     = results,
#   Q                 = 0.01,      # 1% FDR -- recommended starting point
#   n_param           = 4L,        # 4-parameter Hill model (free Hill slope)
#   direction         = "inhibition",
#   ntry_retry        = 3L,        # random restarts on non-convergence
#   min_dynamic_range = 20,        # warn if dynamic range < 20%
#   keep_cytotoxic    = FALSE,     # set TRUE for viability assays
#   seed              = 42L,       # ensures reproducible outlier calls
#   verbose           = TRUE
# )
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
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )

## ----nanobret-step3-nd--------------------------------------------------------
# drc_results <- batch_drc_analysis(
#   batch_results    = results_clean,
#   normalize        = TRUE,
#   nd_if_activation = TRUE   # activation curves -> N/D
# )

## ----nanobret-step3-access----------------------------------------------------
# # Final summary table for plate 1 (IC50, pIC50, R2, Hill slope, ...)
# summary_table <- drc_results$drc_results$plate_01$drc_result$final_summary_table
# 
# # Quality metrics table
# quality_table <- drc_results$drc_results$plate_01$drc_result$curve_quality_table

## ----nanobret-step4-all-------------------------------------------------------
# batch_save_all_drc_plots(
#   batch_drc_results    = drc_results,
#   output_dir           = "DRC_Plots",
#   verbose              = TRUE,
#   y_axis_title         = NULL,          # auto-detects assay type
#   y_limits             = NULL,          # NULL auto-scales; use c(0, 100) for fixed scale
#   organize_by          = "plate",       # "plate" or "compound"
#   format               = "png",
#   width                = 10,
#   height               = 10,
#   dpi                  = 600,
#   compounds_to_plot    = NULL,          # character vector to restrict compounds
#   plates_to_plot       = NULL,          # character vector to restrict plates
#   save_panel           = TRUE,          # also save a combined panel per plate
#   panel_ncol           = 4L,
#   panel_width_per_col  = 6,
#   panel_height_per_row = 6,
#   panel_spacing        = 1,
#   subplot_title        = "auto"         # "auto", "compound", "construct", or fixed string
# )

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
# plot_multiple_compounds(drc_results, target_compound = "KinaseA")
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
#   transparent_background = FALSE,
#   # Consistent legend widths across comparison plots
#   legend_width           = NULL          # NULL = no padding; "auto" = two-pass measure + pad
# )

## ----compare-legend-width, eval=FALSE-----------------------------------------
# compare_plates_drc(drc_results, legend_width = "auto")

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
#   nLuc_orientation = c(KinaseA = "C", KinaseB = "N"),
#   tracer           = c(KinaseA = "Tracer K10", KinaseB = "Tracer 236"),
#   tracer_kd_app    = c(KinaseA = -7.5, KinaseB = -8.0),
#   tracer_concentration_used = c(KinaseA = -7.5, KinaseB = -8.0),
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
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )

## ----viability-step4----------------------------------------------------------
# batch_save_all_drc_plots(drc_results,
#   verbose      = TRUE,
#   y_axis_title = "% Cell Survival",  # optional override; auto-detects "Cell Viability (%)" by default
#   y_limits     = NULL   # NULL auto-scales each plot; use c(0, 100) for a fixed scale
# )
# 
# plot_multiple_compounds(drc_results,
#   compound_indices = 1:6,
#   color_palette    = "colorblind",
#   y_limits         = c(0, 100),   # or NULL to auto-scale to the data
#   y_axis_title     = "% Cell Survival",  # optional override; auto-detects "Cell Viability (%)" by default
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

## ----viability-step5----------------------------------------------------------
# scarab_viability(
#   results_list        = via_results,
#   drc_results_list    = drc_results,
#   plate_name          = "plate_01",
#   date                = "260323",
#   experimenter_abbrev = "TL",
#   cell_line           = "HeLa",
#   cell_type           = "cervical adenocarcinoma",
#   treatment_time      = "72h",
#   measurement_method  = "CellTiter-Glo",
#   assay_volume        = "40",
#   plate_format        = "384-well",
#   decimal_separator   = ","   # use "." for English format
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
#   output_dir    = "./drc_results",
#   verbose       = TRUE
# )
# 
# # 4. Save all individual plots
# batch_save_all_drc_plots(drc_results,
#   verbose      = TRUE,
#   y_axis_title = "% Cell Survival"  # optional override; auto-detects "Cell Viability (%)" by default
# )
# 
# # 5. Overlay selected compounds
# plot_multiple_compounds(drc_results,
#   compound_indices = 1:6,
#   color_palette    = "colorblind",
#   y_limits         = c(0, 100),   # or NULL to auto-scale
#   y_axis_title     = "% Cell Survival",  # optional override; auto-detects "Cell Viability (%)" by default
#   save_plot        = TRUE
# )
# 
# # 6. Export Scarab viability table
# scarab_viability(
#   results_list        = via_results,
#   drc_results_list    = drc_results,
#   plate_name          = "plate_01",
#   date                = "260323",
#   experimenter_abbrev = "TL",
#   cell_line           = "HeLa",
#   treatment_time      = "72h",
#   measurement_method  = "CellTiter-Glo",
#   decimal_separator   = ","
# )

## ----single-plate-ratio-------------------------------------------------------
# # Load a single raw data frame
# raw_data  <- readxl::read_excel("NanoBRET_plate_1.xlsx", sheet = 1, col_names = FALSE)
# info      <- readxl::read_excel("info_tables.xlsx", sheet = "Sheet1")
# 
# # Compute BRET ratios for one plate
# plate_result <- ratio_dose_response_v2(
#   data             = raw_data,
#   control_0perc    = 16,
#   control_100perc  = c(12, 24),
#   info_table       = info,
#   split_replicates = TRUE,
#   verbose          = TRUE
# )
# 
# # Detect outliers for one plate
# plate_clean <- rout_outliers(
#   data      = plate_result$modified_ratio_table,
#   Q         = 0.01,
#   seed      = 42L
# )
# 
# # Visualise outlier-flagged curves for one plate
# plot_outliers_curves(
#   rout_output = plate_clean,
#   title       = "Plate 1"
# )

## ----single-plate-viability---------------------------------------------------
# raw_via  <- readxl::read_excel("viability_plate_1.xlsx", sheet = 1, col_names = FALSE)
# info_via <- readxl::read_excel("info_tables.xlsx", sheet = "Sheet1")
# 
# plate_via <- process_viability_data(
#   data             = raw_via,
#   control_0perc    = 1,
#   control_100perc  = 24,
#   info_table       = info_via,
#   split_replicates = TRUE,
#   verbose          = TRUE
# )

## ----single-plate-plot--------------------------------------------------------
# # Fit a single plate
# fit <- fit_drc_3pl(data = plate_clean$cleaned_table, normalize = TRUE)
# 
# # Plot compound 1
# plot_dose_response(
#   results        = fit,
#   compound_index = 1,
#   y_limits       = c(0, 120),
#   x_limits       = c(-9, -5),     # log10 molar; NULL auto-scales
#   show_ic50_line = TRUE,
#   save_plot      = "compound1_curve.png"
# )

