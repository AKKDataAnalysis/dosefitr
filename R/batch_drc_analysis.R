#' Batch Dose-Response Curve (DRC) Analysis for Multiple Plates
#'
#' @description
#' `batch_drc_analysis()` performs automated dose-response curve (DRC) fitting
#' across multiple plates previously processed by `batch_ratio_analysis()`.
#'
#' It extracts the `modified_ratio_table` from each plate, applies a
#' 3-parameter logistic regression via [`fit_drc_3pl()`] (default) or a
#' 4-parameter logistic regression via [`fit_drc_4pl()`] (selected via the
#' `model` argument), and generates:
#'
#' * Per-plate DRC result files (optional)
#' * A consolidated batch DRC report (`batch_drc_analysis_report.xlsx`)
#'
#' The function computes fitting success rates, R-squared metrics, curve
#' qualities, and compiles them into structured summary sheets.
#'
#'
#' @param batch_results A named list containing results from
#'   [`batch_ratio_analysis()`]. Each element must contain the field
#'   `result$modified_ratio_table`, which will be used for the DRC fit.
#' @param normalize Logical. Whether to normalize responses inside `fit_drc_3pl()` or `fit_drc_4pl()`. Values TRUE or FALSE
#' @param enforce_bottom_threshold Logical. If `TRUE`, forces the lower plateau of
#'   the fitted curve to stay above `bottom_threshold`.
#' @param bottom_threshold Numeric. Minimum acceptable bottom asymptote value.
#' @param r_sqr_threshold Minimum acceptable R-squared for accepting a curve fit.
#' @param model Character. Which dose-response model to use for fitting.
#'   `"3pl"` (default) calls [`fit_drc_3pl()`] (3-parameter logistic, Hill
#'   slope fixed at \eqn{\pm 1}).  `"4pl"` calls [`fit_drc_4pl()`]
#'   (4-parameter logistic, Hill slope freely estimated).
#' @param output_dir Directory where individual plate results and consolidated
#'   batch reports will be saved. Defaults to the working directory.
#' @param generate_reports Logical. If `TRUE` (default), generates a consolidated
#'   Excel report summarizing all DRC fits.
#' @param nd_if_activation Logical. If `FALSE` (default), only flat curves have
#'   their IC50 and pIC50 replaced with `"N/D"` (not determined) in the
#'   Pharmacology_Summary table. If `TRUE`, curves classified as going up
#'   (activation) are also treated as N/D. Curve type is determined inside the
#'   fitting function: a curve is flat when the difference between the
#'   mean of the first three and last three responses
#'   is below `max(15, 15\% of the response range)`; it is activation when the
#'   final responses are higher than the initial ones by that same threshold.
#' @param verbose Logical. If `TRUE`, prints progress details.
#' @param hook_effect Controls automatic hook-effect detection at the two highest
#'   concentration points of inhibition curves. A hook is declared when the mean
#'   response at a candidate concentration exceeds the fitted curve by more than
#'   `hook_threshold` times the residual standard error (Syx) of the first-pass fit.
#'   Accepted values:
#'   \itemize{
#'     \item `FALSE` (default) — disabled; no hook detection.
#'     \item `TRUE` — applied to every compound in the batch.
#'     \item A character vector of `"Construct:Compound"` names — applied only to
#'       those compounds (e.g. `c("BRPF1B:LW501", "BRPF1B:LW454")`).
#'   }
#'   Detected hook points are excluded from fitting but still plotted as red
#'   \eqn{\times} symbols in `plot_dose_response()`.
#' @param hook_threshold Positive numeric. Number of residual standard errors
#'   (Syx) above the fitted curve that a candidate point must exceed to be
#'   classified as a hook. Default `2` (approximately a 2-SE exceedance).
#'   Only used when `hook_effect` is not `FALSE`.
#'
#'
#' @details
#' For each plate in `batch_results`, the function:
#'
#' 1. Extracts the `modified_ratio_table` produced during ratio normalization.
#' 2. Ensures the table contains valid data (non-empty and with column names).
#' 3. Performs a dose-response fit via [`fit_drc_3pl()`] (default) or
#'    [`fit_drc_4pl()`] (when `model = "4pl"`), where:
#'    * The first column is assumed to be log(inhibitor concentration)
#'    * Remaining columns are compound responses
#' 4. Stores:
#'    * Full DRC results for each compound
#'    * Summary tables, final tables, and curve quality metrics
#'    * Optional Excel result files per plate
#'
#' After processing all plates, if `generate_reports = TRUE`, a consolidated
#' Excel file (`batch_drc_analysis_report.xlsx`) is created, containing:
#'
#' * **Summary** sheet - plate-level statistics
#' * **All_Results** - merged results for all compounds across plates
#' * **Curve_Quality** - QC-oriented summary (R2, slope, quality flag)
#' * A **_summary** and **_final_summary** sheet per plate
#'
#'
#' @return
#' A named list of DRC analysis results for each plate.
#' Each entry contains:
#'
#' \describe{
#'   \item{plate_info}{A list with metadata: `data_file`, `info_sheet`, `sheet_number`.}
#'   \item{drc_result}{Full result object from `fit_drc_3pl()`  or `fit_drc_4pl()` (depending on `model`),
#'   including fitted parameters, summary tables, final summary tables, QC metrics, and counts.}
#'   \item{timestamp}{Time when the plate was processed.}
#' }
#'
#' If `generate_reports = TRUE`, also saves:
#'
#' * `batch_drc_analysis_report.xlsx`
#' * Individual Excel files: `drc_results_<plate>.xlsx`
#'
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_drc")
#' dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
#' invisible(file.copy(
#'   list.files(extdata_dir, pattern = "^nanobret_", full.names = TRUE),
#'   work_dir, overwrite = TRUE
#' ))
#'
#' ratio_res <- batch_ratio_analysis(
#'   directory        = work_dir,
#'   info_file        = "nanobret_info.xlsx",
#'   data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
#'   control_0perc    = "1",
#'   control_100perc  = "24",
#'   selected_columns = 1:24,
#'   generate_reports = FALSE,
#'   output_dir       = tempdir(),
#'   verbose          = FALSE
#' )
#'
#' drc_res <- batch_drc_analysis(
#'   batch_results    = ratio_res,
#'   model            = "4pl",
#'   normalize        = FALSE,
#'   output_dir       = tempdir(),
#'   generate_reports = FALSE,
#'   verbose          = FALSE
#' )
#' drc_res$drc_results$plate_01$drc_result$successful_fits
#' }
#' @seealso
#' * [`fit_drc_3pl()`] - Fits the dose-response curve for a single plate.
#' * [`fit_drc_4pl()`] - 4-parameter logistic fit for a single plate.
#' * [`batch_ratio_analysis()`] - Preprocessing step generating modified ratio tables.
#' * `openxlsx` - Excel manipulation used for reporting.
#'
#'
#' @export

batch_drc_analysis <- function(batch_results,
                               normalize = FALSE,
                               enforce_bottom_threshold = FALSE,
                               bottom_threshold = 60,
                               r_sqr_threshold = 0.8,
                               output_dir = NULL,
                               generate_reports = TRUE,
                               model = "3pl",
                               nd_if_activation = FALSE,
                               verbose = TRUE,
                               hook_effect    = FALSE,
                               hook_threshold = 2) {
  
  # ============================================================================
  # 1. SETUP & DEPENDENCIES
  # ============================================================================
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (generate_reports && !requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required.")
  
  # Helper: Safe fallback operator
  `%||%` <- function(a, b) {
    if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
  }
  
  # Helper: Filename sanitization
  sanitize_filename <- function(name, max_len = 50) {
    clean <- gsub("[^A-Za-z0-9_.-]", "_", name)
    clean <- gsub("_+", "_", clean)
    if (nchar(clean) > max_len) clean <- substr(clean, 1, max_len)
    return(clean)
  }
  
  # Model validation
  model <- tolower(model)
  if (!model %in% c("3pl", "4pl"))
    stop("model must be either '3pl' or '4pl'.")

  # hook_effect / hook_threshold validation
  if (!isFALSE(hook_effect) && !isTRUE(hook_effect) && !is.character(hook_effect))
    stop("hook_effect must be FALSE, TRUE, or a character vector of compound names.")
  if (is.character(hook_effect) && length(hook_effect) == 0L)
    stop("hook_effect character vector must contain at least one compound name.")
  if (!is.numeric(hook_threshold) || length(hook_threshold) != 1L ||
      !is.finite(hook_threshold) || hook_threshold <= 0)
    stop("hook_threshold must be a single positive finite number.")
  
  # Input validation
  if (!is.list(batch_results) || length(batch_results) == 0) {
    stop("batch_results must be a non-empty list.")
  }
  
  # Detect the assay type from the source attribute stamped by
  # batch_ratio_analysis() ("nanobret") or batch_viability_analysis()
  # ("viability"). Falls back to "nanobret" with a warning for manually
  # constructed lists that carry no such attribute.
  detected_assay <- attr(batch_results, "assay_source")
  if (!is.null(detected_assay)) {
    assay_type <- detected_assay
    if (verbose)
      message("Assay type auto-detected: \"", assay_type, "\".")
  } else {
    warning("Could not auto-detect assay type from batch_results ",
            "(no 'assay_source' attribute found); defaulting to \"nanobret\". ",
            "Ensure batch_results comes from batch_ratio_analysis() or ",
            "batch_viability_analysis().")
    assay_type <- "nanobret"
  }
  
  # Directory setup
  if (is.null(output_dir)) output_dir <- getwd()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  detailed_dir <- file.path(output_dir, "Detailed_Reports")
  if (generate_reports && !dir.exists(detailed_dir)) {
    dir.create(detailed_dir, recursive = TRUE)
  }
  
  # ============================================================================
  # 2. INTERNAL FUNCTIONS
  # ============================================================================
  
  # ---------------------------------------------------------------------------
  # Hook-effect detection helpers
  # ---------------------------------------------------------------------------

  # Determine whether hook detection should be applied to a given compound.
  # base_name: the Construct:Compound string (replicate suffix already stripped).
  .apply_hook <- function(base_name) {
    if (isFALSE(hook_effect)) return(FALSE)
    if (isTRUE(hook_effect))  return(TRUE)
    # character vector: check membership
    base_name %in% hook_effect
  }

  # Warn about hook_effect names that match no compound in the batch.
  # Called once after all plates are processed.
  .warn_unmatched_hook_names <- function(seen_base_names) {
    if (!is.character(hook_effect) || length(hook_effect) == 0L) return(invisible(NULL))
    unmatched <- setdiff(hook_effect, seen_base_names)
    if (length(unmatched) > 0L)
      warning("hook_effect: the following compound names were not found in any plate: ",
              paste(unmatched, collapse = ", "), call. = FALSE)
  }

  # Detect and exclude hook-effect points for one compound result.
  # Returns a list:
  #   $result      - updated detailed_result (second-pass fit if hooks found)
  #   $data_full   - data.frame(log_inhibitor, response, excluded)
  #   $excluded_conc - numeric vector of excluded log10 concentrations (may be length 0)
  .detect_and_refit_hook <- function(result, data_table, model, normalize,
                                     enforce_bottom_threshold, bottom_threshold,
                                     r_sqr_threshold, assay_type) {

    # Build full long-format data frame for this compound (all replicates)
    comp_base  <- sub("\\.\\d+$", "", result$compound %||% "")
    col_names  <- colnames(data_table)[-1]
    base_names <- sub("\\.\\d+$", "", col_names)
    comp_cols  <- which(base_names == comp_base) + 1L  # +1 for conc col

    if (length(comp_cols) == 0L) {
      # Fallback: use result$data if column matching fails
      df_full <- result$data
      if (is.null(df_full)) {
        return(list(result = result, data_full = NULL, excluded_conc = numeric(0)))
      }
    } else {
      df_full <- stats::na.omit(data.frame(
        log_inhibitor = rep(data_table[[1]], length(comp_cols)),
        response      = as.numeric(unlist(data_table[, comp_cols, drop = FALSE]))
      ))
    }

    # First-pass Syx
    syx <- result$goodness_of_fit$Syx %||% NA_real_
    if (is.na(syx) || syx <= 0 || !is.finite(syx)) {
      df_full$excluded <- FALSE
      return(list(result = result, data_full = df_full, excluded_conc = numeric(0)))
    }

    # Identify the two highest unique concentration levels
    conc_levels <- sort(unique(df_full$log_inhibitor), decreasing = TRUE)
    candidates  <- conc_levels[seq_len(min(2L, length(conc_levels)))]

    # Compute mean residual at each candidate concentration
    model_obj <- result$model
    if (is.null(model_obj)) {
      df_full$excluded <- FALSE
      return(list(result = result, data_full = df_full, excluded_conc = numeric(0)))
    }

    excluded_conc <- numeric(0)
    for (k in seq_along(candidates)) {
      c_val    <- candidates[k]
      mean_obs <- mean(df_full$response[df_full$log_inhibitor == c_val], na.rm = TRUE)
      fitted_v <- tryCatch(
        stats::predict(model_obj, newdata = data.frame(log_inhibitor = c_val)),
        error = function(e) NA_real_
      )
      if (is.na(fitted_v)) break
      residual <- mean_obs - fitted_v
      if (residual > hook_threshold * syx) {
        excluded_conc <- c(excluded_conc, c_val)
      } else {
        # Stop at first non-hook (working inward from highest conc)
        break
      }
    }

    # Annotate full data frame
    df_full$excluded <- df_full$log_inhibitor %in% excluded_conc

    if (length(excluded_conc) == 0L) {
      return(list(result = result, data_full = df_full, excluded_conc = numeric(0)))
    }

    # Second-pass fit: mask excluded rows in data_table
    masked_table <- data_table
    log_col      <- data_table[[1]]
    for (j in comp_cols) {
      above <- !is.na(log_col) & log_col %in% excluded_conc
      masked_table[above, j] <- NA_real_
    }

    second_pass <- tryCatch({
      if (model == "3pl") {
        fit_drc_3pl(
          data                     = masked_table[, c(1L, comp_cols), drop = FALSE],
          output_file              = NULL,
          normalize                = normalize,
          verbose                  = FALSE,
          enforce_bottom_threshold = enforce_bottom_threshold,
          bottom_threshold         = bottom_threshold,
          r_sqr_threshold          = r_sqr_threshold
        )
      } else {
        fit_drc_4pl(
          data                     = masked_table[, c(1L, comp_cols), drop = FALSE],
          output_file              = NULL,
          normalize                = normalize,
          verbose                  = FALSE,
          enforce_bottom_threshold = enforce_bottom_threshold,
          bottom_threshold         = bottom_threshold,
          r_sqr_threshold          = r_sqr_threshold,
          assay_type               = assay_type
        )
      }
    }, error = function(e) NULL)

    if (is.null(second_pass) ||
        length(second_pass$detailed_results) == 0L ||
        !isTRUE(second_pass$detailed_results[[1]]$success)) {
      warning("hook_effect: second-pass fit failed for '", comp_base,
              "' after excluding ", length(excluded_conc), " point(s). ",
              "Keeping first-pass result.", call. = FALSE)
      df_full$excluded <- FALSE
      return(list(result = result, data_full = df_full, excluded_conc = numeric(0)))
    }

    new_result <- second_pass$detailed_results[[1]]
    new_result$hook_excluded_conc  <- excluded_conc
    new_result$hook_first_pass_syx <- syx

    list(result = result, data_full = df_full, excluded_conc = excluded_conc,
         new_result = new_result)
  }

  # Accumulator for all compound base names seen (for unmatched-name warning)
  .seen_hook_names <- character(0)

  generate_drc_batch_report <- function(drc_results, batch_results, main_dir, sub_dir, verbose = TRUE) {
    
    # Helper for safe Excel sheet names
    get_safe_sheet_name <- function(base_name, suffix = "", existing_names = c()) {
      max_len <- 31 - nchar(suffix) - 3
      clean_base <- gsub("[^A-Za-z0-9_]", "_", base_name)
      candidate_base <- substr(clean_base, 1, max_len)
      
      final_name <- paste0(candidate_base, suffix)
      counter <- 1
      while (final_name %in% existing_names) {
        final_name <- paste0(candidate_base, "_", counter, suffix)
        counter <- counter + 1
      }
      return(final_name)
    }
    
    path_pharma <- file.path(main_dir, "Pharmacology_Summary.xlsx")
    path_details <- file.path(sub_dir, "Batch_Analysis_Details.xlsx")
    
    wb_pharma <- openxlsx::createWorkbook()
    wb_details <- openxlsx::createWorkbook()
    
    # Accumulator lists
    summary_list <- list()
    all_results_list <- list()
    quality_list <- list()
    pharm_list <- list()
    pharm_ci_list <- list()
    
    used_sheet_names <- c("Summary", "All_Results", "Curve_Quality")
    
    # --- REPORT GENERATION LOOP ---
    for (plate_name in names(drc_results)) {
      plate_res_obj <- drc_results[[plate_name]]
      
      if (is.null(plate_res_obj$drc_result) || is.null(plate_res_obj$drc_result$summary_table)) next
      drc_summary <- plate_res_obj$drc_result$summary_table
      if (nrow(drc_summary) == 0) next
      
      # --- Access full data table (with replicates) ---
      # Use the table used for fitting (Raw or Normalized) to maintain scale consistency
      used_norm <- plate_res_obj$data_tables$fitting_normalization %||% FALSE
      
      full_data_df <- if (used_norm) {
        plate_res_obj$data_tables$normalized_data
      } else {
        plate_res_obj$data_tables$raw_data
      }
      
      # 1. Standard Summary
      n_compounds <- nrow(drc_summary)
      successful_fits <- sum(!is.na(drc_summary$IC50))
      
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Plate_Name = plate_name,
        Data_File = plate_res_obj$plate_info$data_file,
        N_Compounds = n_compounds,
        Successful_Fits = successful_fits,
        Success_Rate = round(successful_fits / n_compounds * 100, 1),
        Avg_R_squared = round(mean(drc_summary$R_squared, na.rm = TRUE), 3),
        Good_Curves = sum(grepl("Good curve", drc_summary$Curve_Quality, ignore.case = TRUE)),
        stringsAsFactors = FALSE
      )
      
      # 2. All Results
      drc_summary_copy <- drc_summary
      drc_summary_copy$Plate <- plate_name
      cols_order <- c("Plate", setdiff(names(drc_summary_copy), "Plate"))
      drc_summary_copy <- drc_summary_copy[, cols_order, drop = FALSE]
      # Coerce mixed-type columns to character so bind_rows works across plates
      # (plates with N/D compounds have character; plates without have numeric)
      mixed_cols <- c("LogIC50", "IC50",
                      "LogIC50_Lower_95CI", "LogIC50_Upper_95CI",
                      "IC50_Lower_95CI", "IC50_Upper_95CI")
      for (.col in intersect(mixed_cols, names(drc_summary_copy))) {
        drc_summary_copy[[.col]] <- as.character(drc_summary_copy[[.col]])
      }
      all_results_list[[length(all_results_list) + 1]] <- drc_summary_copy
      
      # 3. Quality List
      desired_cols <- c("Plate", "Compound", "Curve_Quality", "R_squared", "Max_Slope",
                        "Ideal_Hill_Slope", "HillSlope",
                        "HillSlope_Lower_95CI", "HillSlope_Upper_95CI",
                        "Bottom", "Top", "LogIC50", "IC50")
      existing_cols <- intersect(desired_cols, names(drc_summary_copy))
      quality_list[[length(quality_list) + 1]] <- drc_summary_copy[, existing_cols, drop = FALSE]
      
      # 4. PHARMACOLOGY SUMMARY
      res_root <- plate_res_obj$drc_result
      detailed_res <- res_root$detailed_results %||% res_root$curve_results %||% res_root$fits
      if (!is.list(detailed_res)) detailed_res <- list()
      
      # Build a per-compound outlier count lookup for this plate.
      # outliers_replaced is an attribute on modified_ratio_table set by
      # rout_outliers_batch(). It lives on batch_results (the input to
      # batch_drc_analysis), not on drc_results (the DRC fit output).
      # Its $column field holds the full column name (e.g. "KinaseA:Cpd1.2");
      # we strip the replicate suffix to get the base compound name and count
      # rows per base name. When outlier detection was not run -> all 0.
      mrt_for_outliers <- batch_results[[plate_name]]$result$modified_ratio_table
      or_attr <- attr(mrt_for_outliers, "outliers_replaced")
      outlier_counts <- if (!is.null(or_attr) && nrow(or_attr) > 0L) {
        base_names <- sub("\\.\\d+$", "", or_attr$column)
        tapply(seq_len(nrow(or_attr)), base_names, length)
      } else {
        integer(0)
      }
      
      if (length(detailed_res) > 0) {
        
        for (i in seq_along(detailed_res)) {
          res <- detailed_res[[i]]
          if (is.null(res$success) || !isTRUE(res$success)) next
          
          # Name parsing
          clean_name <- gsub("\\.\\d+$", "", res$compound %||% "Unknown")
          parts <- strsplit(clean_name, " \\| |:")[[1]]
          parts <- trimws(parts)
          construct_name <- parts[1]
          compound_name  <- if (length(parts) > 1) parts[2] else parts[1]
          
          # --- Highest tested concentration (from concentration column of data table) ---
          # Rounded to the nearest integer for clean display (e.g. 24.55 uM -> 25 uM).
          highest_conc_uM <- NA_real_
          if (!is.null(full_data_df) && nrow(full_data_df) >= 2) {
            log_concs <- suppressWarnings(as.numeric(full_data_df[, 1]))
            log_concs <- log_concs[!is.na(log_concs)]
            if (length(log_concs) > 0)
              highest_conc_uM <- round(max(10^log_concs * 1e6))
          }
          
          # --- pIC50 ---
          log_ic50 <- NA_real_
          if (!is.null(res$parameters) && length(res$parameters$Value) >= 3) {
            log_ic50 <- res$parameters$Value[3]
          }
          pic50    <- if (!is.na(log_ic50)) -log_ic50 else NA_real_
          ic50_uM  <- if (!is.na(log_ic50)) 10^log_ic50 * 1e6  else NA_real_
          ic50_nM  <- if (!is.na(log_ic50)) 10^log_ic50 * 1e9  else NA_real_
          
          # --- N/D: replace IC50 and pIC50 with "N/D" for flat (always) or
          #         activation (when nd_if_activation = TRUE) curves ---
          res_curve_type <- res$curve_type %||% "unknown"
          is_nd <- (res_curve_type == "flat") ||
            (nd_if_activation && res_curve_type == "activation")
          
          # --- IC50 display: replace with ">highest" if IC50 exceeds tested range ---
          ic50_above_range <- !is.na(ic50_uM) && !is.na(highest_conc_uM) &&
            ic50_uM > highest_conc_uM
          ic50_uM_display <- if (ic50_above_range) {
            sprintf(">%g", highest_conc_uM)
          } else if (!is.na(ic50_uM)) {
            as.character(round(ic50_uM, 3))
          } else {
            NA_character_
          }
          ic50_nM_display <- if (ic50_above_range) {
            sprintf(">%g", highest_conc_uM * 1e3)
          } else if (!is.na(ic50_nM)) {
            as.character(round(ic50_nM, 3))
          } else {
            NA_character_
          }
          
          # --- CI ---
          ci_log_lower_bound <- NA_real_
          ci_log_upper_bound <- NA_real_
          if (!is.null(res$confidence_intervals) && !is.null(res$confidence_intervals$LogIC50)) {
            ci_log_lower_bound <- res$confidence_intervals$LogIC50[1]
            ci_log_upper_bound <- res$confidence_intervals$LogIC50[2]
          }
          
          pic50_diff_upper <- NA_real_
          pic50_diff_lower <- NA_real_
          if (!is.na(pic50) && !is.na(ci_log_lower_bound) && !is.na(ci_log_upper_bound)) {
            abs_pic50_upper <- -ci_log_lower_bound
            abs_pic50_lower <- -ci_log_upper_bound
            pic50_diff_upper <- abs_pic50_upper - pic50
            pic50_diff_lower <- pic50 - abs_pic50_lower
          }
          
          # --- NORMALIZED SPAN CALCULATION ---
          span_ratio <- NA_real_
          
          if (!is.null(full_data_df) && nrow(full_data_df) >= 2) {
            
            # Locate all replicate columns for this compound by base name.
            # This works for any number of replicates (2, 4, 6, ...) produced
            # by merge_plate_replicates() or a standard single-plate run.
            data_col_names <- colnames(full_data_df)[-1]          # exclude conc col
            base_col_names <- sub("\\.\\d+$", "", data_col_names)  # strip .2/.3/.4
            compound_base  <- res$compound %||% ""
            rep_positions  <- which(base_col_names == compound_base) + 1L  # +1 for conc col
            
            if (length(rep_positions) > 0) {
              
              # 1. Mean across ALL replicates at the FIRST concentration row
              mean_start_row <- mean(
                as.numeric(unlist(full_data_df[1, rep_positions])), na.rm = TRUE)
              
              # 2. Mean across ALL replicates at the LAST concentration row
              last_idx <- nrow(full_data_df)
              mean_end_row <- mean(
                as.numeric(unlist(full_data_df[last_idx, rep_positions])), na.rm = TRUE)
              
              # 3. Fitted Span (Parameter 5)
              fit_span <- NA_real_
              if (length(res$parameters$Value) >= 5) fit_span <- res$parameters$Value[5]
              
              # 4. Calculation: abs(fit_span) / abs(mean_end_row - mean_start_row)
              diff_window <- abs(mean_end_row - mean_start_row)
              
              if (!is.na(diff_window) && diff_window > 1e-6 && !is.na(fit_span)) {
                span_ratio <- abs(fit_span) / diff_window
              }
            }
          }
          
          # For 3PL: ideal_hill_slope is stored directly on the result.
          # For 4PL: it is NULL; fall back to the estimated HillSlope (Value[4]).
          ideal_hill <- if (!is.null(res$ideal_hill_slope)) {
            res$ideal_hill_slope
          } else if (!is.null(res$parameters) && length(res$parameters$Value) >= 4) {
            res$parameters$Value[4]
          } else {
            NA_real_
          }
          
          # --- WARNING AND EXCLUSION FLAGS ---
          warning_collector <- character()
          exclusion_collector <- character()
          
          # CI Analysis
          if (is.na(pic50_diff_lower) || is.na(pic50_diff_upper)) {
            exclusion_collector <- c(exclusion_collector, "Undefined CI")
          } else {
            # 3-fold = log10(3) = 0.47712
            # 5-fold = log10(5) = 0.69897
            ci_warnings <- character()
            ci_exclusions <- character()
            
            if (pic50_diff_lower > 0.69897) {
              ci_exclusions <- c(ci_exclusions, sprintf("Lower CI >5-fold (%.3f)", pic50_diff_lower))
            } else if (pic50_diff_lower > 0.47712) {
              ci_warnings <- c(ci_warnings, sprintf("Lower CI >3-fold (%.3f)", pic50_diff_lower))
            }
            
            if (pic50_diff_upper > 0.69897) {
              ci_exclusions <- c(ci_exclusions, sprintf("Upper CI >5-fold (%.3f)", pic50_diff_upper))
            } else if (pic50_diff_upper > 0.47712) {
              ci_warnings <- c(ci_warnings, sprintf("Upper CI >3-fold (%.3f)", pic50_diff_upper))
            }
            
            if (length(ci_warnings) > 0) {
              warning_collector <- c(warning_collector, paste(ci_warnings, collapse = "; "))
            }
            if (length(ci_exclusions) > 0) {
              exclusion_collector <- c(exclusion_collector, paste(ci_exclusions, collapse = "; "))
            }
          }
          
          # Hill Slope Analysis
          if (!is.na(ideal_hill)) {
            curve_type <- res$curve_type %||% "unknown"
            hill_message <- ""
            
            if (curve_type == "activation") {
              if (ideal_hill < 0.5 || ideal_hill > 1.5) {
                hill_message <- sprintf("Hill Slope (expected 0.5-1.5): %.3f", ideal_hill)
              }
            } else {
              if (ideal_hill > -0.5 || ideal_hill < -1.5) {
                hill_message <- sprintf("Hill Slope (expected -1.5 to -0.5): %.3f", ideal_hill)
              }
            }
            
            if (nchar(hill_message) > 0) {
              warning_collector <- c(warning_collector, hill_message)
            }
          }
          
          # Normalized Span Analysis
          if (!is.na(span_ratio)) {
            if (span_ratio < 0.5) {
              exclusion_collector <- c(exclusion_collector, sprintf("Norm Span < 0.5 (%.2f)", span_ratio))
            } else if (span_ratio > 1.5) {
              exclusion_collector <- c(exclusion_collector, sprintf("Norm Span > 1.5 (%.2f)", span_ratio))
            }
          }
          
          # Assay Quality Metrics (from batch_results interval_means)
          # interval_means is stored transposed: rows = metrics, cols = constructs.
          # If any of the three quality comments is "insufficient", flag it.
          {
            im <- batch_results[[plate_name]]$result$interval_means
            if (!is.null(im) && construct_name %in% colnames(im)) {
              insuf_metrics <- character()
              luc_q <- if ("Luciferase_signal_comment" %in% rownames(im)) im["Luciferase_signal_comment", construct_name] else NA
              aw_q  <- if ("Assay_window_Comment"      %in% rownames(im)) im["Assay_window_Comment",      construct_name] else NA
              zs_q  <- if ("Assay_z_Comment"           %in% rownames(im)) im["Assay_z_Comment",           construct_name] else NA
              if (!is.na(luc_q) && startsWith(as.character(luc_q), "insufficient"))
                insuf_metrics <- c(insuf_metrics, "Luciferase signal insufficient")
              if (!is.na(aw_q) && startsWith(as.character(aw_q), "insufficient"))
                insuf_metrics <- c(insuf_metrics, "Assay window insufficient")
              if (!is.na(zs_q) && startsWith(as.character(zs_q), "insufficient"))
                insuf_metrics <- c(insuf_metrics, "Z' insufficient")
              if (length(insuf_metrics) > 0)
                exclusion_collector <- c(exclusion_collector, paste(insuf_metrics, collapse = "; "))
            }
          }
          
          # IC50 above tested range -> add to exclusion
          if (ic50_above_range)
            exclusion_collector <- c(exclusion_collector,
                                     sprintf("IC50 above tested range (>%g uM)", highest_conc_uM))
          
          # Set "OK" for empty collectors
          final_warnings <- if (length(warning_collector) > 0) paste(warning_collector, collapse = "; ") else "OK"
          final_exclusions <- if (length(exclusion_collector) > 0) paste(exclusion_collector, collapse = "; ") else "OK"
          
          # Apply N/D for flat (always) or activation (if nd_if_activation = TRUE)
          ic50_uM_final <- if (is_nd) "N/D" else ic50_uM_display
          ic50_nM_final <- if (is_nd) "N/D" else ic50_nM_display
          pic50_final   <- if (is_nd) "N/D" else as.character(round(pic50, 3))
          
          # Count outliers removed for this compound (base name match)
          compound_base_name <- res$compound %||% ""
          n_outliers_removed <- if (length(outlier_counts) > 0L &&
                                    compound_base_name %in% names(outlier_counts)) {
            as.integer(outlier_counts[[compound_base_name]])
          } else {
            0L
          }
          
          # Hook-effect: extract excluded concentrations for this compound
          .hook_excl <- res$hook_excluded_conc %||% numeric(0)
          hook_excluded_str <- if (length(.hook_excl) > 0L) {
            paste(round(.hook_excl, 3), collapse = ", ")
          } else {
            ""
          }

          pharm_list[[length(pharm_list) + 1]] <- data.frame(
            Plate = plate_name,
            Construct = construct_name,
            Compound = compound_name,
            `IC50 (uM)` = ic50_uM_final,
            `IC50 (nM)` = ic50_nM_final,
            pIC50 = pic50_final,
            check.names = FALSE,
            CI_95_Upper = as.character(round(pic50_diff_upper, 3)),
            CI_95_Lower = as.character(round(pic50_diff_lower, 3)),
            Ideal_Hill_Slope = as.character(round(ideal_hill, 3)),
            Normalized_Span = as.character(round(span_ratio, 3)),
            Outliers_Removed = n_outliers_removed,
            Warning = final_warnings,
            Exclusion = final_exclusions,
            Hook_Excluded_Conc = hook_excluded_str,
            stringsAsFactors = FALSE
          )

          # --- CI95% SUMMARY ROW ---
          # Convert LogIC50 CI bounds back to uM/nM for display
          ci_uM_lower <- if (!is.na(ci_log_lower_bound)) round(10^ci_log_lower_bound * 1e6, 3) else NA_real_
          ci_uM_upper <- if (!is.na(ci_log_upper_bound)) round(10^ci_log_upper_bound * 1e6, 3) else NA_real_
          ci_nM_lower <- if (!is.na(ci_log_lower_bound)) round(10^ci_log_lower_bound * 1e9, 3) else NA_real_
          ci_nM_upper <- if (!is.na(ci_log_upper_bound)) round(10^ci_log_upper_bound * 1e9, 3) else NA_real_

          # Format IC50 (uM) with CI: "0.128 (0.080 - 0.204)"
          # Note: lower LogIC50 -> higher IC50, so min(ci_uM_lower, ci_uM_upper) is the true lower bound
          ic50_uM_ci_display <- if (is_nd) {
            "N/D"
          } else if (ic50_above_range) {
            sprintf(">%g (>%g - >%g)", highest_conc_uM, highest_conc_uM, highest_conc_uM)
          } else if (!is.na(ic50_uM) && !is.na(ci_uM_lower) && !is.na(ci_uM_upper)) {
            sprintf("%s (%s - %s)",
                    as.character(round(ic50_uM, 3)),
                    as.character(min(ci_uM_lower, ci_uM_upper)),
                    as.character(max(ci_uM_lower, ci_uM_upper)))
          } else if (!is.na(ic50_uM)) {
            sprintf("%s (NA)", as.character(round(ic50_uM, 3)))
          } else {
            NA_character_
          }

          # Format IC50 (nM) with CI
          ic50_nM_ci_display <- if (is_nd) {
            "N/D"
          } else if (ic50_above_range) {
            sprintf(">%g (>%g - >%g)", highest_conc_uM * 1e3, highest_conc_uM * 1e3, highest_conc_uM * 1e3)
          } else if (!is.na(ic50_nM) && !is.na(ci_nM_lower) && !is.na(ci_nM_upper)) {
            sprintf("%s (%s - %s)",
                    as.character(round(ic50_nM, 3)),
                    as.character(min(ci_nM_lower, ci_nM_upper)),
                    as.character(max(ci_nM_lower, ci_nM_upper)))
          } else if (!is.na(ic50_nM)) {
            sprintf("%s (NA)", as.character(round(ic50_nM, 3)))
          } else {
            NA_character_
          }

          pharm_ci_list[[length(pharm_ci_list) + 1]] <- data.frame(
            Plate = plate_name,
            Construct = construct_name,
            Compound = compound_name,
            `IC50 (uM) (CI95%)` = ic50_uM_ci_display,
            `IC50 (nM) (CI95%)` = ic50_nM_ci_display,
            check.names = FALSE,
            Ideal_Hill_Slope = as.character(round(ideal_hill, 3)),
            Outliers_Removed = n_outliers_removed,
            Warning = final_warnings,
            Exclusion = final_exclusions,
            Hook_Excluded_Conc = hook_excluded_str,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    # Consolidation
    summary_data <- if (length(summary_list) > 0) dplyr::bind_rows(summary_list) else data.frame()
    all_results_combined <- if (length(all_results_list) > 0) dplyr::bind_rows(all_results_list) else data.frame()
    quality_combined <- if (length(quality_list) > 0) dplyr::bind_rows(quality_list) else data.frame()
    pharm_combined <- if (length(pharm_list) > 0) dplyr::bind_rows(pharm_list) else data.frame()
    pharm_ci_combined <- if (length(pharm_ci_list) > 0) dplyr::bind_rows(pharm_ci_list) else data.frame()
    
    # Remove placeholder rows where Construct or Compound is exactly NA, NA_2,
    # NA_3, etc. (artifacts from unnamed columns). Uses a strict regex so that
    # real gene names containing "NA" (e.g. NAGA, CANAL) are never removed.
    na_pattern <- "^NA(_\\d+)?$"
    if (nrow(pharm_combined) > 0) {
      is_na_construct <- is.na(pharm_combined$Construct) |
        grepl(na_pattern, pharm_combined$Construct)
      is_na_compound  <- is.na(pharm_combined$Compound)  |
        grepl(na_pattern, pharm_combined$Compound)
      pharm_combined  <- pharm_combined[!(is_na_construct | is_na_compound), , drop = FALSE]
    }
    if (nrow(pharm_ci_combined) > 0) {
      is_na_construct_ci <- is.na(pharm_ci_combined$Construct) |
        grepl(na_pattern, pharm_ci_combined$Construct)
      is_na_compound_ci  <- is.na(pharm_ci_combined$Compound)  |
        grepl(na_pattern, pharm_ci_combined$Compound)
      pharm_ci_combined  <- pharm_ci_combined[!(is_na_construct_ci | is_na_compound_ci), , drop = FALSE]
    }
    
    # --- EXCEL WRITING ---
    openxlsx::addWorksheet(wb_pharma, "Pharmacology_Summary")
    if (nrow(pharm_combined) > 0) {
      openxlsx::writeData(wb_pharma, "Pharmacology_Summary", pharm_combined)
    } else {
      openxlsx::writeData(wb_pharma, "Pharmacology_Summary", data.frame(Note = "No pharmacology data available"))
    }
    openxlsx::addWorksheet(wb_pharma, "Pharmacology_Summary_CI")
    if (nrow(pharm_ci_combined) > 0) {
      openxlsx::writeData(wb_pharma, "Pharmacology_Summary_CI", pharm_ci_combined)
    } else {
      openxlsx::writeData(wb_pharma, "Pharmacology_Summary_CI", data.frame(Note = "No pharmacology data available"))
    }
    openxlsx::saveWorkbook(wb_pharma, path_pharma, overwrite = TRUE)
    
    openxlsx::addWorksheet(wb_details, "Summary")
    openxlsx::writeData(wb_details, "Summary", summary_data)
    openxlsx::addWorksheet(wb_details, "All_Results")
    if (nrow(all_results_combined) > 0) openxlsx::writeData(wb_details, "All_Results", all_results_combined)
    openxlsx::addWorksheet(wb_details, "Curve_Quality")
    if (nrow(quality_combined) > 0) openxlsx::writeData(wb_details, "Curve_Quality", quality_combined)
    
    for (plate_name in names(drc_results)) {
      plate_res_obj <- drc_results[[plate_name]]
      if (is.null(plate_res_obj$drc_result$summary_table)) next
      
      sheet_sum <- get_safe_sheet_name(plate_name, "_sum", used_sheet_names)
      used_sheet_names <- c(used_sheet_names, sheet_sum)
      openxlsx::addWorksheet(wb_details, sheet_sum)
      openxlsx::writeData(wb_details, sheet_sum, plate_res_obj$drc_result$summary_table)
      
      if (!is.null(plate_res_obj$data_tables$raw_data)) {
        sheet_raw <- get_safe_sheet_name(plate_name, "_raw", used_sheet_names)
        used_sheet_names <- c(used_sheet_names, sheet_raw)
        openxlsx::addWorksheet(wb_details, sheet_raw)
        openxlsx::writeData(wb_details, sheet_raw, plate_res_obj$data_tables$raw_data)
      }
      
      if (!is.null(plate_res_obj$data_tables$normalized_data)) {
        sheet_norm <- get_safe_sheet_name(plate_name, "_norm", used_sheet_names)
        used_sheet_names <- c(used_sheet_names, sheet_norm)
        openxlsx::addWorksheet(wb_details, sheet_norm)
        openxlsx::writeData(wb_details, sheet_norm, plate_res_obj$data_tables$normalized_data)
      }
    }
    
    openxlsx::saveWorkbook(wb_details, path_details, overwrite = TRUE)
    
    if (verbose) {
      message("Reports generated successfully:")
      message("  1. ", path_pharma)
      message("  2. ", path_details)
    }
  }
  
  # Robust data extraction function
  extract_data_for_drc <- function(plate_result) {
    search_locations <- list()
    if (!is.null(plate_result$result)) {
      search_locations <- c(search_locations, list(
        plate_result$result$modified_ratio_table,
        plate_result$result$processed_data,
        plate_result$result$ratio_table,
        plate_result$result$normalized_data,
        plate_result$result$data
      ))
    }
    search_locations <- c(search_locations, list(
      plate_result$modified_ratio_table,
      plate_result$processed_data,
      plate_result$ratio_table,
      plate_result$normalized_data,
      plate_result$data
    ))
    if (is.data.frame(plate_result) && nrow(plate_result) > 0) {
      search_locations <- c(search_locations, list(plate_result))
    }
    for (dt in search_locations) {
      if (!is.null(dt) && is.data.frame(dt) && nrow(dt) > 0 && ncol(dt) >= 3) {
        first_col <- dt[[1]]
        if (is.numeric(first_col) ||
            (is.character(first_col) &&
             all(!is.na(suppressWarnings(as.numeric(first_col[!is.na(first_col)])))))) {
          return(dt)
        }
      }
    }
    return(NULL)
  }
  
  # ============================================================================
  # 3. MAIN EXECUTION
  # ============================================================================
  
  if (verbose) {
    message("==========================================================")
    message("STARTING BATCH DOSE-RESPONSE ANALYSIS")
    message("==========================================================")
    message("Model: ", toupper(model), " (", if (model == "3pl") "Hill slope fixed at +/-1" else "Hill slope freely estimated", ")")
    message("Assay type: ", assay_type)
    message("Main Output: ", output_dir)
  }
  
  drc_results <- list()
  failed_plates <- character()
  total_plates <- length(batch_results)
  
  for (i in seq_along(batch_results)) {
    plate_name <- names(batch_results)[i]
    if (verbose) message(sprintf("\nProcessing %d/%d: %s", i, total_plates, plate_name))
    
    tryCatch({
      proc_start <- Sys.time()
      
      data_table <- extract_data_for_drc(batch_results[[plate_name]])
      if (is.null(data_table)) stop("No valid data table found in plate result.")
      
      output_file <- NULL
      if (generate_reports) {
        clean_name <- sanitize_filename(plate_name)
        output_file <- file.path(detailed_dir, paste0("drc_", clean_name, ".xlsx"))
      }
      
      plate_drc_result <- tryCatch({
        if (model == "3pl") {
          fit_drc_3pl(
            data                     = data_table,
            output_file              = NULL,   # Excel written below after N/D post-processing
            normalize                = normalize,
            verbose                  = FALSE,
            enforce_bottom_threshold = enforce_bottom_threshold,
            bottom_threshold         = bottom_threshold,
            r_sqr_threshold          = r_sqr_threshold
          )
        } else {
          fit_drc_4pl(
            data                     = data_table,
            output_file              = NULL,   # Excel written below after N/D post-processing
            normalize                = normalize,
            verbose                  = FALSE,
            enforce_bottom_threshold = enforce_bottom_threshold,
            bottom_threshold         = bottom_threshold,
            r_sqr_threshold          = r_sqr_threshold,
            assay_type               = assay_type
          )
        }
      }, error = function(e) {
        return(list(
          successful_fits = 0,
          n_compounds = ncol(data_table) - 1,
          summary_table = data.frame(),
          detailed_results = list(),
          error = e$message
        ))
      })
      
      # -- Normalise summary_table columns for downstream compatibility ------
      # 3PL produces `Ideal_Hill_Slope`; 4PL produces `HillSlope`.
      # Add an `Ideal_Hill_Slope` alias in the 4PL case so that the report
      # generator (which references that column by name) works for both models.
      if (model == "4pl" &&
          !is.null(plate_drc_result$summary_table) &&
          nrow(plate_drc_result$summary_table) > 0 &&
          "HillSlope" %in% names(plate_drc_result$summary_table) &&
          !"Ideal_Hill_Slope" %in% names(plate_drc_result$summary_table)) {
        plate_drc_result$summary_table$Ideal_Hill_Slope <-
          plate_drc_result$summary_table$HillSlope
      }
      
      # -- Apply N/D to Summary and Final_Summary for flat (always) and
      # activation (when nd_if_activation = TRUE) curves --------------------
      if (!is.null(plate_drc_result$detailed_results) &&
          length(plate_drc_result$detailed_results) > 0 &&
          !is.null(plate_drc_result$summary_table) &&
          nrow(plate_drc_result$summary_table) > 0) {
        
        # Build named logical vector: is_nd per compound
        is_nd_vec <- vapply(plate_drc_result$detailed_results, function(res) {
          ct <- res$curve_type %||% "unknown"
          (ct == "flat") || (nd_if_activation && ct == "activation")
        }, logical(1L))
        names(is_nd_vec) <- vapply(plate_drc_result$detailed_results, function(res) {
          strsplit(res$compound %||% "", " \\| ")[[1]][1]
        }, character(1L))
        
        # Columns to set to "N/D" in summary_table
        nd_cols <- c("LogIC50", "IC50",
                     "LogIC50_Lower_95CI", "LogIC50_Upper_95CI",
                     "IC50_Lower_95CI",    "IC50_Upper_95CI")
        nd_cols_present <- intersect(nd_cols, names(plate_drc_result$summary_table))
        
        if (length(nd_cols_present) > 0 && any(is_nd_vec)) {
          for (col in nd_cols_present) {
            plate_drc_result$summary_table[[col]] <-
              as.character(plate_drc_result$summary_table[[col]])
          }
          for (cpd in names(is_nd_vec)[is_nd_vec]) {
            row_idx <- which(plate_drc_result$summary_table$Compound == cpd)
            if (length(row_idx) > 0) {
              plate_drc_result$summary_table[row_idx, nd_cols_present] <- "N/D"
            }
          }
        }
        
        # Rebuild final_summary_table from updated summary_table
        st <- plate_drc_result$summary_table
        if (nrow(st) > 0) {
          t_data <- as.data.frame(t(st[, -1, drop = FALSE]))
          colnames(t_data) <- st$Compound
          plate_drc_result$final_summary_table <- t_data
        }
      }

      # -- Hook-effect detection and second-pass refitting ----------------
      # Runs after N/D post-processing so curve_type is already set.
      if (!isFALSE(hook_effect) &&
          !is.null(plate_drc_result$detailed_results) &&
          length(plate_drc_result$detailed_results) > 0L) {

        for (.hi in seq_along(plate_drc_result$detailed_results)) {
          .res <- plate_drc_result$detailed_results[[.hi]]

          # Only process successful inhibition fits
          if (!isTRUE(.res$success)) next
          if ((.res$curve_type %||% "unknown") != "inhibition") next

          .base <- sub("\\.\\d+$", "", .res$compound %||% "")
          .seen_hook_names <<- unique(c(.seen_hook_names, .base))

          if (!.apply_hook(.base)) next

          if (verbose)
            message("  [hook] Checking '", .base, "' for hook effect...")

          .hook_out <- .detect_and_refit_hook(
            result                   = .res,
            data_table               = data_table,
            model                    = model,
            normalize                = normalize,
            enforce_bottom_threshold = enforce_bottom_threshold,
            bottom_threshold         = bottom_threshold,
            r_sqr_threshold          = r_sqr_threshold,
            assay_type               = assay_type
          )

          # Attach data_full to the (possibly updated) result
          if (!is.null(.hook_out$new_result)) {
            # Hook(s) found and second pass succeeded
            .hook_out$new_result$data_full <- .hook_out$data_full
            plate_drc_result$detailed_results[[.hi]] <- .hook_out$new_result
            if (verbose)
              message("  [hook] Excluded ", length(.hook_out$excluded_conc),
                      " point(s) at log10[M] = ",
                      paste(round(.hook_out$excluded_conc, 3), collapse = ", "),
                      " for '", .base, "'")
          } else {
            # No hook detected (or second pass failed) - attach data_full anyway
            plate_drc_result$detailed_results[[.hi]]$data_full <- .hook_out$data_full
          }
        }

        # Rebuild summary_table and final_summary_table after potential refits
        if (!is.null(plate_drc_result$summary_table) &&
            nrow(plate_drc_result$summary_table) > 0L) {
          for (.hi in seq_along(plate_drc_result$detailed_results)) {
            .res2 <- plate_drc_result$detailed_results[[.hi]]
            if (!isTRUE(.res2$success)) next
            if (length(.res2$hook_excluded_conc %||% numeric(0)) == 0L) next
            .cpd_name <- strsplit(.res2$compound %||% "", " \\| ")[[1]][1]
            .row_idx  <- which(plate_drc_result$summary_table$Compound == .cpd_name)
            if (length(.row_idx) == 0L) next
            # Update key columns from second-pass parameters
            .p <- .res2$parameters$Value
            if (!is.null(.p) && length(.p) >= 5L) {
              plate_drc_result$summary_table$Bottom[.row_idx]   <- round(.p[1], 3)
              plate_drc_result$summary_table$Top[.row_idx]      <- round(.p[2], 3)
              plate_drc_result$summary_table$LogIC50[.row_idx]  <- round(.p[3], 3)
              plate_drc_result$summary_table$IC50[.row_idx]     <- format(.p[4], scientific = TRUE)
              plate_drc_result$summary_table$Span[.row_idx]     <- round(.p[5], 3)
            }
            .gof <- .res2$goodness_of_fit
            if (!is.null(.gof)) {
              plate_drc_result$summary_table$R_squared[.row_idx] <- round(.gof$R_squared %||% NA_real_, 3)
              plate_drc_result$summary_table$Syx[.row_idx]       <- round(.gof$Syx %||% NA_real_, 3)
            }
          }
          # Rebuild transposed final_summary_table
          .st2 <- plate_drc_result$summary_table
          if (nrow(.st2) > 0L) {
            .td2 <- as.data.frame(t(.st2[, -1, drop = FALSE]))
            colnames(.td2) <- .st2$Compound
            plate_drc_result$final_summary_table <- .td2
          }
        }
      }  # end hook_effect block
      
      # -- Write per-plate Excel file with N/D-corrected tables -------------
      if (!is.null(output_file) && requireNamespace("openxlsx", quietly = TRUE)) {
        tryCatch({
          wb_plate <- openxlsx::createWorkbook()
          
          # Sheet 1: Final_Summary (transposed)
          openxlsx::addWorksheet(wb_plate, "Final_Summary")
          fst <- plate_drc_result$final_summary_table
          if (!is.null(fst) && nrow(fst) > 0) {
            out_final <- cbind(
              data.frame(Parameter = rownames(fst), stringsAsFactors = FALSE),
              fst
            )
            openxlsx::writeData(wb_plate, "Final_Summary", out_final)
          } else {
            openxlsx::writeData(wb_plate, "Final_Summary", "No data")
          }
          
          # Sheet 2: Summary (one row per compound)
          openxlsx::addWorksheet(wb_plate, "Summary")
          openxlsx::writeData(wb_plate, "Summary", plate_drc_result$summary_table)
          
          # Sheet 3: Normalized_Data
          openxlsx::addWorksheet(wb_plate, "Normalized_Data")
          norm_d <- plate_drc_result$normalized_data
          if (!is.null(norm_d)) openxlsx::writeData(wb_plate, "Normalized_Data", norm_d)
          
          # Sheet 4: Original_Data
          openxlsx::addWorksheet(wb_plate, "Original_Data")
          orig_d <- plate_drc_result$original_data
          if (!is.null(orig_d)) openxlsx::writeData(wb_plate, "Original_Data", orig_d)
          
          openxlsx::saveWorkbook(wb_plate, output_file, overwrite = TRUE)
        }, error = function(e) {
          warning("Could not write per-plate Excel for '", plate_name, "': ", e$message)
        })
      }
      
      d_file <- batch_results[[plate_name]]$data_file %||% "unknown"
      i_sheet <- batch_results[[plate_name]]$info_sheet %||% "unknown"
      s_num <- batch_results[[plate_name]]$sheet_number %||% "unknown"
      
      drc_results[[plate_name]] <- list(
        plate_info = list(original_name = plate_name, data_file = d_file, info_sheet = i_sheet, sheet_number = s_num),
        drc_result = plate_drc_result,
        data_tables = list(
          raw_data = plate_drc_result$original_data %||% data_table,
          normalized_data = plate_drc_result$normalized_data,
          fitting_normalization = plate_drc_result$used_normalized_data %||% normalize
        ),
        processing_time = as.numeric(difftime(Sys.time(), proc_start, units = "secs"))
      )
      
      if (verbose) {
        succ <- plate_drc_result$successful_fits %||% 0
        tot <- plate_drc_result$n_compounds %||% 0
        message(sprintf("  -> Success: %d/%d compounds (%.1f sec)", succ, tot, difftime(Sys.time(), proc_start, units = "secs")))
      }
      
    }, error = function(e) {
      warning(sprintf("Failed to process plate '%s': %s", plate_name, e$message))
      failed_plates <<- c(failed_plates, plate_name)
    })
  }
  
  # ============================================================================
  # Warn about any hook_effect compound names that were never matched
  .warn_unmatched_hook_names(.seen_hook_names)

  # 4. REPORTING & RETURN
  # ============================================================================
  
  report_info <- NULL
  if (length(drc_results) > 0 && generate_reports) {
    if (verbose) {
      message("\n", paste(rep("=", 50), collapse = ""))
      message("Generating consolidated reports...")
    }
    tryCatch({
      report_info <- generate_drc_batch_report(drc_results, batch_results, output_dir, detailed_dir, verbose)
    }, error = function(e) {
      warning("Failed to generate master reports: ", e$message)
    })
  }
  
  if (verbose) {
    message("\n", paste(rep("=", 50), collapse = ""))
    message("BATCH ANALYSIS COMPLETE")
    message(paste(rep("=", 50), collapse = ""))
  }
  
  return(invisible(list(
    drc_results = drc_results,
    metadata = list(
      total_plates = total_plates,
      success_count = length(drc_results),
      failed_plates = failed_plates,
      output_directory = output_dir,
      model = model,
      assay_type = assay_type,
      normalize = normalize,
      timestamp = Sys.time()
    ),
    report_info = report_info
  )))
}

