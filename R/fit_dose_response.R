#' Fit a 3-Parameter Logistic Dose-Response Model
#'
#' The `fit_drc_3pl()` function fits a three-parameter logistic (3PL) dose-response
#' model to experimental data. It supports duplicate measurements, normalization,
#' and optional filtering of low bottom values. The function computes fitted parameters,
#' model statistics, and diagnostic metrics, providing a concise summary of model performance.
#'
#' @param data A data frame containing dose-response data. The first column should
#'   contain log10(inhibitor concentration) values, followed by pairs of response
#'   columns (duplicate measurements).
#' @param output_file Optional character string specifying the output file path.
#'   Supports both Excel (.xlsx) and CSV formats.
#' @param normalize Logical indicating whether to normalize response data.
#'   If TRUE, responses are normalized from 0-100\% based on the first and last values.
#'   Default is FALSE.
#' @param verbose Logical indicating whether to display progress messages and
#'   summary statistics. Default is TRUE
#' @param enforce_bottom_threshold Logical indicating whether to exclude IC50 values
#'   for curves where the Bottom parameter exceeds a threshold. Default is FALSE.
#' @param bottom_threshold Numeric value specifying the Bottom threshold when
#'   `enforce_bottom_threshold = TRUE`. Curves with Bottom <= this value will have
#'   IC50 values set to NA. Default is 60.
#' @param r_sqr_threshold Numeric value between 0 and 1 specifying the R2 cutoff
#'   for curve quality assessment. Curves with R2 below this value will be flagged
#'   as "Low R2" in the quality assessment (default: 0.8).
#' @param assay_type Character string specifying the assay type, used to set
#'   biologically plausible parameter limits in \code{check_biological_plausibility}.
#'   \code{"nanobret"} (default) uses limits appropriate for BRET ratio data
#'   (Bottom: -100 to 600, Top: 0 to 700), applied regardless of
#'   \code{normalize}.
#'   \code{"viability"} with \code{normalize = TRUE} uses limits appropriate
#'   for a normalized 0-100\% scale (Bottom: -20 to 60, Top: 50 to 130),
#'   allowing a small tolerance for over/undershoot.
#'   \code{"viability"} with \code{normalize = FALSE} disables the check
#'   entirely -- raw count data has no meaningful absolute limits since the
#'   scale is instrument-dependent.
#'   Any other value also disables the check entirely.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{summary_table}: Data frame with fitted parameters for all compounds
#'   \item \code{detailed_results}: List containing detailed fitting results for each compound
#'   \item \code{n_compounds}: Number of compounds analyzed
#'   \item \code{successful_fits}: Number of successful curve fits
#'   \item \code{normalized}: Logical indicating if data was normalized
#'   \item \code{original_data}: Original input data
#'   \item \code{processed_data}: Processed data (normalized if requested)
#'   \item \code{threshold_settings}: List with threshold settings and affected compounds (if applicable)
#' }
#'
#' @details
#' This function performs comprehensive dose-response analysis using the following steps:
#'
#' \strong{Data Structure:}
#' - First column: log10(inhibitor concentration) values
#' - Subsequent columns: Pairs of response measurements (technical duplicates)
#' - Example: Column1 = log_conc, Column2 = Response1, Column3 = Response2, Column4 = Response3, Column5 = Response4, etc.
#'
#' \strong{Model Fitting:}
#' - Uses 3-parameter logistic model with fixed Hill slope of -1
#' - Implements multiple starting value strategies for robust convergence
#' - Includes biological plausibility checks and automatic parameter corrections
#' - Calculates ideal Hill slope using 4-parameter model for comparison
#'
#' \strong{Output Parameters:}
#' - Bottom: Minimum response asymptote
#' - Top: Maximum response asymptote
#' - LogIC50: log10(IC50) value
#' - IC50: Half-maximal inhibitory concentration
#' - Span: Top - Bottom (response range)
#' - Confidence intervals for all parameters
#' - Goodness-of-fit metrics (R2, Syx, etc.)
#' - Curve quality assessment
#' - Maximum slope and ideal Hill slope
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_3pl")
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
#' mrt  <- ratio_res$plate_01$result$modified_ratio_table
#' fit3 <- fit_drc_3pl(data = mrt, normalize = FALSE, verbose = FALSE)
#' head(fit3$summary_table[, c("Compound", "LogIC50", "R_squared")])
#' }
#' @section Data Format:
#' The input data should be structured as follows:
#' \preformatted{
#'   log_conc  Targ/Comp1_rep1  Targ/Comp1_rep2  Targ/Comp2_rep1  Targ/Comp2_rep2 ...
#'     NA      25.1              25.1            25.1            25.1
#'   -8.0      30.5              31.2            22.1            23.4
#'   -7.0      45.2              46.8            35.6            36.9
#'   ...       ...            ...             ...             ...
#' }
#'
#' @section Curve Quality Assessment:
#' Curves are automatically assessed for quality based on:
#' \itemize{
#'   \item Maximum slope (shallow curves flagged)
#'   \item Response span (<20\% flagged)
#'   \item R-squared values (<0.5 flagged)
#'   \item Biological plausibility of parameters
#'   \item logIC50 range (>0.666 flagged)
#' }
#'
#' @export
#' @seealso
#' Useful related functions:
#' \itemize{
#'   \item [plot_dose_response()] for plotting results
#'   \item plot_multiple_compounds() for exporting results
#'   \item save_multiple_sheets() for saving results to Excel
#'   \item \code{\link[stats]{nls}} for nonlinear regression details
#' }
#'
#' @export
#'
#' @references
#' For methodological details on dose-response analysis:
#' \itemize{
#'   \item GraphPad Prism Curve Fitting Guide
#'   \item NIH Assay Guidance Manual
#' }
#'
#' @importFrom stats median predict sd setNames
#' @importFrom dplyr %>% mutate across
#' @importFrom stats lm coef as.formula
#' @importFrom utils head tail txtProgressBar setTxtProgressBar
#' @importFrom grDevices colorRampPalette rainbow heat.colors terrain.colors topo.colors cm.colors
#' @export
fit_drc_3pl <- function(data, output_file = NULL, normalize = FALSE, verbose = TRUE,
                        enforce_bottom_threshold = FALSE, bottom_threshold = 60,
                        r_sqr_threshold    = 0.8,
                        assay_type         = "nanobret") {
  # ============================================================================
  # 1. DEPENDENCY CHECK
  # ============================================================================
  required_packages <- c("dplyr", "stats", "grDevices")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages are not installed: ", paste(missing_packages, collapse = ", "),
         "\nPlease install using: install.packages(c(",
         paste0("\"", missing_packages, "\"", collapse = ", "), "))")
  }
  
  # Constants & Helpers
  PARAM_NAMES <- c("Bottom", "Top", "LogIC50", "IC50", "Span")
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  
  # ============================================================================
  # 2. MODEL DEFINITIONS
  # ============================================================================
  
  model_funcs <- list(
    inhibition = function(log_inhibitor, Bottom, Top, LogIC50) {
      Bottom + (Top - Bottom) / (1 + 10^(log_inhibitor - LogIC50))
    },
    activation = function(log_inhibitor, Bottom, Top, LogIC50) {
      Bottom + (Top - Bottom) / (1 + 10^(LogIC50 - log_inhibitor))
    }
  )
  
  four_param_model <- function(log_inhibitor, Bottom, Top, LogIC50, HillSlope) {
    Bottom + (Top - Bottom) / (1 + 10^((log_inhibitor - LogIC50) * HillSlope))
  }
  
  # ============================================================================
  # 3. HELPER FUNCTIONS
  # ============================================================================
  
  detect_curve_type <- function(df_clean) {
    if (nrow(df_clean) < 4) return("unknown")
    sorted_df <- df_clean[order(df_clean$log_inhibitor), ]
    resps <- sorted_df$response
    
    initial_avg <- mean(head(resps, 3), na.rm = TRUE)
    final_avg   <- mean(tail(resps, 3), na.rm = TRUE)
    
    # Scale-robust classifier: a curve is non-flat only if the difference
    # between the means of the first and last 3 responses (along the
    # log-concentration axis) exceeds BOTH:
    #   (a) 15% of the observed response range (catches strong curves at
    #       any scale, including raw BRET ratios with range ~10 units), and
    #   (b) k * SE(diff) where SE(diff) = sqrt((var(head3) + var(tail3))/3),
    #       k = 3 (~3-sigma signal-to-noise test, suppresses random-noise
    #       false positives on truly flat curves).
    # Threshold = max((a), (b)). Both relative and noise-aware -> works on
    # normalized 0-100% data AND raw plate-reader counts without an
    # absolute floor.
    range_val <- diff(range(resps, na.rm = TRUE))
    head3 <- head(resps, 3); tail3 <- tail(resps, 3)
    var_head <- if (length(head3) >= 2) stats::var(head3, na.rm = TRUE) else 0
    var_tail <- if (length(tail3) >= 2) stats::var(tail3, na.rm = TRUE) else 0
    if (is.na(var_head)) var_head <- 0
    if (is.na(var_tail)) var_tail <- 0
    se_diff <- sqrt((var_head + var_tail) / 3)
    thr_rel   <- range_val * 0.15
    thr_noise <- if (se_diff > 0) 3 * se_diff else 0
    threshold <- max(thr_rel, thr_noise)
    
    if (initial_avg > final_avg + threshold) return("inhibition")
    if (final_avg > initial_avg + threshold) return("activation")
    return("flat")
  }
  
  correct_hill_slope <- function(hill_slope, df_clean) {
    if (is.na(hill_slope)) return(NA_real_)
    ctype <- detect_curve_type(df_clean)
    if (ctype == "inhibition" && hill_slope > 0) return(-abs(hill_slope))
    if (ctype == "activation" && hill_slope < 0) return(abs(hill_slope))
    return(hill_slope)
  }
  
  create_empty_result <- function(comp_name, reason = "Model failed") {
    list(
      parameters = data.frame(Parameter = PARAM_NAMES, Value = rep(NA_real_, 5), stringsAsFactors = FALSE),
      confidence_intervals = list(
        Bottom = c(NA_real_, NA_real_), Top = c(NA_real_, NA_real_),
        LogIC50 = c(NA_real_, NA_real_), IC50 = c(NA_real_, NA_real_),
        Bottom_Lower = NA_real_, Bottom_Upper = NA_real_,
        Top_Lower = NA_real_, Top_Upper = NA_real_
      ),
      goodness_of_fit = list(
        R_squared = NA_real_, Syx = NA_real_, Sum_of_Squares = NA_real_,
        Degrees_of_Freedom = NA_real_
      ),
      curve_quality = paste("Fit", tolower(reason)),
      max_slope = NA_real_, ideal_hill_slope = NA_real_, model = NULL,
      success = FALSE, compound = comp_name
    )
  }
  
  check_biological_plausibility <- function(params, df_clean) {
    exp_min <- min(df_clean$response, na.rm = TRUE)
    exp_max <- max(df_clean$response, na.rm = TRUE)
    ctype   <- detect_curve_type(df_clean)
    
    # Parameter limits depend on assay type AND whether data are normalized.
    #
    # "nanobret"            : BRET ratio scale; wide limits to accommodate the
    #                         full dynamic range of NanoBRET assays. Applied
    #                         regardless of normalize, since NanoBRET data are
    #                         always on a ratio scale.
    # "viability" + normalized=TRUE  : 0-100% scale with small over/undershoot
    #                         tolerance (Bottom: -20 to 60, Top: 50 to 130).
    # "viability" + normalized=FALSE : raw counts -- no absolute limits are
    #                         meaningful (scale is instrument-dependent), so
    #                         the check is disabled entirely.
    # anything else         : check disabled.
    if (assay_type == "viability" && !normalize)
      return(list(needs_correction = FALSE))
    
    if (!assay_type %in% c("nanobret", "viability"))
      return(list(needs_correction = FALSE))
    
    if (assay_type == "nanobret") {
      bottom_limits  <- if (ctype == "activation") c(-100, Inf) else c(-100, 600)
      top_limits     <- c(0, 700)
    } else {
      # viability + normalize=TRUE: normalized scale, small over/undershoot tolerance
      bottom_limits  <- c(-20, 60)
      top_limits     <- c(50, 130)
    }
    logIC50_limits <- c(-20, 5)
    
    corrections <- list()
    reasons     <- list()
    
    if (params[1] < bottom_limits[1] || params[1] > bottom_limits[2] || !is.finite(params[1])) {
      corrections$Bottom <- max(bottom_limits[1], min(bottom_limits[2], exp_min))
      reasons$Bottom <- sprintf("Biologically implausible (%.2f). Using min: %.2f", params[1], exp_min)
    }
    if (params[2] < top_limits[1] || params[2] > top_limits[2] || !is.finite(params[2])) {
      corrections$Top <- max(top_limits[1], min(top_limits[2], exp_max))
      reasons$Top <- sprintf("Biologically implausible (%.2f). Using max: %.2f", params[2], exp_max)
    }
    if (params[3] < logIC50_limits[1] || params[3] > logIC50_limits[2] || !is.finite(params[3])) {
      fallback <- median(df_clean$log_inhibitor, na.rm = TRUE)
      corrections$LogIC50 <- max(logIC50_limits[1], min(logIC50_limits[2], fallback))
      reasons$LogIC50 <- sprintf("Biologically implausible (%.2f). Using median: %.2f", params[3], fallback)
    }
    
    if (length(corrections) > 0) {
      new_p <- params
      if ("Bottom"  %in% names(corrections)) new_p[1] <- corrections$Bottom
      if ("Top"     %in% names(corrections)) new_p[2] <- corrections$Top
      if ("LogIC50" %in% names(corrections)) new_p[3] <- corrections$LogIC50
      
      return(list(
        corrected_params    = c(new_p[1:3], 10^new_p[3], new_p[2] - new_p[1]),
        corrections_applied = corrections,
        correction_reasons  = reasons,
        needs_correction    = TRUE
      ))
    }
    list(needs_correction = FALSE)
  }
  
  normalize_dataframe <- function(df) {
    df |> dplyr::mutate(dplyr::across(-1, ~ {
      v <- suppressWarnings(as.numeric(as.character(.x)))
      vc <- v[!is.na(v)]
      if (length(vc) < 2) return(rep(NA_real_, length(v)))
      first <- vc[1]; last <- vc[length(vc)]
      if (!is.finite(first) || !is.finite(last) || abs(last - first) < .Machine$double.eps) {
        return(rep(NA_real_, length(v)))
      }
      (v - first) / (last - first) * 100
    }))
  }
  
  # ============================================================================
  # 4. SINGLE PAIR ANALYSIS LOGIC
  # ============================================================================
  
  analyze_single_pair <- function(pair_data, comp_name) {
    # 4.1 Data Preparation
    n_rep_cols <- ncol(pair_data) - 1L
    df_clean <- stats::na.omit(data.frame(
      log_inhibitor = rep(pair_data[, 1], n_rep_cols),
      response = as.numeric(unlist(pair_data[, -1]))
    ))
    
    if (nrow(df_clean) < 4 || stats::sd(df_clean$response, na.rm = TRUE) < 1e-6) {
      if (verbose) warning("Data validation failed for ", comp_name)
      return(create_empty_result(comp_name, "Validation failed"))
    }
    
    # Initial values
    min_resp <- min(df_clean$response, na.rm = TRUE)
    max_resp <- max(df_clean$response, na.rm = TRUE)
    approx_ic50 <- tryCatch({
      stats::approx(df_clean$response, df_clean$log_inhibitor,
                    xout = (min_resp + max_resp)/2, rule = 2)$y
    }, error = function(e) NULL) %||% stats::median(df_clean$log_inhibitor, na.rm = TRUE)
    
    # 4.2 Model Selection
    curve_type <- detect_curve_type(df_clean)
    model_func <- if (curve_type == "activation") model_funcs$activation else model_funcs$inhibition
    
    # 4.3 Fit Attempt Loop
    control_defs <- stats::nls.control(maxiter = 500, tol = 1e-04, minFactor = 1/4096, warnOnly = TRUE)
    control_relax <- stats::nls.control(maxiter = 1000, tol = 1e-03, minFactor = 1/1024, warnOnly = TRUE)
    
    start_strategies <- list(
      list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50),
      list(Bottom = min_resp * 0.8, Top = max_resp * 1.2, LogIC50 = approx_ic50),
      list(Bottom = min_resp * 1.2, Top = max_resp * 0.8, LogIC50 = approx_ic50),
      list(Bottom = max(0, min_resp - 10), Top = min(150, max_resp + 10), LogIC50 = approx_ic50)
    )
    
    fit <- NULL
    for (start_vals in start_strategies) {
      fit <- tryCatch({
        stats::nls(response ~ model_func(log_inhibitor, Bottom, Top, LogIC50),
                   data = df_clean, start = start_vals, control = control_defs, algorithm = "port")
      }, error = function(e) NULL)
      if (!is.null(fit)) break
    }
    if (is.null(fit)) {
      fit <- tryCatch({
        stats::nls(response ~ model_func(log_inhibitor, Bottom, Top, LogIC50),
                   data = df_clean, start = start_strategies[[1]], control = control_relax, algorithm = "port")
      }, error = function(e) NULL)
    }
    
    # 4.4 Handle Failure
    if (is.null(fit)) {
      if (verbose) warning("Could not fit model for ", comp_name)
      res <- create_empty_result(comp_name, "Approximate (model failed)")
      res$parameters$Value <- c(min_resp, max_resp, approx_ic50,
                                if(!is.na(approx_ic50)) 10^approx_ic50 else NA_real_,
                                max_resp - min_resp)
      res$success <- TRUE
      res$curve_type <- curve_type
      return(res)
    }
    
    # 4.5 Process Success
    params <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(params) || any(is.na(params))) {
      if (verbose) warning("Error extracting coefficients for ", comp_name)
      return(create_empty_result(comp_name, "Coefficient extraction failed"))
    }
    
    initial_params <- c(params, 10^params[3], params[2] - params[1])
    
    # Goodness of Fit
    gof_results <- tryCatch({
      resids <- stats::resid(fit)
      tss <- sum((df_clean$response - mean(df_clean$response, na.rm = TRUE))^2, na.rm = TRUE)
      rss <- sum(resids^2, na.rm = TRUE)
      dof <- max(1, nrow(df_clean) - length(params))
      list(
        R_squared = if(tss > 1e-10) 1 - rss/tss else 0,
        Syx = if(dof > 0) sqrt(rss/dof) else NA_real_,
        Sum_of_Squares = rss,
        Degrees_of_Freedom = dof
      )
    }, error = function(e) {
      list(R_squared = NA_real_, Syx = NA_real_, Sum_of_Squares = NA_real_, Degrees_of_Freedom = NA_real_)
    })
    
    # Confidence Intervals
    ci_results <- tryCatch({
      # Try Profile Likelihood first
      tryCatch({
        ci_b <- unname(stats::confint(stats::profile(fit, "Bottom"), level = 0.95))
        ci_t <- unname(stats::confint(stats::profile(fit, "Top"), level = 0.95))
        ci_l <- unname(stats::confint(stats::profile(fit, "LogIC50"), level = 0.95))
        
        # Safe IC50 calculation from log values
        ic50_ci <- if (!any(is.na(ci_l))) 10^ci_l else c(NA_real_, NA_real_)
        
        list(
          Bottom = ci_b, Top = ci_t, LogIC50 = ci_l, IC50 = ic50_ci,
          Bottom_Lower = ci_b[1], Bottom_Upper = ci_b[2],
          Top_Lower = ci_t[1], Top_Upper = ci_t[2]
        )
      }, error = function(e) {
        # Fallback to Asymptotic (Standard Errors)
        fit_sum <- tryCatch(summary(fit), error = function(e2) NULL)
        
        if (is.null(fit_sum)) {
          return(list(
            Bottom = c(NA_real_, NA_real_), Top = c(NA_real_, NA_real_),
            LogIC50 = c(NA_real_, NA_real_), IC50 = c(NA_real_, NA_real_),
            Bottom_Lower = NA_real_, Bottom_Upper = NA_real_,
            Top_Lower = NA_real_, Top_Upper = NA_real_
          ))
        }
        
        se_mat <- fit_sum$coefficients
        
        get_se <- function(pname) {
          if (!is.null(se_mat) && pname %in% rownames(se_mat)) {
            se_mat[pname, "Std. Error"]
          } else {
            NA_real_
          }
        }
        
        z <- 1.96
        
        calc_int <- function(val, se_val) {
          if (is.na(val) || is.na(se_val)) return(c(NA_real_, NA_real_))
          c(val - z * se_val, val + z * se_val)
        }
        
        ci_b <- calc_int(params["Bottom"], get_se("Bottom"))
        ci_t <- calc_int(params["Top"], get_se("Top"))
        ci_l <- calc_int(params["LogIC50"], get_se("LogIC50"))
        
        ic50_ci <- if (!any(is.na(ci_l))) 10^ci_l else c(NA_real_, NA_real_)
        
        list(
          Bottom = ci_b, Top = ci_t, LogIC50 = ci_l, IC50 = ic50_ci,
          Bottom_Lower = ci_b[1], Bottom_Upper = ci_b[2],
          Top_Lower = ci_t[1], Top_Upper = ci_t[2]
        )
      })
    }, error = function(e) {
      list(
        Bottom = c(NA_real_, NA_real_), Top = c(NA_real_, NA_real_),
        LogIC50 = c(NA_real_, NA_real_), IC50 = c(NA_real_, NA_real_),
        Bottom_Lower = NA_real_, Bottom_Upper = NA_real_,
        Top_Lower = NA_real_, Top_Upper = NA_real_
      )
    })
    
    # Plausibility Check
    check <- check_biological_plausibility(params, df_clean)
    final_params <- if (check$needs_correction) check$corrected_params else initial_params
    
    # Quality & Hill Slope
    ideal_hs <- tryCatch({
      start_4pl <- list(
        Bottom = final_params[1], Top = final_params[2],
        LogIC50 = final_params[3],
        HillSlope = if(curve_type == "activation") 1 else -1
      )
      fit_4 <- stats::nls(
        response ~ four_param_model(log_inhibitor, Bottom, Top, LogIC50, HillSlope),
        data = df_clean, start = start_4pl,
        control = stats::nls.control(maxiter = 200, warnOnly = TRUE),
        algorithm = "port"
      )
      correct_hill_slope(round(stats::coef(fit_4)["HillSlope"], 3), df_clean)
    }, error = function(e) NA_real_)
    
    # Curve Quality Metrics
    span_val <- final_params[5]
    max_slope <- if (!is.na(span_val)) -span_val * log(10) / 4 else NA_real_
    q_flags <- character()
    if (!is.na(max_slope) && abs(max_slope) < 5) q_flags <- c(q_flags, "Very shallow slope")
    else if (!is.na(max_slope) && abs(max_slope) < 15) q_flags <- c(q_flags, "Shallow slope")
    if (!is.na(span_val) && abs(span_val) < 20) q_flags <- c(q_flags, "Small span")
    if (!is.na(gof_results$R_squared) && gof_results$R_squared < r_sqr_threshold) {
      q_flags <- c(q_flags, "Low R2")
    }
    if (check$needs_correction) q_flags <- c(q_flags, "Params corrected")
    
    # Return complete result
    list(
      parameters = data.frame(Parameter = PARAM_NAMES, Value = final_params, stringsAsFactors = FALSE),
      confidence_intervals = ci_results,
      goodness_of_fit = gof_results,
      curve_quality = if(length(q_flags) == 0) "Good curve" else paste(q_flags, collapse = "; "),
      max_slope = max_slope,
      ideal_hill_slope = ideal_hs,
      model = fit,
      success = TRUE,
      compound = comp_name,
      biological_plausibility_check = check,
      data = df_clean,
      curve_type = curve_type
    )
  }
  
  # ============================================================================
  # 5. DATA SETUP
  # ============================================================================
  
  if (ncol(data) < 3) stop("Data must have at least 3 columns")
  
  original_data <- data
  
  if (verbose) message("Generating normalized data (0-100%)...")
  normalized_data <- normalize_dataframe(data)
  
  if (normalize) {
    if (verbose) message("Analysis will use: NORMALIZED data.")
    analysis_data <- normalized_data
  } else {
    if (verbose) message("Analysis will use: ORIGINAL data.")
    analysis_data <- original_data
  }
  
  # ============================================================================
  # 6. MAIN ANALYSIS LOOP
  # ============================================================================
  
  if (verbose) message("Starting dose-response analysis...")
  
  # Group data columns by base compound name (strip trailing .2, .3, ... suffix).
  # This supports any number of replicates: 2 replicates (original behaviour),
  # 4 replicates (two merged plates), 6 replicates (three merged plates), etc.
  c_names   <- colnames(analysis_data)[-1]          # data column names only
  base_names <- sub("\\.\\d+$", "", c_names)         # strip .2 / .3 / .4 ...
  compound_order <- unique(base_names)               # preserve first-appearance order
  n_compounds <- length(compound_order)
  
  all_results <- lapply(seq_len(n_compounds), function(i) {
    base  <- compound_order[i]
    # All column indices (1-based in analysis_data) belonging to this compound
    data_cols <- which(base_names == base) + 1L      # +1 for the conc column offset
    
    comp_name <- base
    
    if (verbose) cat(sprintf("\rProcessing compound %d/%d: %s", i, n_compounds, substr(comp_name, 1, 30)))
    analyze_single_pair(analysis_data[, c(1L, data_cols), drop = FALSE], comp_name)
  })
  if (verbose) cat("\n")
  
  # ============================================================================
  # 7. SUMMARY GENERATION
  # ============================================================================
  
  # Formatting functions
  fmt_sci <- function(x) {
    if (is.null(x) || is.na(x)) return(NA_character_)
    format(x, scientific = TRUE)
  }
  
  fmt_rd <- function(x, d = 3) {
    if (is.null(x) || is.na(x)) return(NA_real_)
    round(x, d)
  }
  
  summary_list <- lapply(all_results, function(res) {
    # Safe compound name extraction
    safe_compound_name <- tryCatch({
      parts <- strsplit(res$compound, " \\| ")[[1]]
      if (length(parts) > 0) parts[1] else res$compound
    }, error = function(e) res$compound)
    
    # Empty row for failed fits
    empty_row <- data.frame(
      Compound = safe_compound_name,
      Bottom = NA_real_, Top = NA_real_, LogIC50 = NA_real_, IC50 = NA_character_,
      Bottom_Lower_95CI = NA_real_, Bottom_Upper_95CI = NA_real_,
      Top_Lower_95CI = NA_real_, Top_Upper_95CI = NA_real_,
      LogIC50_Lower_95CI = NA_real_, LogIC50_Upper_95CI = NA_real_,
      IC50_Lower_95CI = NA_character_, IC50_Upper_95CI = NA_character_,
      Span = NA_real_,
      R_squared = NA_real_, Syx = NA_real_, Sum_of_Squares = NA_real_,
      Degrees_of_Freedom = NA_real_,
      Max_Slope = NA_real_, Ideal_Hill_Slope = NA_real_,
      Curve_Quality = res$curve_quality %||% "Fit failed",
      stringsAsFactors = FALSE
    )
    
    if (!isTRUE(res$success)) return(empty_row)
    p <- res$parameters$Value
    if (all(is.na(p))) return(empty_row)
    
    # Threshold Logic
    apply_thresh <- FALSE
    if (enforce_bottom_threshold && !is.na(p[1]) && p[1] >= bottom_threshold) {
      if ((res$curve_type %||% "unknown") %in% c("inhibition", "flat")) apply_thresh <- TRUE
    }
    
    ci <- res$confidence_intervals
    gof <- res$goodness_of_fit
    
    # Build summary row
    data.frame(
      Compound = safe_compound_name,
      Bottom = fmt_rd(p[1]),
      Top = fmt_rd(p[2]),
      LogIC50 = if (!apply_thresh) fmt_rd(p[3]) else NA_real_,
      IC50 = if (!apply_thresh) fmt_sci(p[4]) else NA_character_,
      
      Bottom_Lower_95CI = fmt_rd(ci$Bottom_Lower),
      Bottom_Upper_95CI = fmt_rd(ci$Bottom_Upper),
      Top_Lower_95CI = fmt_rd(ci$Top_Lower),
      Top_Upper_95CI = fmt_rd(ci$Top_Upper),
      
      LogIC50_Lower_95CI = if (!apply_thresh) fmt_rd(ci$LogIC50[1]) else NA_real_,
      LogIC50_Upper_95CI = if (!apply_thresh) fmt_rd(ci$LogIC50[2]) else NA_real_,
      IC50_Lower_95CI = if (!apply_thresh) fmt_sci(ci$IC50[1]) else NA_character_,
      IC50_Upper_95CI = if (!apply_thresh) fmt_sci(ci$IC50[2]) else NA_character_,
      
      Span = fmt_rd(p[5]),
      R_squared = fmt_rd(gof$R_squared),
      Syx = fmt_rd(gof$Syx),
      Sum_of_Squares = fmt_rd(gof$Sum_of_Squares),
      Degrees_of_Freedom = gof$Degrees_of_Freedom,
      Max_Slope = fmt_rd(res$max_slope),
      Ideal_Hill_Slope = fmt_rd(res$ideal_hill_slope),
      Curve_Quality = res$curve_quality,
      stringsAsFactors = FALSE
    )
  })
  
  # Consolidated Summary Table
  summary_table <- dplyr::bind_rows(summary_list)
  
  # Final Summary (Transposed)
  if (nrow(summary_table) == 0) {
    final_summary_table <- data.frame()
  } else {
    t_data <- as.data.frame(t(summary_table[, -1]))
    colnames(t_data) <- summary_table$Compound
    final_summary_table <- t_data
  }
  
  # Identify compounds affected by threshold
  threshold_affected <- character()
  if (enforce_bottom_threshold) {
    for (result in all_results) {
      if (isTRUE(result$success) &&
          !is.na(result$parameters$Value[1]) &&
          result$parameters$Value[1] >= bottom_threshold &&
          result$curve_type == "inhibition") {
        comp_name_parts <- strsplit(result$compound, " \\| ")[[1]]
        comp_name <- if (length(comp_name_parts) > 0) comp_name_parts[1] else result$compound
        threshold_affected <- c(threshold_affected, comp_name)
      }
    }
  }
  
  # ============================================================================
  # 8. REPORTING & EXPORT
  # ============================================================================
  
  if (verbose && nrow(summary_table) > 0) {
    succ_n <- sum(!is.na(summary_table$IC50))
    tot_n  <- nrow(summary_table)
    cat("\n", strrep("=", 50), "\n", sep = "")
    cat("DOSE-RESPONSE ANALYSIS COMPLETE\n")
    cat(strrep("=", 50), "\n")
    cat("  . Compounds analyzed: ", tot_n, "\n")
    cat("  . Successful fits: ", succ_n, " (", round(succ_n/tot_n*100, 1), "%)\n", sep = "")
    cat("  . Data used for fit: ", if(normalize) "Normalized (0-100%)" else "Original (Raw)", "\n")
    
    if (enforce_bottom_threshold && length(threshold_affected) > 0) {
      excl <- sum(is.na(summary_table$IC50) & !is.na(summary_table$Bottom))
      cat("  . IC50 values excluded (Bottom >= ", bottom_threshold, "): ", excl, "\n", sep = "")
      
      if (verbose > 1) {
        cat("\nCOMPOUNDS WITH EXCLUDED IC50 VALUES:\n")
        for (i in seq_along(threshold_affected)) {
          bottom_val <- summary_table$Bottom[summary_table$Compound == threshold_affected[i]]
          cat("  . ", threshold_affected[i], " (Bottom = ", round(bottom_val, 1), ")\n", sep = "")
        }
      }
    }
    
    if ("Curve_Quality" %in% names(summary_table)) {
      cat("\nCURVE QUALITY DISTRIBUTION:\n")
      quality_counts <- table(summary_table$Curve_Quality)
      for (quality in names(sort(quality_counts, decreasing = TRUE))) {
        count <- quality_counts[quality]
        cat("  . ", sprintf("%-35s: %d (%.1f%%)", quality, count, count/tot_n*100), "\n")
      }
    }
    cat(strrep("=", 50), "\n")
  }
  
  # Save results to file
  if (!is.null(output_file)) {
    f_ext <- tolower(tools::file_ext(output_file))
    
    if (f_ext == "xlsx" && requireNamespace("openxlsx", quietly = TRUE)) {
      wb <- openxlsx::createWorkbook()
      
      # Sheet 1: Final_Summary (Transposed)
      openxlsx::addWorksheet(wb, "Final_Summary")
      if (nrow(final_summary_table) > 0) {
        out_final <- data.frame(
          Parameter = rownames(final_summary_table),
          final_summary_table,
          check.names = FALSE
        )
        openxlsx::writeData(wb, "Final_Summary", out_final)
      } else {
        openxlsx::writeData(wb, "Final_Summary", "No data")
      }
      
      # Sheet 2: Summary (Detailed)
      openxlsx::addWorksheet(wb, "Summary")
      openxlsx::writeData(wb, "Summary", summary_table)
      
      # Sheet 3: Normalized_Data
      openxlsx::addWorksheet(wb, "Normalized_Data")
      openxlsx::writeData(wb, "Normalized_Data", normalized_data)
      
      # Sheet 4: Original_Data
      openxlsx::addWorksheet(wb, "Original_Data")
      openxlsx::writeData(wb, "Original_Data", original_data)
      
      openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
      
      if(verbose) {
        message("Results saved to Excel: ", output_file)
        message("  1. Final_Summary (Transposed)")
        message("  2. Summary (Detailed)")
        message("  3. Normalized_Data (0-100%)")
        message("  4. Original_Data (Raw)")
      }
      
    } else {
      # Fallback to CSV
      if (f_ext == "xlsx") {
        output_file <- sub("\\.xlsx$", ".csv", output_file, ignore.case = TRUE)
        if(verbose) warning("openxlsx package missing. Falling back to CSV format.")
      }
      utils::write.csv(summary_table, output_file, row.names = FALSE)
      if(verbose) message("Saved Summary CSV to: ", output_file)
    }
  }
  
  # ============================================================================
  # 9. RETURN RESULTS
  # ============================================================================
  
  list(
    summary_table = summary_table,
    final_summary_table = final_summary_table,
    detailed_results = all_results,
    original_data = original_data,
    normalized_data = normalized_data,
    used_normalized_data = normalize,
    n_compounds = n_compounds,
    successful_fits = sum(!is.na(summary_table$IC50)),
    normalized = normalize,
    threshold_settings = if (enforce_bottom_threshold) {
      list(
        enforce = TRUE,
        threshold = bottom_threshold,
        affected_compounds = threshold_affected
      )
    } else NULL
  )
}
