#' Dose-Response Curve Analysis using 4-Parameter Logistic Model
#'
#' Fits 4-parameter logistic (4PL) models to dose-response data for both inhibition
#' and activation curves. Automatically detects curve type, applies biological
#' plausibility checks, and provides comprehensive curve quality assessment.
#'
#' @param data A data.frame where the first column contains log inhibitor concentrations
#' and subsequent columns contain response values in duplicate pairs.
#' @param output_file Optional character string specifying path to save results.
#' Supports .csv and .xlsx formats.
#' @param normalize Logical indicating whether to normalize response data to 0-100%
#' based on first and last data points. Default is FALSE.
#' @param verbose Logical indicating whether to display detailed progress messages.
#' Default is FALSE.
#' @param enforce_bottom_threshold Logical indicating whether to exclude IC50 values
#' for inhibition curves where Bottom parameter exceeds threshold. Default is FALSE.
#' @param bottom_threshold Numeric value (0-100) for Bottom parameter above which
#' IC50 values are excluded for inhibition curves. Default is 60.
#' @param r_sqr_threshold Numeric value (0-1) for minimum R-squared to consider
#' curve fit acceptable. Default is 0.5.
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
#' @return A list containing:
#' \itemize{
#' \item \code{summary_table}: Data.frame with fitted parameters, confidence intervals,
#' and quality metrics for all compounds
#' \item \code{detailed_results}: List of detailed fitting results for each compound
#' \item \code{n_compounds}: Number of compounds analyzed
#' \item \code{successful_fits}: Number of successful model fits
#' \item \code{normalized}: Logical indicating if data was normalized
#' \item \code{original_data}: Original input data
#' \item \code{processed_data}: Processed data used for analysis
#' \item \code{threshold_settings}: Settings for bottom threshold enforcement
#' \item \code{parameter_order_corrections}: Number of parameter order corrections applied
#' }
#'
#' @details
#' The function uses a 4-parameter logistic model:
#' \deqn{Response = Bottom + \frac{Top - Bottom}{1 + 10^{(log_{10}(inhibitor) - LogIC50) \times HillSlope}}}
#'
#' Key features:
#' \itemize{
#' \item Automatic detection of inhibition (HillSlope < 0) vs activation (HillSlope > 0) curves
#' \item Biological plausibility checks with parameter corrections when needed
#' \item Multiple starting value strategies for robust convergence
#' \item Comprehensive curve quality assessment based on multiple metrics
#' \item Confidence interval calculation via profiling or normal approximation
#' }
#'
#' @section Data Format:
#' Input data should be structured as:
#' \preformatted{
#' log_inhibitor | Compound1_rep1 | Compound1_rep2 | Compound2_rep1 | Compound2_rep2 | ...
#' -3.0 | 45.2 | 47.8 | 32.1 | 30.9 | ...
#' -2.0 | 38.7 | 40.1 | 28.5 | 29.2 | ...
#' -1.0 | 25.4 | 23.9 | 45.6 | 47.2 | ...
#' 0.0 | 15.2 | 16.8 | 62.3 | 60.9 | ...
#' }
#'
#'
#' @section Curve Quality Assessment:
#' Curves are classified based on:
#' \itemize{
#' \item R-squared < threshold: "Low R2"
#' \item Maximum slope < 5: "Very shallow slope"
#' \item Maximum slope < 15: "Shallow slope"
#' \item Span < 20: "Small span"
#' \item Parameter corrections: "Parameters corrected"
#' \item logIC50 range (>0.666 flagged)
#' }
#'
#' @examples
#' # Example with dummy data
#' data <- data.frame(
#'   log_inhibitor = c(-2, -1, 0, 1, 2),
#'   compound1_rep1 = c(100, 80, 50, 30, 10),
#'   compound1_rep2 = c(102, 82, 48, 28, 12)
#' )
#' result <- fit_drc_4pl(data, normalize = TRUE, verbose = TRUE)
#' head(result$summary_table)
#'
#' @export


fit_drc_4pl <- function(data, output_file = NULL, normalize = FALSE, verbose = TRUE,
                        enforce_bottom_threshold = FALSE, bottom_threshold = 60,
                        r_sqr_threshold = 0.8,
                        assay_type = "nanobret") {

  # Check and load required packages
  required_packages <- c("dplyr", "stats", "graphics", "grDevices")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Required packages are not installed: ", paste(missing_packages, collapse = ", "),
         "\nPlease install using: install.packages(c(",
         paste0("\"", missing_packages, "\"", collapse = ", "), "))")
  }

  library(dplyr, warn.conflicts = FALSE)

  # Constants for 4-parameter model
  PARAM_NAMES <- c("Bottom", "Top", "LogIC50", "HillSlope", "IC50", "Span")

  # NULL coalescing operator
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

  # 4-parameter logistic model
  four_param_model <- function(log_inhibitor, Bottom, Top, LogIC50, HillSlope) {
    Bottom + (Top - Bottom) / (1 + 10^((log_inhibitor - LogIC50) * HillSlope))
  }

  # Detect curve type based on data pattern.
  # Uses a relative threshold (max(15, range * 0.15)) on the head/tail
  # mean difference -- identical to fit_drc_3pl. This scales correctly
  # for both raw count data (large dynamic range) and normalized 0-100%
  # data, avoiding the fixed-threshold misclassification of flat curves.
  detect_curve_type <- function(data) {
    if (nrow(data) < 4) return("unknown")
    sorted_df <- data[order(data$log_inhibitor), ]
    resps <- sorted_df$response[!is.na(sorted_df$response)]
    if (length(resps) < 4) return("unknown")

    initial_avg <- mean(head(resps, 3), na.rm = TRUE)
    final_avg   <- mean(tail(resps, 3), na.rm = TRUE)

    range_val <- diff(range(resps, na.rm = TRUE))
    threshold <- max(15, range_val * 0.15)

    if (initial_avg > final_avg + threshold) return("inhibition")
    if (final_avg > initial_avg + threshold) return("activation")
    return("flat")
  }

  # Create empty result structure for failed fits
  create_empty_result <- function(comp_name, reason = "Model failed") {
    list(
      parameters = data.frame(Parameter = PARAM_NAMES, Value = rep(NA, 6), stringsAsFactors = FALSE),
      confidence_intervals = list(
        Bottom = c(NA, NA), Top = c(NA, NA), LogIC50 = c(NA, NA),
        HillSlope = c(NA, NA), IC50 = c(NA, NA),
        Bottom_Lower = NA, Bottom_Upper = NA, Top_Lower = NA, Top_Upper = NA
      ),
      goodness_of_fit = list(
        R_squared = NA, Syx = NA, Sum_of_Squares = NA,
        Total_Sum_of_Squares = NA, Regression_Sum_of_Squares = NA, Degrees_of_Freedom = NA
      ),
      curve_quality = paste("Fit", tolower(reason)),
      max_slope = NA,
      model = NULL,
      success = FALSE,
      compound = comp_name
    )
  }

  # Recalculate dependent parameters after corrections
  recalculate_dependent_params <- function(params, corrections) {
    corrected_params <- params

    if ("Bottom" %in% names(corrections)) corrected_params[1] <- corrections$Bottom
    if ("Top" %in% names(corrections)) corrected_params[2] <- corrections$Top
    if ("LogIC50" %in% names(corrections)) corrected_params[3] <- corrections$LogIC50
    if ("HillSlope" %in% names(corrections)) corrected_params[4] <- corrections$HillSlope

    c(corrected_params[1:4], 10^corrected_params[3], corrected_params[2] - corrected_params[1])
  }

  # Recalculate confidence intervals after parameter corrections
  recalculate_ci <- function(original_ci, corrections) {
    ci <- original_ci

    if ("Bottom" %in% names(corrections)) {
      ci$Bottom <- ci$Bottom_Lower <- ci$Bottom_Upper <- NA
    }
    if ("Top" %in% names(corrections)) {
      ci$Top <- ci$Top_Lower <- ci$Top_Upper <- NA
    }
    if ("LogIC50" %in% names(corrections)) {
      ci$LogIC50 <- ci$IC50 <- c(NA, NA)
    }
    if ("HillSlope" %in% names(corrections)) {
      ci$HillSlope <- c(NA, NA)
    }

    ci
  }

  # Correct parameter order based on curve type
  correct_parameter_order <- function(params, data, curve_type) {
    bottom <- params[1]
    top <- params[2]
    log_ic50 <- params[3]
    hill_slope <- params[4]

    # Check consistency between curve type and Hill Slope
    if (curve_type == "inhibition" && !is.na(hill_slope) && hill_slope > 0) {
      corrected_hill_slope <- -abs(hill_slope)

      if (bottom > top) {
        corrected_params <- c(top, bottom, log_ic50, corrected_hill_slope)
      } else {
        corrected_params <- c(bottom, top, log_ic50, corrected_hill_slope)
      }
      return(list(
        corrected_params = corrected_params,
        was_corrected = TRUE,
        correction_reason = "Hill Slope inconsistent with curve type (inhibition)"
      ))
    }

    if (curve_type == "activation" && !is.na(hill_slope) && hill_slope < 0) {
      corrected_hill_slope <- abs(hill_slope)

      if (bottom > top) {
        corrected_params <- c(top, bottom, log_ic50, corrected_hill_slope)
      } else {
        corrected_params <- c(bottom, top, log_ic50, corrected_hill_slope)
      }
      return(list(
        corrected_params = corrected_params,
        was_corrected = TRUE,
        correction_reason = "Hill Slope inconsistent with curve type (activation)"
      ))
    }

    # Check for Bottom/Top inversion only
    if (!is.na(bottom) && !is.na(top) && bottom > top) {
      corrected_params <- c(top, bottom, log_ic50, hill_slope)
      return(list(
        corrected_params = corrected_params,
        was_corrected = TRUE,
        correction_reason = "Bottom and Top inverted"
      ))
    }

    list(
      corrected_params = params,
      was_corrected = FALSE,
      correction_reason = NA
    )
  }

  # Check biological plausibility and apply corrections if needed.
  # Limits depend on assay_type and normalize, matching fit_drc_3pl behaviour:
  #   "nanobret"                  : BRET ratio scale; wide limits always applied.
  #   "viability" + normalize=TRUE: 0-100% scale with small over/undershoot tolerance.
  #   "viability" + normalize=FALSE: raw counts -- no absolute limits, check disabled.
  #   anything else               : check disabled.
  check_biological_plausibility <- function(params, data) {
    if (assay_type == "viability" && !normalize)
      return(list(needs_correction = FALSE))

    if (!assay_type %in% c("nanobret", "viability"))
      return(list(needs_correction = FALSE))

    responses <- data$response[!is.na(data$response)]
    exp_min <- min(responses, na.rm = TRUE)
    exp_max <- max(responses, na.rm = TRUE)

    curve_type <- detect_curve_type(data)

    if (assay_type == "nanobret") {
      bottom_limits <- if (curve_type == "activation") c(-100, Inf) else c(-100, 600)
      top_limits    <- c(0, 700)
    } else {
      # viability + normalize=TRUE
      bottom_limits <- c(-20, 60)
      top_limits    <- c(50, 130)
    }
    logIC50_limits    <- c(-20, 5)
    hill_slope_limits <- c(-5, 5)

    corrections <- list()
    reasons     <- list()

    # Check Bottom
    if (params[1] < bottom_limits[1] || params[1] > bottom_limits[2] || !is.finite(params[1])) {
      corrections$Bottom <- max(bottom_limits[1], min(bottom_limits[2], exp_min))
      reasons$Bottom <- sprintf("Biologically implausible (%.2f). Using experimental minimum: %.2f",
                                params[1], exp_min)
    }

    # Check Top
    if (params[2] < top_limits[1] || params[2] > top_limits[2] || !is.finite(params[2])) {
      corrections$Top <- max(top_limits[1], min(top_limits[2], exp_max))
      reasons$Top <- sprintf("Biologically implausible (%.2f). Using experimental maximum: %.2f",
                             params[2], exp_max)
    }

    # Check LogIC50
    if (params[3] < logIC50_limits[1] || params[3] > logIC50_limits[2] || !is.finite(params[3])) {
      fallback <- median(data$log_inhibitor, na.rm = TRUE)
      corrections$LogIC50 <- max(logIC50_limits[1], min(logIC50_limits[2], fallback))
      reasons$LogIC50 <- sprintf("Biologically implausible (%.2f). Using median concentration: %.2f",
                                 params[3], fallback)
    }

    # Check HillSlope
    if (params[4] < hill_slope_limits[1] || params[4] > hill_slope_limits[2] || !is.finite(params[4])) {
      default_hill <- if (curve_type == "activation") 1 else -1
      corrections$HillSlope <- max(hill_slope_limits[1], min(hill_slope_limits[2], default_hill))
      reasons$HillSlope <- sprintf("Biologically implausible (%.2f). Using default: %.2f",
                                   params[4], default_hill)
    }

    if (length(corrections) > 0) {
      return(list(
        corrected_params = recalculate_dependent_params(params, corrections),
        corrections_applied = corrections,
        correction_reasons = reasons,
        needs_correction = TRUE
      ))
    }

    list(needs_correction = FALSE)
  }

  # Prepare and validate data for analysis
  prepare_and_validate_data <- function(pair_data) {
    log_inhibitor <- pair_data[, 1]
    response <- as.numeric(unlist(pair_data[, -1]))
    n_rep_cols <- ncol(pair_data) - 1L

    df <- data.frame(log_inhibitor = rep(log_inhibitor, n_rep_cols), response = response)
    df_clean <- stats::na.omit(df)

    if (nrow(df_clean) < 5 || stats::sd(df_clean$response, na.rm = TRUE) < 1e-6) {
      return(list(valid = FALSE, df_clean = df_clean))
    }

    min_resp <- min(df_clean$response, na.rm = TRUE)
    max_resp <- max(df_clean$response, na.rm = TRUE)

    # Calculate approximate IC50
    approx_ic50 <- tryCatch({
      suppressWarnings(
        stats::approx(df_clean$response, df_clean$log_inhibitor,
                      xout = (min_resp + max_resp) / 2)$y
      )
    }, error = function(e) stats::median(df_clean$log_inhibitor, na.rm = TRUE))

    list(valid = TRUE, df_clean = df_clean,
         min_response = min_resp, max_response = max_resp, approx_ic50 = approx_ic50)
  }

  # Robust nonlinear fitting with multiple strategies
  try_robust_fit <- function(df_clean, min_resp, max_resp, approx_ic50, curve_type) {
    if (nrow(df_clean) < 5) return(NULL)

    control_configs <- list(
      default = stats::nls.control(maxiter = 500, tol = 1e-04, minFactor = 1/4096, warnOnly = TRUE),
      relaxed = stats::nls.control(maxiter = 1000, tol = 1e-03, minFactor = 1/1024, warnOnly = TRUE)
    )

    # Define start strategies based on curve type
    if (curve_type == "inhibition") {
      start_strategies <- list(
        list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50, HillSlope = -1),
        list(Bottom = min_resp * 0.8, Top = max_resp * 1.2, LogIC50 = approx_ic50, HillSlope = -1.5),
        list(Bottom = min_resp * 1.2, Top = max_resp * 0.8, LogIC50 = approx_ic50, HillSlope = -2),
        list(Bottom = max(0, min_resp - 10), Top = min(150, max_resp + 10), LogIC50 = approx_ic50, HillSlope = -0.8),
        list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50, HillSlope = -0.5)
      )
    } else if (curve_type == "activation") {
      start_strategies <- list(
        list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50, HillSlope = 1),
        list(Bottom = min_resp * 0.8, Top = max_resp * 1.2, LogIC50 = approx_ic50, HillSlope = 1.5),
        list(Bottom = min_resp * 1.2, Top = max_resp * 0.8, LogIC50 = approx_ic50, HillSlope = 2),
        list(Bottom = max(0, min_resp - 10), Top = min(150, max_resp + 10), LogIC50 = approx_ic50, HillSlope = 0.8),
        list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50, HillSlope = 0.5)
      )
    } else {
      # For flat or unknown curves: try both directions
      start_strategies <- list(
        list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50, HillSlope = -1),
        list(Bottom = min_resp, Top = max_resp, LogIC50 = approx_ic50, HillSlope = 1),
        list(Bottom = min_resp * 0.8, Top = max_resp * 1.2, LogIC50 = approx_ic50, HillSlope = -1.5),
        list(Bottom = min_resp * 0.8, Top = max_resp * 1.2, LogIC50 = approx_ic50, HillSlope = 1.5)
      )
    }

    # Try all strategies with default control
    for (start_vals in start_strategies) {
      fit <- tryCatch({
        stats::nls(
          response ~ four_param_model(log_inhibitor, Bottom, Top, LogIC50, HillSlope),
          data = df_clean, start = start_vals, control = control_configs$default, algorithm = "port"
        )
      }, error = function(e) NULL)

      if (!is.null(fit)) return(fit)
    }

    # Final attempt with relaxed parameters
    tryCatch({
      stats::nls(
        response ~ four_param_model(log_inhibitor, Bottom, Top, LogIC50, HillSlope),
        data = df_clean, start = start_strategies[[1]], control = control_configs$relaxed, algorithm = "port"
      )
    }, error = function(e) NULL)
  }

  # Calculate goodness of fit metrics
  calculate_goodness_of_fit <- function(fit, df_clean) {
    tryCatch({
      residuals <- stats::resid(fit)
      mean_resp <- mean(df_clean$response, na.rm = TRUE)
      total_ss <- sum((df_clean$response - mean_resp)^2)
      residual_ss <- sum(residuals^2)
      regression_ss <- total_ss - residual_ss
      n_obs <- nrow(df_clean)
      dof <- max(1, n_obs - length(stats::coef(fit)))

      list(
        R_squared = if (total_ss > 0) regression_ss / total_ss else 0,
        Syx = sqrt(residual_ss / dof),
        Sum_of_Squares = residual_ss,
        Degrees_of_Freedom = dof
      )
    }, error = function(e) {
      list(R_squared = NA, Syx = NA, Sum_of_Squares = NA, Degrees_of_Freedom = NA)
    })
  }

  # Calculate confidence intervals safely
  calculate_ci <- function(fit) {
    safe_ci <- function(param_name) {
      tryCatch({
        prof <- stats::profile(fit, which = param_name, alpha = 0.05)
        unname(stats::confint(prof, level = 0.95))
      }, error = function(e) {
        tryCatch({
          if (requireNamespace("lmtest", quietly = TRUE)) {
            unname(lmtest::confint2(fit, level = 0.95)[param_name, ])
          } else {
            se <- summary(fit)$coefficients[param_name, "Std. Error"]
            est <- stats::coef(fit)[param_name]
            z <- stats::qnorm(0.975)
            c(est - z * se, est + z * se)
          }
        }, error = function(e2) c(NA, NA))
      })
    }

    bottom_ci <- safe_ci("Bottom")
    top_ci <- safe_ci("Top")
    logIC50_ci <- safe_ci("LogIC50")
    hill_slope_ci <- safe_ci("HillSlope")

    list(
      Bottom = bottom_ci, Top = top_ci, LogIC50 = logIC50_ci, HillSlope = hill_slope_ci,
      IC50 = 10^logIC50_ci,
      Bottom_Lower = bottom_ci[1], Bottom_Upper = bottom_ci[2],
      Top_Lower = top_ci[1], Top_Upper = top_ci[2]
    )
  }

  # Calculate curve quality metrics
  calculate_curve_quality <- function(params, gof_results, plausibility_check = NULL, logIC50_ci = NULL) {
    tryCatch({
      span <- params[6]
      hill_slope <- params[4]
      max_slope <- -span * abs(hill_slope) * log(10) / 4

      quality_flags <- character()

      # Use same criteria as before for consistency
      if (abs(max_slope) < 5) quality_flags <- c(quality_flags, "Very shallow slope")
      else if (abs(max_slope) < 15) quality_flags <- c(quality_flags, "Shallow slope")
      if (abs(span) < 20) quality_flags <- c(quality_flags, "Small span")
      if (gof_results$R_squared < r_sqr_threshold) quality_flags <- c(quality_flags, "Low R2")

      # NOVA VERIFICACAO: Intervalo de confianca do LogIC50 muito amplo
      if (!is.null(logIC50_ci) && !any(is.na(logIC50_ci))) {
        ci_range <- abs(logIC50_ci[2] - logIC50_ci[1])
        if (ci_range > 0.666) {
          quality_flags <- c(quality_flags, "Wide logIC50 CI range")
        }
      }

      if (!is.null(plausibility_check) && plausibility_check$needs_correction) {
        quality_flags <- c(quality_flags, "Parameters corrected (biologically implausible)")
      }

      # Flag extreme Hill slopes, direction-aware:
      #   Inhibition curves: HillSlope is negative (e.g. -1 is standard)
      #     steep = HillSlope < -3  |  shallow = HillSlope > -0.3
      #   Activation curves: HillSlope is positive (e.g. +1 is standard)
      #     steep = HillSlope > 3   |  shallow = HillSlope < 0.3
      if (!is.na(hill_slope)) {
        if (hill_slope < 0) {
          # Inhibition curve
          if (hill_slope < -3)  quality_flags <- c(quality_flags, "Steep Hill slope")
          else if (hill_slope > -0.3) quality_flags <- c(quality_flags, "Shallow Hill slope")
        } else {
          # Activation curve
          if (hill_slope > 3)   quality_flags <- c(quality_flags, "Steep Hill slope")
          else if (hill_slope < 0.3)  quality_flags <- c(quality_flags, "Shallow Hill slope")
        }
      }

      list(
        quality = if (length(quality_flags) == 0) "Good curve" else paste(quality_flags, collapse = "; "),
        max_slope = max_slope
      )
    }, error = function(e) {
      list(quality = "Error in quality assessment", max_slope = NA)
    })
  }

  # Main analysis function for each compound pair
  analyze_single_pair <- function(pair_data, comp_name) {
    prepared <- prepare_and_validate_data(pair_data)
    if (!prepared$valid) {
      if (verbose) warning("Data validation failed for ", comp_name)
      return(create_empty_result(comp_name, "Validation failed"))
    }

    # Detect curve type once
    curve_type <- detect_curve_type(prepared$df_clean)

    fit <- try_robust_fit(prepared$df_clean, prepared$min_response,
                          prepared$max_response, prepared$approx_ic50, curve_type)

    if (is.null(fit)) {
      if (verbose) warning("Could not fit model for ", comp_name)
      result <- create_empty_result(comp_name, "Approximate (model failed)")
      result$parameters$Value <- c(prepared$min_response, prepared$max_response,
                                   prepared$approx_ic50, -1, 10^prepared$approx_ic50,
                                   prepared$max_response - prepared$min_response)
      result$curve_quality <- "Model failed - approximate parameters"
      result$success <- TRUE
      result$curve_type <- curve_type
      return(result)
    }

    params <- tryCatch(unname(stats::coef(fit)), error = function(e) NULL)
    if (is.null(params) || any(is.na(params))) {
      if (verbose) warning("Error extracting coefficients for ", comp_name)
      return(create_empty_result(comp_name, "Coefficient extraction failed"))
    }

    # Pass curve_type to avoid redundant detection
    order_correction <- correct_parameter_order(params, prepared$df_clean, curve_type)
    if (order_correction$was_corrected) {
      params <- order_correction$corrected_params
    }

    # Process results
    initial_params <- c(params, 10^params[3], params[2] - params[1])
    ci_results <- calculate_ci(fit)

    # If correct_parameter_order() flipped the HillSlope sign, the CI bounds
    # were computed on the original (wrong-sign) fit and must be negated.
    # The interval width is still valid; only the sign and order change.
    if (order_correction$was_corrected &&
        grepl("Hill Slope inconsistent", order_correction$correction_reason %||% "")) {
      hs_ci <- ci_results$HillSlope
      if (!any(is.na(hs_ci))) {
        ci_results$HillSlope <- sort(-hs_ci)  # negate and re-sort [lower, upper]
      }
    }

    gof_results <- calculate_goodness_of_fit(fit, prepared$df_clean)
    plausibility_check <- check_biological_plausibility(params, prepared$df_clean)

    # Apply corrections if needed
    if (plausibility_check$needs_correction) {
      if (verbose) {
        cat("Biological plausibility correction for ", comp_name, ":\n", sep = "")
        for (param in names(plausibility_check$correction_reasons)) {
          cat("  - ", param, ": ", plausibility_check$correction_reasons[[param]], "\n", sep = "")
        }
      }
      final_params <- plausibility_check$corrected_params
      ci_results <- recalculate_ci(ci_results, plausibility_check$corrections_applied)
    } else {
      final_params <- initial_params
    }

    curve_quality_info <- calculate_curve_quality(final_params, gof_results, plausibility_check, ci_results$LogIC50)

    list(
      parameters = data.frame(Parameter = PARAM_NAMES, Value = final_params, stringsAsFactors = FALSE),
      confidence_intervals = ci_results,
      goodness_of_fit = gof_results,
      curve_quality = curve_quality_info$quality,
      max_slope = curve_quality_info$max_slope,
      model = fit,
      success = TRUE,
      compound = comp_name,
      biological_plausibility_check = plausibility_check,
      parameter_order_correction = order_correction,
      data = prepared$df_clean,
      curve_type = curve_type
    )
  }

  # Normalize a data frame to 0-100% (same logic as fit_drc_3pl)
  normalize_dataframe <- function(df) {
    df %>% dplyr::mutate(dplyr::across(-1, ~ {
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

  # Data validation and normalization
  if (ncol(data) < 3) stop("Data must have at least 3 columns")

  original_data <- data

  if (verbose) message("Generating normalized data (0-100%)...")
  normalized_data <- normalize_dataframe(data)

  if (normalize) {
    if (verbose) message("Analysis will use: NORMALIZED data.")
    data <- normalized_data
  } else {
    if (verbose) message("Analysis will use: ORIGINAL data.")
  }

  # Main analysis loop
  # Group data columns by base compound name (strip trailing .2, .3, ... suffix).
  # This supports any number of replicates: 2 replicates (original behaviour),
  # 4 replicates (two merged plates), 6 replicates (three merged plates), etc.
  if (verbose) cat("Starting 4-parameter dose-response analysis...\n")

  c_names_data  <- colnames(data)[-1]               # data column names only
  base_names    <- sub("\\.\\d+$", "", c_names_data) # strip .2 / .3 / .4 ...
  compound_order <- unique(base_names)               # preserve first-appearance order
  n_compounds   <- length(compound_order)

  all_results <- lapply(seq_len(n_compounds), function(i) {
    if (verbose) cat(sprintf("\rProcessing compound %d/%d", i, n_compounds))

    base      <- compound_order[i]
    data_cols <- which(base_names == base) + 1L      # +1 for the conc column offset
    comp_name <- base

    analyze_single_pair(data[, c(1L, data_cols), drop = FALSE], comp_name)
  })
  if (verbose) cat("\n")

  # Create summary table
  summary_table <- do.call(rbind, lapply(all_results, function(result) {
    if (isTRUE(result$success)) {
      params <- unname(result$parameters$Value)
      gof <- result$goodness_of_fit
      ci <- result$confidence_intervals
      curve_type <- result$curve_type

      # Apply threshold for inhibition AND flat/unknown curves
      apply_threshold <- FALSE
      if (enforce_bottom_threshold && !is.na(params[1]) && params[1] >= bottom_threshold) {
        if (curve_type == "inhibition" || curve_type == "flat" || curve_type == "unknown") {
          apply_threshold <- TRUE
        }
      }

      data.frame(
        Compound = strsplit(result$compound, " \\| ")[[1]][1],
        Bottom = round(params[1], 3),
        Top = round(params[2], 3),
        LogIC50 = if (!apply_threshold) round(params[3], 3) else NA,
        HillSlope = if (!apply_threshold) round(params[4], 3) else NA,
        IC50 = if (!apply_threshold) format(params[5], scientific = TRUE) else NA,
        Bottom_Lower_95CI = if (!is.na(ci$Bottom_Lower)) round(ci$Bottom_Lower, 3) else NA_real_,
        Bottom_Upper_95CI = if (!is.na(ci$Bottom_Upper)) round(ci$Bottom_Upper, 3) else NA_real_,
        Top_Lower_95CI = if (!is.na(ci$Top_Lower)) round(ci$Top_Lower, 3) else NA_real_,
        Top_Upper_95CI = if (!is.na(ci$Top_Upper)) round(ci$Top_Upper, 3) else NA_real_,
        LogIC50_Lower_95CI = if (!apply_threshold && !is.na(ci$LogIC50[1])) round(ci$LogIC50[1], 3) else NA,
        LogIC50_Upper_95CI = if (!apply_threshold && !is.na(ci$LogIC50[2])) round(ci$LogIC50[2], 3) else NA,
        HillSlope_Lower_95CI = if (!apply_threshold && !is.na(ci$HillSlope[1])) round(ci$HillSlope[1], 3) else NA,
        HillSlope_Upper_95CI = if (!apply_threshold && !is.na(ci$HillSlope[2])) round(ci$HillSlope[2], 3) else NA,
        IC50_Lower_95CI = if (!apply_threshold && !is.na(ci$IC50[1])) format(ci$IC50[1], scientific = TRUE) else NA,
        IC50_Upper_95CI = if (!apply_threshold && !is.na(ci$IC50[2])) format(ci$IC50[2], scientific = TRUE) else NA,
        Span = round(params[6], 3),
        R_squared = round(gof$R_squared, 3),
        Syx = round(gof$Syx, 3),
        Sum_of_Squares = round(gof$Sum_of_Squares, 3),
        Degrees_of_Freedom = gof$Degrees_of_Freedom,
        Max_Slope = round(result$max_slope %||% NA, 3),
        Curve_Quality = result$curve_quality %||% "Not assessed",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Compound = strsplit(result$compound, " \\| ")[[1]][1],
        Bottom = NA, Top = NA, LogIC50 = NA, HillSlope = NA, IC50 = NA,
        Bottom_Lower_95CI = NA, Bottom_Upper_95CI = NA,
        Top_Lower_95CI = NA, Top_Upper_95CI = NA,
        LogIC50_Lower_95CI = NA, LogIC50_Upper_95CI = NA,
        HillSlope_Lower_95CI = NA, HillSlope_Upper_95CI = NA,
        IC50_Lower_95CI = NA, IC50_Upper_95CI = NA, Span = NA,
        R_squared = NA, Syx = NA, Sum_of_Squares = NA, Degrees_of_Freedom = NA,
        Max_Slope = NA, Curve_Quality = "Fit failed",
        stringsAsFactors = FALSE
      )
    }
  }))

  # Create final_summary_table (transposed version)
  if (nrow(summary_table) > 0) {
    compound_names <- summary_table$Compound
    transposed_data <- as.data.frame(t(summary_table[, -1]))
    colnames(transposed_data) <- compound_names
    final_summary_table <- transposed_data

  } else {
    final_summary_table <- data.frame()
  }

  # Identify compounds affected by threshold
  threshold_affected <- character()
  if (enforce_bottom_threshold) {
    for (result in all_results) {
      if (isTRUE(result$success) &&
          !is.na(result$parameters$Value[1]) &&
          result$parameters$Value[1] >= bottom_threshold) {
        curve_type <- result$curve_type
        if (curve_type == "inhibition" || curve_type == "flat" || curve_type == "unknown") {
          comp_name <- strsplit(result$compound, " \\| ")[[1]][1]
          threshold_affected <- c(threshold_affected, comp_name)
        }
      }
    }
  }

  # Count order corrections
  order_corrections <- sum(sapply(all_results, function(x) {
    if (!is.null(x$parameter_order_correction)) {
      x$parameter_order_correction$was_corrected
    } else {
      FALSE
    }
  }))

  # Print summary statistics
  if (verbose) {
    successful <- sum(!is.na(summary_table$IC50))
    total <- nrow(summary_table)
    success_rate <- round(successful / total * 100, 1)
    threshold_count <- length(threshold_affected)

    cat("\n", strrep("=", 50), "\n", sep = "")
    cat("4-PARAMETER DOSE-RESPONSE ANALYSIS COMPLETED SUCCESSFULLY!\n")
    cat(strrep("=", 50), "\n")
    cat("SUMMARY STATISTICS:\n")
    cat("  . Compounds analyzed: ", total, "\n")
    cat("  . Successful fits: ", successful, " (", success_rate, "%)\n", sep = "")
    cat("  . Failed fits: ", total - successful, "\n")

    if (order_corrections > 0) {
      cat("  . Parameter order corrections: ", order_corrections, "\n")
    }

    if (enforce_bottom_threshold && threshold_count > 0) {
      cat("  . IC50 values excluded (Bottom <=", bottom_threshold, "): ", threshold_count, "\n", sep = "")

      cat("\nCOMPOUNDS WITH EXCLUDED IC50 VALUES:\n")
      for (i in seq_along(threshold_affected)) {
        bottom_val <- summary_table$Bottom[summary_table$Compound == threshold_affected[i]]
        cat("  . ", threshold_affected[i], " (Bottom = ",
            round(bottom_val, 1), ")\n", sep = "")
      }
    }

    cat("  . Problematic curves: ", sum(grepl("shallow|Low R2|Small span|Wide CI range", summary_table$Curve_Quality)), "\n")

    if ("Curve_Quality" %in% names(summary_table)) {
      cat("\nCURVE QUALITY DISTRIBUTION:\n")
      quality_counts <- table(summary_table$Curve_Quality)
      for (quality in names(sort(quality_counts, decreasing = TRUE))) {
        count <- quality_counts[quality]
        cat("  . ", sprintf("%-35s: %d (%.1f%%)", quality, count, count/total*100), "\n")
      }
    }
    cat(strrep("=", 50), "\n\n")
  }

  # Save results to file
  if (!is.null(output_file)) {
    file_ext <- tolower(tools::file_ext(output_file))

    if (file_ext == "xlsx" && requireNamespace("openxlsx", quietly = TRUE)) {
      wb <- openxlsx::createWorkbook()

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

      openxlsx::addWorksheet(wb, "Summary")
      openxlsx::writeData(wb, "Summary", summary_table)

      openxlsx::addWorksheet(wb, "Normalized_Data")
      openxlsx::writeData(wb, "Normalized_Data", normalized_data)

      openxlsx::addWorksheet(wb, "Original_Data")
      openxlsx::writeData(wb, "Original_Data", original_data)

      openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

      if (verbose) {
        message("Results saved to Excel: ", output_file)
        message("  1. Final_Summary (Transposed)")
        message("  2. Summary (Detailed)")
        message("  3. Normalized_Data (0-100%)")
        message("  4. Original_Data (Raw)")
      }
    } else {
      if (file_ext == "xlsx") {
        output_file <- sub("\\.xlsx$", ".csv", output_file, ignore.case = TRUE)
        if (verbose) warning("Falling back to CSV format...")
      }
      utils::write.csv(summary_table, output_file, row.names = FALSE)
      if (verbose) cat("Results saved to:", output_file, "\n")
    }
  }

  # Return results object - MODIFICADA
  list(
    summary_table = summary_table,
    final_summary_table = final_summary_table,  # NOVA TABELA
    detailed_results = all_results,
    n_compounds = n_compounds,
    successful_fits = sum(!is.na(summary_table$IC50)),
    normalized = normalize,
    original_data = original_data,
    normalized_data = normalized_data,
    used_normalized_data = normalize,
    processed_data = data,
    threshold_settings = if (enforce_bottom_threshold) {
      list(
        enforce_bottom_threshold = enforce_bottom_threshold,
        bottom_threshold = bottom_threshold,
        affected_compounds = threshold_affected
      )
    } else {
      NULL
    },
    parameter_order_corrections = order_corrections
  )
}

