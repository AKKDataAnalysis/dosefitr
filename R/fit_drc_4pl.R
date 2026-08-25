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
#' @param normalize Logical indicating whether to normalize response data to 0-100\%
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
#' @param bottom_limits_nanobret_inhibition Numeric length-2 vector
#'   \code{c(lower, upper)} with the plausible Bottom range for
#'   \code{assay_type = "nanobret"} inhibition/flat/unknown curves.
#'   Default \code{c(-100, 600)}.
#' @param bottom_limits_nanobret_activation Numeric length-2 vector with the
#'   plausible Bottom range for \code{assay_type = "nanobret"} activation
#'   curves. Default \code{c(-100, Inf)}.
#' @param top_limits_nanobret_inhibition Numeric length-2 vector with the
#'   plausible Top range for \code{assay_type = "nanobret"}
#'   inhibition/flat/unknown curves. Default \code{c(0, 700)}.
#' @param top_limits_nanobret_activation Numeric length-2 vector with the
#'   plausible Top range for \code{assay_type = "nanobret"} activation curves.
#'   Default \code{c(0, 700)}.
#' @param bottom_limits_viability_inhibition Numeric length-2 vector with the
#'   plausible Bottom range for \code{assay_type = "viability"} (only applied
#'   when \code{normalize = TRUE}) inhibition/flat/unknown curves.
#'   Default \code{c(-20, 60)}.
#' @param bottom_limits_viability_activation Numeric length-2 vector with the
#'   plausible Bottom range for \code{assay_type = "viability"} (only applied
#'   when \code{normalize = TRUE}) activation curves. Default \code{c(-20, 60)}.
#' @param top_limits_viability_inhibition Numeric length-2 vector with the
#'   plausible Top range for \code{assay_type = "viability"} (only applied
#'   when \code{normalize = TRUE}) inhibition/flat/unknown curves.
#'   Default \code{c(50, 130)}.
#' @param top_limits_viability_activation Numeric length-2 vector with the
#'   plausible Top range for \code{assay_type = "viability"} (only applied
#'   when \code{normalize = TRUE}) activation curves. Default \code{c(50, 130)}.
#' @param logIC50_limits Numeric length-2 vector with the plausible LogIC50
#'   range, shared by all assay types. Fits outside this range have LogIC50
#'   set to \code{NA}. Default \code{c(-20, 5)}.
#' @param hill_slope_limits Numeric length-2 vector with the plausible
#'   HillSlope range, shared by all assay types. Fits outside this range are
#'   reset to the default sign convention (\eqn{\pm 1}) clamped to these
#'   limits. Default \code{c(-5, 5)}.
#'
#'   All limit arguments must be numeric vectors of length 2 with
#'   \code{lower <= upper}; \code{-Inf}/\code{Inf} are allowed as open bounds.
#'   Invalid values cause an immediate error before any fitting.
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
#' \deqn{Response = Bottom + \frac{Top - Bottom}{1 + 10^{(LogIC50 - log_{10}(inhibitor)) \times HillSlope}}}
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
#' \item Parameter corrections: the quality string states the reason directly:
#'   "Implausible fit" (converged but out-of-limits), "No dose response
#'   (inactive compound)" (diverged/unidentifiable fit on flat data), or
#'   "Partial response at highest doses only" (diverged/unidentifiable fit on
#'   a response that only begins at the highest doses).
#' \item Constrained refit: for the "partial response" case the fit is
#'   retried with constraints (plateau bounds, then Bottom fixed at 0 for
#'   inhibition curves).  A successful refit is flagged "IC50 from
#'   constrained fit (plateaus bounded)" or "IC50 from constrained fit
#'   (Bottom fixed at 0)" — the IC50 is real but rests on the constraint
#'   assumption.
#' \item logIC50 range (>0.666 flagged)
#' }
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_4pl")
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
#' fit4 <- fit_drc_4pl(data = mrt, normalize = FALSE, verbose = FALSE)
#' head(fit4$summary_table[, c("Compound", "LogIC50", "R_squared")])
#' }
#' @export


fit_drc_4pl <- function(data, output_file = NULL, normalize = FALSE, verbose = TRUE,
                        enforce_bottom_threshold = FALSE, bottom_threshold = 60,
                        r_sqr_threshold = 0.8,
                        assay_type = "nanobret",
                        bottom_limits_nanobret_inhibition = c(-100, 600),
                        bottom_limits_nanobret_activation = c(-100, Inf),
                        top_limits_nanobret_inhibition = c(0, 700),
                        top_limits_nanobret_activation = c(0, 700),
                        bottom_limits_viability_inhibition = c(-20, 60),
                        bottom_limits_viability_activation = c(-20, 60),
                        top_limits_viability_inhibition = c(50, 130),
                        top_limits_viability_activation = c(50, 130),
                        logIC50_limits = c(-20, 5),
                        hill_slope_limits = c(-5, 5)) {

  # Fail fast on malformed plausibility-limit arguments (before any fitting).
  .validate_all_limit_args(environment())

  # Check and load required packages
  required_packages <- c("dplyr", "stats", "grDevices")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages are not installed: ", paste(missing_packages, collapse = ", "),
         "\nPlease install using: install.packages(c(", 
         paste0("\"", missing_packages, "\"", collapse = ", "), "))")
  }
  
  # Constants for 4-parameter model
  PARAM_NAMES <- c("Bottom", "Top", "LogIC50", "HillSlope", "IC50", "Span")
  
  # NULL coalescing operator
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b
  
  # 4-parameter logistic model
  four_param_model <- function(log_inhibitor, Bottom, Top, LogIC50, HillSlope) {
    Bottom + (Top - Bottom) / (1 + 10^((LogIC50 - log_inhibitor) * HillSlope))
  }
  
  # Detect curve type based on data pattern.
  # Uses a relative threshold (max(15, range * 0.15)) on the head/tail
  # mean difference -- identical to fit_drc_3pl.
  # Direction classification is delegated to the SHARED noise-aware,
  # scale-robust classifier (single source of truth in shared_fit_strategy.R).
  # This replaces the former hardcoded absolute floor of 15, which mislabeled
  # real raw-scale curves as flat, and guarantees 3PL/4PL agree on direction.
  detect_curve_type <- function(data) detect_curve_type_shared(data)
  
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
    
    c(corrected_params[1:4], if (!is.na(corrected_params[3])) 10^corrected_params[3] else NA_real_, corrected_params[2] - corrected_params[1])
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
  # Limits depend on assay_type and normalize, matching fit_drc_3pl behaviour.
  #   anything else               : check disabled.
  # Select the Bottom/Top plausibility limits for this assay type and curve
  # direction.  Single source of truth: used by check_biological_plausibility()
  # for post-hoc correction AND by try_constrained_fit() as port bounds.
  select_plausibility_limits <- function(curve_type) {
    if (assay_type == "nanobret") {
      bottom_limits <- if (curve_type == "activation") {
        bottom_limits_nanobret_activation
      } else {
        bottom_limits_nanobret_inhibition
      }
      top_limits <- if (curve_type == "activation") {
        top_limits_nanobret_activation
      } else {
        top_limits_nanobret_inhibition
      }
    } else {
      # viability + normalize=TRUE
      bottom_limits <- if (curve_type == "activation") {
        bottom_limits_viability_activation
      } else {
        bottom_limits_viability_inhibition
      }
      top_limits <- if (curve_type == "activation") {
        top_limits_viability_activation
      } else {
        top_limits_viability_inhibition
      }
    }
    list(bottom_limits = bottom_limits, top_limits = top_limits)
  }

  check_biological_plausibility <- function(params, data) {
    if (assay_type == "viability" && !normalize)
      return(list(needs_correction = FALSE))

    if (!assay_type %in% c("nanobret", "viability"))
      return(list(needs_correction = FALSE))

    responses <- data$response[!is.na(data$response)]
    exp_min <- min(responses, na.rm = TRUE)
    exp_max <- max(responses, na.rm = TRUE)

    curve_type <- detect_curve_type(data)

    lims <- select_plausibility_limits(curve_type)
    bottom_limits <- lims$bottom_limits
    top_limits    <- lims$top_limits
    # logIC50_limits and hill_slope_limits are used directly (function args).
    
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
      corrections$LogIC50 <- NA_real_
      reasons$LogIC50 <- sprintf("Biologically implausible (%.2f). Setting to NA.",
                                 params[3])
    }
    
    # Check HillSlope
    if (params[4] < hill_slope_limits[1] || params[4] > hill_slope_limits[2] || !is.finite(params[4])) {
      default_hill <- if (curve_type == "activation") 1 else -1
      corrections$HillSlope <- max(hill_slope_limits[1], min(hill_slope_limits[2], default_hill))
      reasons$HillSlope <- sprintf("Biologically implausible (%.2f). Using default: %.2f",
                                   params[4], default_hill)
    }
    
    if (length(corrections) > 0) {
      # Classify WHY the correction fired, for a more informative quality
      # label.  Two different situations both end up here:
      #   (a) the nls fit converged but a parameter sits just outside the
      #       plausibility limits -> genuine "biologically implausible";
      #   (b) the fit DIVERGED or is unidentifiable (asymptotes exploded far
      #       beyond the data, or LogIC50 landed outside the tested doses),
      #       which typically means the compound is simply inactive, or the
      #       response only begins at the edge of the dose range.
      doses <- data$log_inhibitor[!is.na(data$log_inhibitor)]
      resp_range <- exp_max - exp_min
      asymp_lo <- exp_min - 3 * resp_range
      asymp_hi <- exp_max + 3 * resp_range
      fit_diverged <- !is.finite(params[1]) || !is.finite(params[2]) ||
                      params[1] < asymp_lo || params[1] > asymp_hi ||
                      params[2] < asymp_lo || params[2] > asymp_hi
      ic50_outside_range <- !is.finite(params[3]) ||
                            params[3] < min(doses) - 1 || params[3] > max(doses) + 1

      # Does the data itself show a response at the two highest doses?
      # Threshold combines the observed data range with a minimum 20% change
      # relative to the plateau level, so noisy flat data (small range) does
      # not register as a response.
      hi_doses  <- utils::tail(sort(unique(doses)), 2)
      edge_mean <- mean(data$response[data$log_inhibitor %in% hi_doses], na.rm = TRUE)
      base_mean <- mean(data$response[!data$log_inhibitor %in% hi_doses], na.rm = TRUE)
      min_effect <- max(0.25 * resp_range, 0.20 * abs(base_mean))
      edge_response <- if (curve_type == "activation") {
        is.finite(edge_mean) && is.finite(base_mean) &&
          (edge_mean - base_mean) > min_effect
      } else {
        is.finite(edge_mean) && is.finite(base_mean) &&
          (base_mean - edge_mean) > min_effect
      }

      # Whole-window trend test: does the observed data show a real
      # direction-consistent response between the two lowest and the two
      # highest doses?  Used to distinguish a genuinely flat compound
      # ("No dose response") from one whose response happens OUTSIDE the
      # tested window (e.g. an inhibition that is already complete at the
      # lowest dose -> "IC50 below tested range"), which a diverged or
      # out-of-range fit otherwise mislabels as inactivity.
      lo_doses <- utils::head(sort(unique(doses)), 2)
      lo_mean  <- mean(data$response[data$log_inhibitor %in% lo_doses], na.rm = TRUE)
      trend_response <- if (curve_type == "activation") {
        is.finite(lo_mean) && is.finite(edge_mean) &&
          (edge_mean - lo_mean) > min_effect
      } else {
        is.finite(lo_mean) && is.finite(edge_mean) &&
          (lo_mean - edge_mean) > min_effect
      }

      # Near-flat converged fit guard.  When the fit converged, LogIC50 sits
      # inside the tested range, and the data show no edge response, the
      # firing correction is an asymptote of an essentially flat fit resting
      # just outside the plausibility limits (classic case: a vehicle
      # control fitted with Bottom ~ 100 on a normalized viability scale,
      # whose upper limit is 60).  Clamping that asymptote would MANUFACTURE
      # a sigmoid the data do not support -- the plotted curve would dive
      # away from points it actually fits.  Keep the converged parameters
      # and ask the caller to reclassify the curve as flat (no dose
      # response, so the batch layer reports IC50 as N/D).  Genuine
      # partial/shallow curves are protected twice: by the edge_response
      # test above and by the fitted-change test here (a real 100 -> 70
      # curve moves ~30 points across the tested doses, far above the
      # threshold).
      fitted_change <- abs(
        four_param_model(max(doses), params[1], params[2], params[3], params[4]) -
        four_param_model(min(doses), params[1], params[2], params[3], params[4])
      )
      near_flat_fit <- is.finite(fitted_change) &&
                       fitted_change < max(0.25 * resp_range, 0.10 * abs(base_mean))

      # Tight-fit escape: a converged curve that tracks the data closely
      # (fitted change across the window >= 10x the residual noise AND
      # R2 >= 0.9) is a real, well-resolved response -- merely shallow or
      # sitting on a high baseline (e.g. NanoBRET ratios ~650, where
      # 0.10 x baseline = 67 units dwarfs a genuine 50-unit sigmoid with
      # ~1-unit noise).  Such fits deserve the parameter correction, not a
      # flat reclassification.
      ok_row   <- !is.na(data$log_inhibitor) & !is.na(data$response)
      fv       <- four_param_model(data$log_inhibitor[ok_row],
                                   params[1], params[2], params[3], params[4])
      resid    <- data$response[ok_row] - fv
      resid_sd <- stats::sd(resid)
      ss_res   <- sum(resid^2)
      ss_tot   <- sum((data$response[ok_row] - mean(data$response[ok_row]))^2)
      r2_fit   <- if (ss_tot > 0) 1 - ss_res / ss_tot else 0
      tight_fit <- is.finite(resid_sd) && is.finite(fitted_change) &&
                   fitted_change >= 10 * resid_sd && r2_fit >= 0.9

      if (!fit_diverged && !ic50_outside_range && !edge_response && near_flat_fit && !tight_fit) {
        return(list(
          needs_correction       = FALSE,
          reclassify_flat        = TRUE,
          fit_diverged           = fit_diverged,
          ic50_outside_range     = ic50_outside_range,
          edge_response_detected = edge_response,
          fit_label              = "No dose response (inactive compound)"
        ))
      }

      fit_label <- if ((fit_diverged || ic50_outside_range) && !edge_response) {
        if (ic50_outside_range && trend_response) {
          # Real response visible in the data, but its midpoint lies outside
          # the tested window -- say so instead of claiming inactivity.
          if (is.finite(params[3]) && params[3] < min(doses)) {
            "IC50 below tested range (N/D)"
          } else {
            "IC50 above tested range (N/D)"
          }
        } else {
          "No dose response (inactive compound)"
        }
      } else if ((fit_diverged || ic50_outside_range) && edge_response) {
        "Partial response at highest doses only"
      } else {
        "Implausible fit"
      }

      return(list(
        corrected_params = recalculate_dependent_params(params, corrections),
        corrections_applied = corrections,
        correction_reasons = reasons,
        needs_correction = TRUE,
        fit_diverged = fit_diverged,
        ic50_outside_range = ic50_outside_range,
        edge_response_detected = edge_response,
        fit_label = fit_label,
        # Bug A: a diverged / out-of-range fit with no edge response is labelled
        # "No dose response" but kept its activation/inhibition curve_type, so it
        # escaped the batch N/D marking (which keys on curve_type == "flat").
        # The same applies to the out-of-window N/D labels.
        reclassify_flat = fit_label %in% c("No dose response (inactive compound)",
                                           "IC50 below tested range (N/D)",
                                           "IC50 above tested range (N/D)")
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
  
  # Robust nonlinear fitting with multiple strategies.
  #
  # Starting values come from the SHARED strategy module (make_start_grid),
  # so this OLS/nls engine and the ROUT outlier engine seed their fits from
  # one source of truth. Unlike the previous "first strategy that converges
  # wins" behaviour, this tries ALL starts and keeps the BEST-fitting one
  # (lowest residual sum of squares) -- more robust for steep curves, which
  # matters because these fits sit underneath outlier detection.
  try_robust_fit <- function(df_clean, min_resp, max_resp, approx_ic50, curve_type) {
    if (nrow(df_clean) < 5) return(NULL)

    control_configs <- list(
      default = stats::nls.control(maxiter = 500, tol = 1e-04, minFactor = 1/4096, warnOnly = TRUE),
      relaxed = stats::nls.control(maxiter = 1000, tol = 1e-03, minFactor = 1/1024, warnOnly = TRUE)
    )

    # Build the shared starting-value grid on the drc parameterisation
    # (Bottom, Top, LogIC50, HillSlope). x here is log10 concentration.
    build_grid <- function(direction) {
      make_start_grid(df_clean$log_inhibitor, df_clean$response,
                      direction = direction, param = "drc")
    }
    if (curve_type == "inhibition") {
      start_strategies <- build_grid("inhibition")
    } else if (curve_type == "activation") {
      start_strategies <- build_grid("activation")
    } else {
      # For flat or unknown curves: try BOTH directions.
      start_strategies <- c(build_grid("inhibition"), build_grid("activation"))
    }

    # Fit from every start; keep the converged fit with the lowest RSS.
    fit_one <- function(start_vals, control) {
      tryCatch(
        stats::nls(
          response ~ four_param_model(log_inhibitor, Bottom, Top, LogIC50, HillSlope),
          data = df_clean,
          start = list(Bottom = start_vals[["Bottom"]], Top = start_vals[["Top"]],
                       LogIC50 = start_vals[["LogIC50"]], HillSlope = start_vals[["HillSlope"]]),
          control = control, algorithm = "port"
        ),
        error = function(e) NULL
      )
    }

    best <- NULL; best_rss <- Inf
    for (start_vals in start_strategies) {
      fit <- fit_one(start_vals, control_configs$default)
      if (is.null(fit)) next
      rss <- tryCatch(sum(stats::resid(fit)^2), error = function(e) Inf)
      if (is.finite(rss) && rss < best_rss) { best <- fit; best_rss <- rss }
    }
    if (!is.null(best)) return(best)

    # Final attempt with relaxed control from the first start (fallback).
    fit_one(start_strategies[[1]], control_configs$relaxed)
  }

  # Constrained refit for diverged "edge response" fits.
  #
  # When the unconstrained fit diverges but the data DO respond at the
  # highest doses, the plateau past the edge is unidentifiable and nls
  # wanders to absurd values.  This automates GraphPad's recommended remedy
  # (constrain a parameter), weakest assumption first:
  #   1. "bounds": refit with port box constraints -- Bottom/Top bounded to
  #      the plausibility limits, LogIC50 to the tested doses +/- 1 log,
  #      HillSlope to hill_slope_limits;
  #   2. "bottom_fixed": if the bounded fit still fails the plausibility
  #      check, fix Bottom = 0 (complete killing; inhibition curves only)
  #      and refit the remaining three parameters.
  # Returns list(fit = <nls>, method = "bounds"|"bottom_fixed") when the
  # constrained fit passes check_biological_plausibility(), else NULL.
  try_constrained_fit <- function(df_clean, curve_type) {
    if (nrow(df_clean) < 5) return(NULL)

    lims  <- select_plausibility_limits(curve_type)
    doses <- df_clean$log_inhibitor[!is.na(df_clean$log_inhibitor)]
    resp  <- df_clean$response[!is.na(df_clean$response)]
    resp_range <- max(resp) - min(resp)

    # port bounds must be finite; replace any infinite plausibility limit
    # with a generous data-driven sentinel.
    lo_b <- lims$bottom_limits[1]; hi_b <- lims$bottom_limits[2]
    lo_t <- lims$top_limits[1];    hi_t <- lims$top_limits[2]
    if (!is.finite(lo_b)) lo_b <- min(resp) - 10 * resp_range
    if (!is.finite(hi_b)) hi_b <- max(resp) + 10 * resp_range
    if (!is.finite(lo_t)) lo_t <- min(resp) - 10 * resp_range
    if (!is.finite(hi_t)) hi_t <- max(resp) + 10 * resp_range

    lower4 <- c(Bottom = lo_b, Top = lo_t,
                LogIC50 = min(doses) - 1, HillSlope = hill_slope_limits[1])
    upper4 <- c(Bottom = hi_b, Top = hi_t,
                LogIC50 = max(doses) + 1, HillSlope = hill_slope_limits[2])

    # Same multi-start grid as the unconstrained fit, clamped into bounds
    # (nudged slightly off the exact boundary so port starts feasible).
    clamp_start <- function(s) {
      v <- pmin(pmax(s[names(lower4)], lower4), upper4)
      eps <- 1e-6 * pmax(1, abs(upper4 - lower4))
      pmin(pmax(v, lower4 + eps), upper4 - eps)
    }
    if (curve_type == "inhibition") {
      grids <- make_start_grid(doses, df_clean$response, direction = "inhibition", param = "drc")
    } else if (curve_type == "activation") {
      grids <- make_start_grid(doses, df_clean$response, direction = "activation", param = "drc")
    } else {
      grids <- c(make_start_grid(doses, df_clean$response, direction = "inhibition", param = "drc"),
                 make_start_grid(doses, df_clean$response, direction = "activation", param = "drc"))
    }

    control <- stats::nls.control(maxiter = 500, warnOnly = TRUE)

    fit_bounded <- function(s) {
      tryCatch(
        stats::nls(
          response ~ four_param_model(log_inhibitor, Bottom, Top, LogIC50, HillSlope),
          data = df_clean,
          start = list(Bottom = s[["Bottom"]], Top = s[["Top"]],
                       LogIC50 = s[["LogIC50"]], HillSlope = s[["HillSlope"]]),
          control = control, algorithm = "port", lower = lower4, upper = upper4
        ),
        error = function(e) NULL
      )
    }

    best <- NULL; best_rss <- Inf
    for (s0 in grids) {
      f <- fit_bounded(clamp_start(s0))
      if (is.null(f)) next
      rss <- tryCatch(sum(stats::resid(f)^2), error = function(e) Inf)
      if (is.finite(rss) && rss < best_rss) { best <- f; best_rss <- rss }
    }

    # Accept the bounded fit only if its parameters need no correction.
    if (!is.null(best)) {
      p <- tryCatch(unname(stats::coef(best)), error = function(e) NULL)
      if (!is.null(p) && length(p) == 4 && !any(is.na(p)) &&
          !isTRUE(check_biological_plausibility(p, df_clean)$needs_correction)) {
        return(list(fit = best, method = "bounds"))
      }
    }

    # Step 2: fix Bottom = 0 (inhibition only -- no sensible constant for an
    # unreached upper plateau on activation curves).
    if (curve_type != "inhibition") return(NULL)

    fit_fixed <- function(s) {
      tryCatch(
        stats::nls(
          response ~ four_param_model(log_inhibitor, 0, Top, LogIC50, HillSlope),
          data = df_clean,
          start = list(Top = s[["Top"]], LogIC50 = s[["LogIC50"]],
                       HillSlope = s[["HillSlope"]]),
          control = control, algorithm = "port",
          lower = lower4[c("Top", "LogIC50", "HillSlope")],
          upper = upper4[c("Top", "LogIC50", "HillSlope")]
        ),
        error = function(e) NULL
      )
    }

    best2 <- NULL; best_rss2 <- Inf
    for (s0 in grids) {
      f <- fit_fixed(clamp_start(s0))
      if (is.null(f)) next
      rss <- tryCatch(sum(stats::resid(f)^2), error = function(e) Inf)
      if (is.finite(rss) && rss < best_rss2) { best2 <- f; best_rss2 <- rss }
    }
    if (is.null(best2)) return(NULL)

    p2 <- tryCatch(unname(stats::coef(best2)), error = function(e) NULL)
    if (is.null(p2) || length(p2) != 3 || any(is.na(p2))) return(NULL)
    if (!isTRUE(check_biological_plausibility(c(0, p2), df_clean)$needs_correction)) {
      return(list(fit = best2, method = "bottom_fixed"))
    }
    NULL
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
  # Manual profile-likelihood extension for bounds profile() cannot reach.
  # profile.nls() often stalls on near-singular gradients (its nested refits
  # die before tau crosses the 95% threshold) and returns NA for a bound that
  # actually exists.  When a side comes back NA, redo the profiling manually
  # on that side: fix the parameter on an outward grid, refit the remaining
  # parameters, and interpolate where |tau| = sqrt((RSS - RSS0) / s2) crosses
  # qnorm(0.975).  A side that never crosses within generous limits stays
  # NA -- an honest "unbounded", not a fabricated number.
  extend_profile_side <- function(fit, df_clean, param_name, side) {
    est_all <- stats::coef(fit)
    if (!param_name %in% names(est_all)) return(NA_real_)
    est  <- unname(est_all[param_name])
    rss0 <- sum(stats::resid(fit)^2)
    s2   <- rss0 / stats::df.residual(fit)
    thr  <- stats::qnorm(0.975)

    resp  <- df_clean$response
    doses <- df_clean$log_inhibitor
    rng   <- diff(range(resp, na.rm = TRUE))
    if (!is.finite(rng) || rng <= 0) rng <- max(abs(resp), na.rm = TRUE)

    limit <- switch(param_name,
      Bottom    = ,
      Top       = if (side == "lower") min(resp, na.rm = TRUE) - rng
                  else max(resp, na.rm = TRUE) + rng,
      LogIC50   = if (side == "lower") min(doses, na.rm = TRUE) - 3
                  else max(doses, na.rm = TRUE) + 3,
      HillSlope = if (side == "lower") -max(abs(hill_slope_limits))
                  else max(abs(hill_slope_limits)),
      NA_real_)
    if (is.na(limit) || (side == "upper" && limit <= est) ||
        (side == "lower" && limit >= est)) return(NA_real_)

    all_params  <- c("Bottom", "Top", "LogIC50", "HillSlope")
    free_params <- setdiff(all_params, param_name)

    # tau of the fit with param_name fixed at v (constant baked into formula)
    tau_at <- function(v) {
      rhs <- all_params
      rhs[rhs == param_name] <- format(v, digits = 12)
      fml <- stats::as.formula(
        sprintf("response ~ four_param_model(log_inhibitor, %s)",
                paste(rhs, collapse = ", ")))
      # suppressWarnings: nested-fit convergence notes are diagnostic noise
      # here -- a stalled nested fit just yields a conservative tau via RSS.
      nested <- tryCatch(
        suppressWarnings(
          stats::nls(fml, data = df_clean, start = as.list(est_all[free_params]),
                     algorithm = "port",
                     control = stats::nls.control(maxiter = 500, warnOnly = TRUE))),
        error = function(e) NULL)
      if (is.null(nested)) return(NA_real_)
      sign(v - est) * sqrt(max(sum(stats::resid(nested)^2) - rss0, 0) / s2)
    }

    # Outward scan for the first crossing of +/-thr
    fracs <- c(0.03, 0.06, 0.10, 0.16, 0.24, 0.34, 0.46, 0.60, 0.76, 1.0)
    grid  <- est + (limit - est) * fracs
    taus  <- vapply(grid, tau_at, numeric(1))
    crossed <- which(if (side == "upper") taus >= thr else taus <= -thr)
    if (length(crossed) == 0) return(NA_real_)
    i <- crossed[1]
    v_lo <- if (i == 1) est else grid[i - 1]  # bracket: no crossing .. crossing
    v_hi <- grid[i]

    # Bisection refinement
    for (k in seq_len(6)) {
      v_mid <- (v_lo + v_hi) / 2
      t_mid <- tau_at(v_mid)
      if (is.na(t_mid)) break
      if (if (side == "upper") t_mid >= thr else t_mid <= -thr) v_hi <- v_mid
      else v_lo <- v_mid
    }
    (v_lo + v_hi) / 2
  }

  calculate_ci <- function(fit, df_clean = NULL) {
    safe_ci <- function(param_name) {
      ci <- tryCatch({
        prof <- stats::profile(fit, which = param_name, alpha = 0.05)
        unname(stats::confint(prof, level = 0.95))
      }, error = function(e) {
        tryCatch({
          se  <- summary(fit)$coefficients[param_name, "Std. Error"]
          est <- stats::coef(fit)[param_name]
          z   <- stats::qnorm(0.975)
          c(est - z * se, est + z * se)
        }, error = function(e2) c(NA, NA))
      })
      # Recover bounds profile() dropped (stalled nested refits), when they exist
      if (!is.null(df_clean) && any(is.na(ci))) {
        if (is.na(ci[1])) {
          lo <- extend_profile_side(fit, df_clean, param_name, "lower")
          if (!is.na(lo)) ci[1] <- lo
        }
        if (is.na(ci[2])) {
          hi <- extend_profile_side(fit, df_clean, param_name, "upper")
          if (!is.na(hi)) ci[2] <- hi
        }
      }
      ci
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
  calculate_curve_quality <- function(params, gof_results, plausibility_check = NULL, logIC50_ci = NULL,
                                      curve_type = "unknown") {
    tryCatch({
      span <- params[6]
      hill_slope <- params[4]
      max_slope <- -span * abs(hill_slope) * log(10) / 4

      quality_flags <- character()

      # A curve detected as flat is never a "Good curve": lead with the
      # no-response label regardless of how clean the fit looks, and drop the
      # generic "Implausible fit" (the flat classification is more specific).
      # When the plausibility check supplied a more specific N/D reason
      # (out-of-window response), lead with that instead.
      nd_labels <- c("No dose response (inactive compound)",
                     "IC50 below tested range (N/D)",
                     "IC50 above tested range (N/D)")
      if (identical(curve_type, "flat")) {
        lead_label <- if (!is.null(plausibility_check$fit_label) &&
                          plausibility_check$fit_label %in% nd_labels) {
          plausibility_check$fit_label
        } else {
          "No dose response (inactive compound)"
        }
        quality_flags <- c(lead_label, quality_flags)
      }

      # Use same criteria as before for consistency
      if (abs(max_slope) < 5) quality_flags <- c(quality_flags, "Very shallow slope")
      else if (abs(max_slope) < 15) quality_flags <- c(quality_flags, "Shallow slope")
      if (abs(span) < 20) quality_flags <- c(quality_flags, "Small span")
      if (gof_results$R_squared < r_sqr_threshold) quality_flags <- c(quality_flags, "Low R2")
      
      if (!is.null(logIC50_ci) && !any(is.na(logIC50_ci))) {
        ci_range <- abs(logIC50_ci[2] - logIC50_ci[1])
        if (ci_range > 0.666) {
          quality_flags <- c(quality_flags, "Wide logIC50 CI range")
        }
      }
      
      if (!is.null(plausibility_check) && isTRUE(plausibility_check$needs_correction)) {
        fit_label <- plausibility_check$fit_label
        if (is.null(fit_label)) fit_label <- "Implausible fit"
        # For flat curves the leading N/D label already covers this.
        if (!(identical(curve_type, "flat") && fit_label %in% c("Implausible fit", nd_labels))) {
          quality_flags <- c(quality_flags, fit_label)
        }
      }

      # Constrained-refit success: the IC50 is real but rests on the
      # constraint assumption -- say so explicitly.
      if (!is.null(plausibility_check) && !is.null(plausibility_check$constrained_label)) {
        quality_flags <- c(quality_flags, plausibility_check$constrained_label)
      }
      
      # Flag extreme Hill slopes, direction-aware
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
    ci_results <- calculate_ci(fit, prepared$df_clean)
    
    # If correct_parameter_order() flipped the HillSlope sign, the CI bounds
    # were computed on the original (wrong-sign) fit and must be negated.
    if (order_correction$was_corrected &&
        grepl("Hill Slope inconsistent", order_correction$correction_reason %||% "")) {
      hs_ci <- ci_results$HillSlope
      if (!any(is.na(hs_ci))) {
        ci_results$HillSlope <- sort(-hs_ci)  # negate and re-sort [lower, upper]
      }
    }
    
    gof_results <- calculate_goodness_of_fit(fit, prepared$df_clean)
    plausibility_check <- check_biological_plausibility(params, prepared$df_clean)
    # A converged but near-flat fit whose asymptote crossed a plausibility
    # limit is kept as-is and treated as no dose response (guard inside
    # check_biological_plausibility): this drives the flat quality label and
    # the N/D marking in batch_drc_analysis().
    if (isTRUE(plausibility_check$reclassify_flat)) curve_type <- "flat"

    # Constrained-refit fallback: when the unconstrained fit diverged (or its
    # LogIC50 left the tested range) but the data show a real edge response,
    # retry with constraints (bounds first, then Bottom fixed at 0) instead of
    # patching parameters post-hoc.  Flat/inactive compounds are NOT refit --
    # an inactive compound should not get an IC50.
    constrained_info <- NULL
    if (isTRUE(plausibility_check$needs_correction) &&
        isTRUE(plausibility_check$edge_response_detected) &&
        (isTRUE(plausibility_check$fit_diverged) ||
         isTRUE(plausibility_check$ic50_outside_range))) {
      cf <- try_constrained_fit(prepared$df_clean, curve_type)
      if (!is.null(cf)) {
        cf_params <- unname(stats::coef(cf$fit))
        # bottom_fixed fits have no Bottom coefficient; re-insert the fixed 0
        if (cf$method == "bottom_fixed") cf_params <- c(0, cf_params)
        cf_order <- correct_parameter_order(cf_params, prepared$df_clean, curve_type)
        if (cf_order$was_corrected) cf_params <- cf_order$corrected_params
        cf_check <- check_biological_plausibility(cf_params, prepared$df_clean)
        if (!isTRUE(cf_check$needs_correction)) {
          # Accept the constrained fit: rebuild all derived quantities
          fit <- cf$fit
          params <- cf_params
          order_correction <- cf_order
          initial_params <- c(params, 10^params[3], params[2] - params[1])
          ci_results <- calculate_ci(fit, prepared$df_clean)
          # Same HillSlope CI sign fix as for the unconstrained path
          if (cf_order$was_corrected &&
              grepl("Hill Slope inconsistent", cf_order$correction_reason %||% "")) {
            hs_ci <- ci_results$HillSlope
            if (!any(is.na(hs_ci))) {
              ci_results$HillSlope <- sort(-hs_ci)
            }
          }
          gof_results <- calculate_goodness_of_fit(fit, prepared$df_clean)
          plausibility_check <- cf_check
          constrained_info <- list(
            applied = TRUE,
            method = cf$method,
            label = if (cf$method == "bottom_fixed") {
              "IC50 from constrained fit (Bottom fixed at 0)"
            } else {
              "IC50 from constrained fit (plateaus bounded)"
            }
          )
          plausibility_check$constrained_label <- constrained_info$label
          if (verbose) {
            cat("Constrained refit accepted for ", comp_name, " (", cf$method, ")
", sep = "")
          }
        }
      }
    }

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
    
    curve_quality_info <- calculate_curve_quality(final_params, gof_results, plausibility_check, ci_results$LogIC50, curve_type)
    
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
      constrained_fit = constrained_info,
      data = prepared$df_clean,
      curve_type = curve_type
    )
  }
  
  # Normalize a data frame to 0-100% (same logic as fit_drc_3pl)
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
        LogIC50 = if (!apply_threshold) round(params[3], 3) else NA_real_,
        HillSlope = if (!apply_threshold) round(params[4], 3) else NA_real_,
        IC50 = if (!apply_threshold) format(params[5], scientific = TRUE) else NA_character_,
        Bottom_Lower_95CI = if (!is.na(ci$Bottom_Lower)) round(ci$Bottom_Lower, 3) else NA_real_,
        Bottom_Upper_95CI = if (!is.na(ci$Bottom_Upper)) round(ci$Bottom_Upper, 3) else NA_real_,
        Top_Lower_95CI = if (!is.na(ci$Top_Lower)) round(ci$Top_Lower, 3) else NA_real_,
        Top_Upper_95CI = if (!is.na(ci$Top_Upper)) round(ci$Top_Upper, 3) else NA_real_,
        LogIC50_Lower_95CI = if (!apply_threshold && !is.na(ci$LogIC50[1])) round(ci$LogIC50[1], 3) else NA_real_,
        LogIC50_Upper_95CI = if (!apply_threshold && !is.na(ci$LogIC50[2])) round(ci$LogIC50[2], 3) else NA_real_,
        HillSlope_Lower_95CI = if (!apply_threshold && !is.na(ci$HillSlope[1])) round(ci$HillSlope[1], 3) else NA_real_,
        HillSlope_Upper_95CI = if (!apply_threshold && !is.na(ci$HillSlope[2])) round(ci$HillSlope[2], 3) else NA_real_,
        IC50_Lower_95CI = if (!apply_threshold && !is.na(ci$IC50[1])) format(ci$IC50[1], scientific = TRUE) else NA_character_,
        IC50_Upper_95CI = if (!apply_threshold && !is.na(ci$IC50[2])) format(ci$IC50[2], scientific = TRUE) else NA_character_,
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
        Bottom = NA_real_, Top = NA_real_, LogIC50 = NA_real_, HillSlope = NA_real_, IC50 = NA_character_,
        Bottom_Lower_95CI = NA_real_, Bottom_Upper_95CI = NA_real_,
        Top_Lower_95CI = NA_real_, Top_Upper_95CI = NA_real_,
        LogIC50_Lower_95CI = NA_real_, LogIC50_Upper_95CI = NA_real_,
        HillSlope_Lower_95CI = NA_real_, HillSlope_Upper_95CI = NA_real_,
        IC50_Lower_95CI = NA_character_, IC50_Upper_95CI = NA_character_, Span = NA_real_,
        R_squared = NA_real_, Syx = NA_real_, Sum_of_Squares = NA_real_, Degrees_of_Freedom = NA_real_,
        Max_Slope = NA_real_, Curve_Quality = "Fit failed",
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
  
  # Return results object
  list(
    summary_table = summary_table,
    final_summary_table = final_summary_table,
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

