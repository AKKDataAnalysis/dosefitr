#' Shared dose-response fitting strategy helpers
#'
#' These internal helpers centralise the pieces of fitting *strategy* that
#' behave identically across the package's fitting engines:
#'
#' \enumerate{
#'   \item \code{make_start_grid()}        -- a grid of robust starting-value
#'         sets, mirroring the multi-start strategy in \code{fit_drc_4pl.R}.
#'   \item \code{optim_sign_for()}         -- Hill-slope sign convention helper.
#'   \item \code{detect_curve_type_shared()} -- data-driven inhibition /
#'         activation / flat classification (scale-aware and noise-aware; works
#'         for raw BRET ratios and 0--100\% data), shared by both fitters.
#' }
#'
#' Two parameterisations are supported so the same strategy can seed either engine:
#'   * \strong{drc}      : y = Bottom + (Top-Bottom)/(1 + 10^((LogIC50 - log10x)*Hill))
#'                         x on log10 scale, LogIC50 on log10 scale.
#'   * \strong{optim}    : OptimModel::hill_model, x linear, lec50 = ln(EC50),
#'                         slope m on the SAME numeric scale as drc's Hill
#'                         (verified: m == Hill).
#'
#' Sign convention note (IMPORTANT):
#'   * drc reports an INHIBITION curve with Hill < 0 (e.g. -1).
#'   * OptimModel::hill_model produces a DECREASING (inhibition) curve with m > 0
#'     (e.g. +1). So the two conventions are sign-FLIPPED. All conversions below
#'     handle this explicitly via \code{optim_sign_for()}.
#'
#' @keywords internal
#' @name shared_fit_strategy
NULL


# --- Sign helper --------------------------------------------------------------

#' Hill-slope sign for a given direction and parameterisation
#'
#' @param direction "inhibition" or "activation"/"agonist".
#' @param param "drc" or "optim".
#' @return +1 or -1 (numeric).
#' @keywords internal
optim_sign_for <- function(direction, param = c("optim", "drc")) {
  param <- match.arg(param)
  is_inhib <- direction %in% c("inhibition", "inhib")
  if (param == "drc") {
    if (is_inhib) -1 else 1        # drc: inhibition is negative Hill
  } else {
    if (is_inhib) 1 else -1        # OptimModel: inhibition is positive m
  }
}


# --- Plausibility-limit validation -------------------------------------------

#' Validate a user-supplied plausibility-limit argument
#'
#' Every configurable limit in \code{fit_drc_3pl()}, \code{fit_drc_4pl()} and
#' \code{batch_drc_analysis()} must be a numeric vector of length 2 with
#' \code{lower <= upper}. \code{-Inf}/\code{Inf} are allowed (open bounds);
#' \code{NA}/\code{NaN} are not. Fails fast with an error naming the argument.
#'
#' @param x The argument value supplied by the user.
#' @param name Character string: the argument name, used in the error message.
#' @return \code{x}, invisibly, when valid.
#' @keywords internal
.validate_limit_arg <- function(x, name) {
  if (!is.numeric(x) || length(x) != 2L || any(is.na(x)) || x[1] > x[2]) {
    stop(sprintf(
      "`%s` must be a numeric vector of length 2 with lower <= upper (got: %s).",
      name,
      paste(format(x), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(x)
}

#' Validate the full set of plausibility-limit arguments
#'
#' Convenience wrapper around \code{.validate_limit_arg()} for the ten limit
#' arguments shared by \code{fit_drc_3pl()}, \code{fit_drc_4pl()} and
#' \code{batch_drc_analysis()}. Values are read from \code{env} (the calling
#' function's environment) by name.
#'
#' @param env An environment containing the ten limit arguments.
#' @return Invisibly \code{TRUE} when all are valid.
#' @keywords internal
.validate_all_limit_args <- function(env) {
  limit_arg_names <- c(
    "bottom_limits_nanobret_inhibition", "bottom_limits_nanobret_activation",
    "top_limits_nanobret_inhibition",    "top_limits_nanobret_activation",
    "bottom_limits_viability_inhibition", "bottom_limits_viability_activation",
    "top_limits_viability_inhibition",   "top_limits_viability_activation",
    "logIC50_limits", "hill_slope_limits"
  )
  for (nm in limit_arg_names) .validate_limit_arg(get(nm, envir = env), nm)
  invisible(TRUE)
}


# --- Starting-value grid ------------------------------------------------------

#' Build a grid of robust starting values (mirrors fit_drc_4pl.R strategy)
#'
#' @param x_log Numeric vector of log-concentrations.
#' @param y Numeric vector of responses.
#' @param direction "inhibition" or "activation"/"agonist".
#' @param param "drc" or "optim" -- controls the parameter names/scales returned.
#' @param slopes Numeric vector of |Hill| magnitudes to seed. Default
#'   \code{c(1, 0.5, 2, 3, 4)} spans shallow to steep so that both gentle and
#'   steep curves get an adequate seed (the original \code{fit_drc_4pl.R} grid
#'   maxed at 2, which under-seeds steep NanoBRET/viability curves).
#' @return A list of named numeric vectors. For \code{param="drc"} names are
#'   \code{Bottom,Top,LogIC50,HillSlope}; for \code{param="optim"} names are
#'   \code{emin,emax,lec50,m} (lec50 = ln(EC50), m on the same scale as Hill).
#' @keywords internal
make_start_grid <- function(x_log, y, direction,
                            param = c("drc", "optim"),
                            slopes = c(1, 0.5, 2, 3, 4)) {
  param <- match.arg(param)
  ok <- is.finite(x_log) & is.finite(y)
  x_log <- x_log[ok]; y <- y[ok]

  min_r <- min(y, na.rm = TRUE); max_r <- max(y, na.rm = TRUE)
  mid   <- (min_r + max_r) / 2

  # Approximate the log-concentration at the half-maximal response (log10 or ln,
  # whatever x_log is on). Robust to non-monotone noise via approx().
  approx_ic50_log <- tryCatch(
    suppressWarnings(stats::approx(y, x_log, xout = mid, ties = mean)$y),
    error = function(e) stats::median(x_log, na.rm = TRUE))
  if (!is.finite(approx_ic50_log))
    approx_ic50_log <- stats::median(x_log, na.rm = TRUE)

  sgn <- optim_sign_for(direction, param)

  # Plateau perturbations mirror fit_drc_4pl.R's five strategies.
  plateau_sets <- list(
    c(min_r,       max_r),
    c(min_r * 0.8, max_r * 1.2),
    c(min_r * 1.2, max_r * 0.8),
    c(min_r,       max_r),
    c(min_r,       max_r)
  )

  n <- max(length(slopes), length(plateau_sets))
  slopes       <- rep(slopes,       length.out = n)
  plateau_sets <- rep(plateau_sets, length.out = n)

  lapply(seq_len(n), function(i) {
    bot <- plateau_sets[[i]][1]; top <- plateau_sets[[i]][2]
    hill <- sgn * slopes[i]
    if (param == "drc") {
      # LogIC50 must be on the SAME log scale as x_log passed in.
      c(Bottom = bot, Top = top, LogIC50 = approx_ic50_log, HillSlope = hill)
    } else {
      # OptimModel needs lec50 = ln(EC50). If x_log is log10, convert; if it is
      # already ln, use directly. Caller must tell us via attribute on x_log.
      lec50_val <- if (isTRUE(attr(x_log, "is_ln"))) approx_ic50_log
                   else approx_ic50_log * log(10)      # log10 -> ln
      c(emin = bot, emax = top, lec50 = lec50_val, m = hill)
    }
  })
}


# --- Shared curve-direction classifier ---------------------------------------

#' Classify a dose-response curve as inhibition, activation, or flat
#'
#' Single source of truth for curve-direction detection, shared by
#' \code{fit_drc_3pl} and \code{fit_drc_4pl}. Compares the mean of the three
#' lowest-concentration responses with the three highest, and calls the curve
#' non-flat only when that difference exceeds a threshold that is BOTH
#' scale-relative and noise-aware:
#'   (a) 15%% of the observed response range (works at any scale, including
#'       raw BRET ratios whose full range is only a fraction of a unit), and
#'   (b) 3 * SE(diff), where SE(diff) = sigma / sqrt(3) and sigma is the
#'       trend-robust lag-1 (von Neumann) noise estimate
#'       sigma = sqrt(sum(diff(resps)^2) / (2 * (n - 1))). This is a ~3-sigma
#'       signal-to-noise test that suppresses random-noise false positives.
#' threshold = max((a),(b)). No absolute floor (an absolute floor of 15, as in
#' the previous fit_drc_4pl copy, mislabels real raw-scale curves as flat).
#'
#' The von Neumann estimator is used instead of the raw variance of the first
#' and last three points because those extreme points sit on the shoulders of
#' a sigmoid: their within-window variance is dominated by curve SLOPE, not
#' measurement noise, so \code{3 * sqrt((var(head3)+var(tail3))/3)} explodes on
#' clean-but-steep curves and mislabels genuine inhibition as flat. Successive
#' first differences are nearly insensitive to a smooth monotonic trend, so the
#' lag-1 estimate isolates noise regardless of curve steepness or scale.
#'
#' Rows with a missing (\code{NA}) \code{log_inhibitor} are dropped before
#' ordering: a point with no concentration cannot be assigned to the low- or
#' high-dose end, and (e.g. control rows carrying 0/100 sentinels) would
#' otherwise sort to one end and corrupt the head/tail statistics.
#'
#' Replicates are collapsed to a per-concentration MEAN profile before the
#' head/tail comparison and the noise estimate. Direction is a property of the
#' concentration-response profile, not of the individual wells. When two
#' replicates disagree by a large systematic offset, sorting the pooled long
#' vector by concentration interleaves them (high, low, high, low, ...), which
#' both scrambles the first/last-three means and inflates the lag-1 noise
#' estimate with between-replicate offset, mislabelling a genuine curve as flat.
#' Aggregating first removes that artifact.
#'
#' @param data A data frame with columns \code{log_inhibitor} and
#'   \code{response} (extra columns ignored).
#' @return One of "inhibition", "activation", "flat", or "unknown".
#' @keywords internal
detect_curve_type_shared <- function(data) {
  if (is.null(data) || nrow(data) < 4) return("unknown")
  # Keep only rows with a real concentration AND a real response. A point with
  # no concentration cannot be ordered onto the low/high-dose end.
  keep <- !is.na(data$log_inhibitor) & !is.na(data$response)
  clean_df <- data[keep, , drop = FALSE]
  if (nrow(clean_df) < 4) return("unknown")

  # Collapse replicates to a per-concentration mean profile, sorted by dose.
  agg <- stats::aggregate(response ~ log_inhibitor, data = clean_df, FUN = mean)
  agg <- agg[order(agg$log_inhibitor), , drop = FALSE]
  resps <- agg$response
  n <- length(resps)
  if (n < 4) return("unknown")

  initial_avg <- mean(utils::head(resps, 3))
  final_avg   <- mean(utils::tail(resps, 3))

  range_val <- diff(range(resps))
  thr_rel   <- range_val * 0.15

  # Trend-robust noise estimate from successive (lag-1) differences.
  d <- diff(resps)
  sigma <- if (n >= 2) sqrt(sum(d^2) / (2 * (n - 1))) else 0
  if (!is.finite(sigma)) sigma <- 0
  se_diff   <- sigma / sqrt(3)
  thr_noise <- if (se_diff > 0) 3 * se_diff else 0

  threshold <- max(thr_rel, thr_noise)

  if (initial_avg > final_avg + threshold) return("inhibition")
  if (final_avg > initial_avg + threshold) return("activation")
  "flat"
}

# ============================================================================
# Adaptive number formatting for report tables
# ============================================================================

# Plain notation for values in [1e-3, 1e7), scientific only outside;
# always 4 significant digits. Used for IC50/EC50 and its CI columns.
fmt_adaptive <- function(x, d = 4) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(NA_character_)
  if (abs(x) >= 1e-3 && abs(x) < 1e7) {
    format(signif(x, d), scientific = FALSE, trim = TRUE)
  } else {
    format(signif(x, d), scientific = TRUE, trim = TRUE)
  }
}

# Apply a plain (non-scientific) Excel number format to the numeric columns
# of a data frame already written to a sheet. Excel's default "General"
# format auto-switches to scientific notation for very small/large values.
add_plain_numfmt <- function(wb, sheet, df, numFmt = "0.############") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(invisible(wb))
  num_cols <- which(vapply(df, is.numeric, logical(1)))
  if (length(num_cols) == 0) return(invisible(wb))
  st <- openxlsx::createStyle(numFmt = numFmt)
  openxlsx::addStyle(wb, sheet, style = st, rows = 1 + seq_len(nrow(df)),
                     cols = num_cols, gridExpand = TRUE, stack = TRUE)
  invisible(wb)
}
