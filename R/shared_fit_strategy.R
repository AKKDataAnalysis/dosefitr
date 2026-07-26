#' Shared dose-response fitting strategy helpers
#'
#' These internal helpers centralise the three pieces of fitting *strategy* that
#' should behave identically across the package's two fitting engines:
#'
#' \enumerate{
#'   \item \code{detect_direction()} -- data-driven inhibition / activation / flat
#'         classification (scale-aware; works for BRET ratios and 0--100\% data).
#'   \item \code{make_start_grid()}  -- a grid of robust starting-value sets,
#'         mirroring the multi-start strategy in \code{fit_drc_4pl.R}.
#'   \item \code{choose_model()}     -- the shared 3PL-vs-4PL selection policy
#'         (rsdr-guarded, prefers the simpler 3PL unless 4PL is clearly better).
#' }
#'
#' Two parameterisations are supported so the same strategy can seed either engine:
#'   * \strong{drc}      : y = Bottom + (Top-Bottom)/(1 + 10^((log10x - LogIC50)*Hill))
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


# --- Direction detection ------------------------------------------------------

#' Detect dose-response direction from data (scale-aware)
#'
#' @param x_log Numeric vector of log-concentrations (log10 or ln; only the
#'   ordering is used).
#' @param y Numeric vector of responses (same length as \code{x_log}).
#' @param rel_frac Relative threshold as a fraction of the response range
#'   (default 0.15, matching \code{fit_drc_4pl.R}).
#' @param abs_floor Optional absolute floor on the threshold. If \code{NULL}
#'   (default) it is derived as \code{0.15 * range} with no fixed floor, which
#'   makes the test scale-invariant (safe for BRET ratios ~0.2--1 AND 0--100\%
#'   data). Pass e.g. \code{15} to reproduce \code{fit_drc_4pl.R}'s 0--100\%
#'   behaviour exactly.
#' @return One of \code{"inhibition"}, \code{"activation"}, \code{"flat"},
#'   or \code{"unknown"}.
#' @keywords internal
detect_direction <- function(x_log, y, rel_frac = 0.15, abs_floor = NULL) {
  ok <- is.finite(x_log) & is.finite(y)
  x_log <- x_log[ok]; y <- y[ok]
  if (length(y) < 4L) return("unknown")

  o <- order(x_log); ys <- y[o]
  k <- max(3L, floor(length(ys) / 4L))          # use up to a quarter of points per tail
  k <- min(k, floor(length(ys) / 2L))
  initial_avg <- mean(utils::head(ys, k), na.rm = TRUE)
  final_avg   <- mean(utils::tail(ys, k), na.rm = TRUE)

  range_val <- diff(range(ys, na.rm = TRUE))
  if (!is.finite(range_val) || range_val <= 0) return("flat")
  threshold <- rel_frac * range_val
  if (!is.null(abs_floor)) threshold <- max(abs_floor, threshold)

  if (initial_avg > final_avg + threshold) return("inhibition")
  if (final_avg > initial_avg + threshold) return("activation")
  "flat"
}


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


# --- Shared 3PL vs 4PL selection policy --------------------------------------

#' Shared 3PL-vs-4PL model-selection policy
#'
#' Prefers the simpler 3PL model and only accepts 4PL when it fits meaningfully
#' better, guarding against implausible Hill slopes. This is the single source of
#' truth for the selection rule used across engines.
#'
#' @param rsdr3,rsdr4 Robust (or ordinary) residual SD of the 3PL and 4PL fits.
#'   Either may be \code{NA}/\code{NULL} if that fit failed.
#' @param hill4 Estimated Hill slope (absolute value taken internally) of the 4PL fit.
#' @param improve_frac Fractional improvement in rsdr required to accept 4PL.
#'   4PL is accepted when \code{rsdr4 <= rsdr3 * (1 - improve_frac_tol)} where the
#'   default \code{tol} reproduces \code{rout_outliers()}'s \code{<= rsdr3 * 1.10}
#'   rule; see \code{tol10} argument.
#' @param hill_min,hill_max Plausibility bounds on \code{|hill4|}.
#' @param tol10 If \code{TRUE} (default) use \code{rout_outliers()}'s original
#'   \code{rsdr4 <= rsdr3 * 1.10} acceptance (i.e. accept 4PL unless it is >10\%
#'   WORSE). If \code{FALSE}, require genuine improvement:
#'   \code{rsdr4 <= rsdr3 * (1 - improve_frac)}.
#' @return List with \code{model} ("3PL" or "4PL") and \code{reason}.
#' @keywords internal
choose_model <- function(rsdr3, rsdr4, hill4,
                         improve_frac = 0.05,
                         hill_min = 0.1, hill_max = 5,
                         tol10 = TRUE) {
  has3 <- !is.null(rsdr3) && is.finite(rsdr3)
  has4 <- !is.null(rsdr4) && is.finite(rsdr4) && !is.null(hill4) && is.finite(hill4)

  if (!has4 && !has3) return(list(model = NA_character_, reason = "both fits failed"))
  if (!has4)          return(list(model = "3PL", reason = "4PL fit unavailable"))
  if (!has3)          return(list(model = "4PL", reason = "3PL fit unavailable"))

  mag_ok <- abs(hill4) >= hill_min && abs(hill4) <= hill_max
  if (!mag_ok)
    return(list(model = "3PL",
                reason = sprintf("4PL Hill=%.2f outside [%.1f,%.1f] -- rejected",
                                 hill4, hill_min, hill_max)))

  accept4 <- if (tol10) rsdr4 <= rsdr3 * 1.10
             else       rsdr4 <= rsdr3 * (1 - improve_frac)

  if (accept4)
    list(model = "4PL",
         reason = sprintf("4PL accepted (rsdr %.4f vs 3PL %.4f)", rsdr4, rsdr3))
  else
    list(model = "3PL",
         reason = sprintf("4PL not better (rsdr %.4f > 3PL %.4f) -- using 3PL", rsdr4, rsdr3))
}
