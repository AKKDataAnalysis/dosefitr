#' ROUT-based outlier detection for dose-response curves
#'
#' Detects outliers in dose-response data using the ROUT (Robust regression and
#' Outlier removal) method applied to 3-parameter (3PL) or 4-parameter (4PL)
#' Hill models. Supports replicate-aware filtering and safeguards against
#' model misfit.
#'
#' @param data A data.frame where the first column contains log-transformed
#' concentrations and the remaining columns contain response values (e.g. BRET ratios).
#'
#' @param Q False discovery rate (FDR) threshold for ROUT outlier detection.
#' Must be between 0 and 1 (exclusive). Default is 0.01.
#'
#' @param n_param Number of model parameters: 3 (fixed Hill slope) or 4 (free Hill slope).
#'
#' @param conc_col Index of the concentration column in `data`. Default is 1.
#'
#' @param log_base Logarithm base of the concentration values. Either `"log10"` or `"ln"`.
#'
#' @param direction Direction of the dose-response curve:
#' `"inhibition"` or `"agonist"`. Determines Hill slope sign.
#'
#' @param min_dynamic_range Minimum dynamic range (%) required for reliable fitting.
#'
#' @param ntry_retry Number of retry attempts for model fitting using random starts.
#'
#' @param verbose Logical; if TRUE, prints detailed progress and diagnostics.
#'
#' @details
#' The function performs the following steps:
#'
#' \enumerate{
#'   \item Removes control rows (NA concentrations) from fitting.
#'   \item Fits 3PL and/or 4PL Hill models using robust regression.
#'   \item Applies ROUT outlier detection based on standardized residuals.
#'   \item Applies replicate-consistency filtering.
#'   \item Removes systematic residual patterns (model misfit safeguard).
#'   \item Replaces detected outliers with NA in the cleaned dataset.
#' }
#'
#' Additional safeguards include:
#' \itemize{
#'   \item Hill slope sign and magnitude constraints (4PL only)
#'   \item Dynamic range warnings
#'   \item Minimum data point requirements
#' }
#'
#' @return A list with the following elements:
#'
#' \describe{
#'   \item{results}{Data frame with fitted values, residuals, and outlier flags.}
#'   \item{cleaned_table}{Input data with detected outliers replaced by NA.}
#'   \item{outlier_table}{Summary of detected outliers.}
#'   \item{skipped_table}{Compounds excluded from analysis with reasons.}
#'   \item{cleared_systematic}{Flags removed due to systematic residual patterns.}
#'   \item{params}{List of parameters used in the analysis.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- rout_outliers(data = my_data)
#'
#' result$outlier_table
#' result$cleaned_table
#' }
#'
#' @references
#' Motulsky HJ, Brown RE (2006).
#' Detecting outliers when fitting data with nonlinear regression.
#' BMC Bioinformatics.
#'
#' @seealso
#' \code{\link{rout_outliers_batch}}
#'
#' @export

rout_outliers <- function(data,
                          Q                  = 0.01,
                          n_param            = 4L,
                          conc_col           = 1L,
                          log_base           = "log10",
                          direction          = "inhibition",
                          min_dynamic_range  = 20,
                          ntry_retry         = 3L,
                          verbose            = TRUE) {


  # Internal helpers (not exported)
  # ---------------------------------------------------------------------------

  # .robust_starts: compute robust starting values for the Hill model.
  # Uses the OptimModel hill_model start function, then maps the four
  # parameters (emin, emax, lec50, m) to the (b0, t0, lec50, hill) names
  # used by the rest of rout_outliers().
  #
  # Returns a list with:
  #   $b0     --  bottom (emin): response at x -> 0 (vehicle / 0% control)
  #   $t0     --  top   (emax): response at x -> Inf (saturating compound)
  #   $lec50  --  log(EC50) on the natural-log scale (as used by hill_model)
  #   $hill   --  Hill slope starting value (sign set by hill_fixed)
  .robust_starts <- function(y, x, hill_fixed) {
    sv <- tryCatch(
      attr(OptimModel::hill_model, "start")(x, y),
      error = function(e) {
        c(emin  = min(y, na.rm = TRUE),
          emax  = max(y, na.rm = TRUE),
          lec50 = log(median(x[x > 0], na.rm = TRUE)),
          m     = hill_fixed)
      }
    )
    list(
      b0    = sv[["emin"]],
      t0    = sv[["emax"]],
      lec50 = sv[["lec50"]],
      hill  = hill_fixed
    )
  }

  # .fit_model: fit a 3PL or 4PL Hill model via rout_fitter().
  # Handles the sign/magnitude guard for 4PL and the retry logic.
  #
  # rout_fitter() argument order: rout_fitter(theta0, f.model, x, y, ...)
  #
  # Returns a list with:
  #   $fit      --  rout_fit object, or NULL if fitting failed / guard triggered
  #   $n_param  --  number of parameters actually used (3 or 4)
  #   $f_model  --  the model function used (for curve prediction)
  #   $retried  --  TRUE if a random-restart retry was accepted
  .fit_model <- function(x, y, n_param, hill_fixed, Q, ntry_retry) {
    # Build a 3PL model by fixing the Hill slope to hill_fixed.
    # hill_model uses 4 parameters: (emin, emax, lec50, m).
    # For 3PL we fix m = hill_fixed and optimise only (emin, emax, lec50).
    if (n_param == 3L) {
      hf <- hill_fixed  # capture in closure
      f_model <- function(theta, x) {
        OptimModel::hill_model(c(theta[1L], theta[2L], theta[3L], hf), x)
      }
      attr(f_model, "start") <- function(x, y) {
        sv <- attr(OptimModel::hill_model, "start")(x, y)
        sv[1:3]  # drop the Hill slope
      }
    } else {
      f_model <- OptimModel::hill_model
    }

    fit <- tryCatch(
      OptimModel::rout_fitter(theta0 = NULL, f.model = f_model,
                              x = x, y = y, Q = Q, ntry = 0L),
      error = function(e) NULL
    )

    retried <- FALSE

    # Retry logic: attempt random restarts if initial fit did not converge
    if (!is.null(fit) && !fit$Converge && ntry_retry > 0L) {
      fit_retry <- tryCatch(
        OptimModel::rout_fitter(theta0 = NULL, f.model = f_model,
                                x = x, y = y, Q = Q, ntry = ntry_retry),
        error = function(e) NULL
      )
      # Accept retry only if it converged AND has rsdr <= original rsdr
      if (!is.null(fit_retry) && fit_retry$Converge &&
          fit_retry$rsdr <= fit$rsdr) {
        fit     <- fit_retry
        retried <- TRUE
      }
    }

    # 4PL guard: reject if Hill slope has wrong sign or extreme magnitude
    if (n_param == 4L && !is.null(fit)) {
      hill_est <- fit$par[4L]
      sign_ok  <- (hill_fixed * hill_est) > 0
      mag_ok   <- abs(hill_est) >= 0.1 && abs(hill_est) <= 5
      if (!sign_ok || !mag_ok) return(list(fit=NULL, n_param=n_param,
                                           f_model=f_model, retried=FALSE))
    }

    list(fit = fit, n_param = n_param, f_model = f_model, retried = retried)
  }

  # .clean_flags: ensure outlier.adj is a plain logical vector.
  .clean_flags <- function(fit) {
    fit$outlier.adj <- as.logical(unname(fit$outlier.adj))
    fit$outlier     <- as.logical(unname(fit$outlier))
    fit
  }

  # .detect_systematic_flags: post-hoc check for model-misfit false positives.
  #
  # Minimum 2 flagged points required: a single flagged point cannot be "consecutive"
  # and is not a systematic pattern by definition.
  #
  # Returns a list:
  #   $flags    -- updated logical vector (same length as outlier_flags input)
  #   $cleared  -- TRUE if flags were cleared by this filter, FALSE otherwise
  #   $reason   -- character string describing why flags were cleared (or "")
  .detect_systematic_flags <- function(outlier_flags, std_residuals, x_log_fit) {
    flagged_idx <- which(outlier_flags)

    # Need at least 2 flagged points AND at least 2 unique concentrations.
    # (Two replicates flagged at the same concentration is a true outlier pair,
    # not a systematic model-misfit pattern.)
    if (length(flagged_idx) < 2L)
      return(list(flags = outlier_flags, cleared = FALSE, reason = ""))

    # Condition 1: all flagged residuals have the same sign
    flag_signs <- sign(std_residuals[flagged_idx])
    if (!all(flag_signs == flag_signs[1L]))
      return(list(flags = outlier_flags, cleared = FALSE, reason = ""))

    # Condition 2: flagged points span >= 2 consecutive concentrations.
    # "Consecutive" = adjacent in the sorted unique concentration vector,
    # allowing a gap of 1 (one missing step between two flagged steps).
    flagged_concs  <- sort(unique(x_log_fit[flagged_idx]))
    if (length(flagged_concs) < 2L)
      return(list(flags = outlier_flags, cleared = FALSE, reason = ""))
    all_concs      <- sort(unique(x_log_fit))
    conc_positions <- match(flagged_concs, all_concs)  # positions in full conc ladder

    # Check that the flagged positions form a contiguous run (max gap = 1 step)
    pos_diffs <- diff(conc_positions)
    if (any(pos_diffs > 2L))
      return(list(flags = outlier_flags, cleared = FALSE, reason = ""))

    # Both conditions met: clear all flags
    sign_word <- if (flag_signs[1L] > 0) "positive" else "negative"
    reason <- sprintf(
      "systematic %s residuals at %d consecutive concentration(s) [%.2f to %.2f] -- possible model misfit",
      sign_word, length(flagged_concs),
      min(flagged_concs), max(flagged_concs))

    outlier_flags[flagged_idx] <- FALSE
    list(flags         = outlier_flags,
         cleared       = TRUE,
         n_cleared     = length(flagged_idx),
         residual_sign = sign_word,
         reason        = reason)
  }


  # ---- Input validation ----
  if (!n_param %in% c(3L, 4L))
    stop("n_param must be 3 or 4")
  if (!direction %in% c("inhibition", "agonist"))
    stop('direction must be "inhibition" or "agonist"')
  if (!log_base %in% c("log10", "ln"))
    stop('log_base must be "log10" or "ln"')
  if (!is.numeric(ntry_retry) || length(ntry_retry) != 1 || ntry_retry < 0)
    stop("ntry_retry must be a single non-negative integer")
  ntry_retry <- as.integer(ntry_retry)
  if (!is.numeric(Q) || length(Q) != 1 || Q <= 0 || Q >= 1)
    stop("Q must be a single number between 0 and 1 (exclusive)")

  conc_raw      <- data[[conc_col]]
  ctrl_rows     <- which(is.na(conc_raw))
  dose_rows     <- which(!is.na(conc_raw))
  x_log         <- as.numeric(conc_raw[dose_rows])
  x_raw         <- if (log_base == "log10") 10^x_log else exp(x_log)
  hill_fixed    <- if (direction == "inhibition") 1L else -1L
  # Column name for the log-concentration output column  --  reflects actual log base used
  conc_col_name <- if (log_base == "log10") "log10_conc" else "ln_conc"

  if (verbose) {
    model_desc <- if (n_param == 4L) "4PL (Hill free, sign-constrained)" else "3PL (Hill fixed)"
    cat(sprintf("Control rows     : %s (excluded from outlier detection)\n",
                if (length(ctrl_rows) > 0) paste(ctrl_rows, collapse = ", ") else "none"))
    cat(sprintf("Dose rows used   : %d\n", length(dose_rows)))
    cat(sprintf("Default model    : %s\n", model_desc))
    cat(sprintf("ROUT Q           : %.3f\n\n", Q))
  }

  value_cols <- setdiff(seq_len(ncol(data)), conc_col)
  col_names  <- colnames(data)[value_cols]
  valid_mask <- !grepl("^NA", col_names)
  value_cols <- value_cols[valid_mask]
  col_names  <- col_names[valid_mask]
  # Strip trailing ".N" suffix to group replicates (e.g. "SRP9" and "SRP9.2" -> "SRP9").
  # Limitation: compound names that genuinely end in ".N" (e.g. "Compound.1") will be
  # incorrectly merged with "Compound". Avoid such names in the input data.
  base_names       <- sub("\\.(\\d+)$", "", col_names)
  unique_compounds <- unique(base_names)

  if (verbose) cat(sprintf("Compounds found  : %d\n\n", length(unique_compounds)))

  # collect skipped entries via the lapply return value instead of
  # <<- superassignment. Each element is either a results data.frame (normal),
  # NULL (fit failed  --  already handled below), or a list(skipped=<df>) sentinel.
  # This is safe under parallelisation and avoids shared mutable state.

  results_list <- lapply(unique_compounds, function(cmpd) {

    rep_cols   <- value_cols[base_names == cmpd]
    # Use [[col]][rows] not [rows, col]  --  the latter returns a tibble when data is
    # a tibble (e.g. from readxl), which cannot be coerced to numeric directly.
    y_list     <- lapply(rep_cols, function(col) as.numeric(data[[col]][dose_rows]))
    x_rep      <- rep(x_raw,  length(rep_cols))
    x_log_rep  <- rep(x_log,  length(rep_cols))
    y_rep      <- unlist(y_list)
    rep_labels <- rep(seq_along(rep_cols), each = length(dose_rows))
    row_idx    <- rep(dose_rows, length(rep_cols))
    col_idx    <- rep(rep_cols,  each = length(dose_rows))

    valid     <- !is.na(x_rep) & !is.na(y_rep)
    x_fit     <- x_rep[valid];     x_log_fit <- x_log_rep[valid]
    y_fit     <- y_rep[valid];     rep_fit   <- rep_labels[valid]
    row_fit   <- row_idx[valid];   col_fit   <- col_idx[valid]

    # Minimum 8 points: ensures ~2 observations per parameter for 4PL
    if (length(x_fit) < 8L) {
      return(list(skipped = data.frame(
        compound = cmpd,
        reason   = sprintf("too few valid points (%d, minimum 8)", length(x_fit)),
        n_valid  = length(x_fit),
        dynamic_range_pct = NA_real_,
        stringsAsFactors = FALSE)))
    }

    # ---- Dynamic range (robust, always non-negative) ----
    sv        <- .robust_starts(y_fit, x_fit, hill_fixed)
    dyn_range <- if (abs(sv$b0) < 1e-9) NA_real_
    else abs(100 * (sv$b0 - sv$t0) / sv$b0)

    if (!is.na(dyn_range) && dyn_range < min_dynamic_range && verbose)
      message(sprintf(
        "Warning: %s: low dynamic range (%.1f%% < %.0f%%)  --  curve may be flat or non-responsive",
        cmpd, dyn_range, min_dynamic_range))

    # ---- Fit 3PL (always) and 4PL (when n_param=4) in a single pass ----
    #
    # Both models are fitted upfront so the rsdr comparison requires no extra
    # optimizer call. The decision logic then selects the winner cleanly.
    res3 <- .fit_model(x_fit, y_fit, n_param = 3L, hill_fixed = hill_fixed, Q = Q,
                       ntry_retry = ntry_retry)

    if (n_param == 4L) {
      res4 <- .fit_model(x_fit, y_fit, n_param = 4L, hill_fixed = hill_fixed, Q = Q,
                         ntry_retry = ntry_retry)

      use_4pl <- !is.null(res4$fit) &&
        !is.null(res3$fit) &&
        res4$fit$rsdr <= res3$fit$rsdr * 1.10

      if (use_4pl) {
        res_chosen <- res4
        if (verbose) message(sprintf(
          "%s: using 4PL (Hill=%.3f, rsdr %.1f%% better than 3PL)",
          cmpd, res4$fit$par[4L],
          100 * (1 - res4$fit$rsdr / res3$fit$rsdr)))
      } else {
        res_chosen <- res3
        if (!is.null(res4$fit) && verbose) {
          # 4PL fitted but rsdr guard rejected it
          message(sprintf(
            "%s: 4PL rsdr (%.4f) > 3PL rsdr (%.4f) x 1.10  --  using 3PL",
            cmpd, res4$fit$rsdr, res3$fit$rsdr))
        } else if (is.null(res4$fit) && verbose) {
          message(sprintf(
            "%s: 4PL sign/magnitude guard triggered  --  using 3PL", cmpd))
        }
      }
    } else {
      res_chosen <- res3
    }

    if (is.null(res_chosen$fit)) {
      return(list(skipped = data.frame(
        compound = cmpd, reason = "fit failed (rout_fitter error)",
        n_valid = length(x_fit),
        dynamic_range_pct = if (is.na(dyn_range)) NA_real_ else round(dyn_range, 1),
        stringsAsFactors = FALSE)))
    }

    fit        <- .clean_flags(res_chosen$fit)
    used_model <- if (res_chosen$n_param == 4L) "4PL" else "3PL"
    was_retried <- isTRUE(res_chosen$retried)

    # Verbose retry notification  --  distinct from the final convergence summary
    if (was_retried && verbose) {
      if (fit$Converge) {
        message(sprintf("%s: converged after random-restart retry", cmpd))
      } else {
        message(sprintf("%s: did not converge even after retry  --  outlier calls may be unreliable", cmpd))
      }
    }

    # ---- Replicate-agreement filter ----
    # Clears a flag if the paired replicate at the same concentration is NOT
    # flagged AND |flagged - paired| < 2 * rsdr (both reps agree within noise).
    if (any(fit$outlier.adj)) {
      rsdr_val      <- fit$rsdr
      outlier_flags <- fit$outlier.adj
      # Works correctly for >2 replicates:
      #   - any(outlier_flags[paired_idx]): if ANY other rep at this conc is also
      #     flagged, we keep the flag (both reps are suspicious  --  likely a real effect).
      #   - all(abs(...)): ALL other reps must agree within 2*rsdr for the flag to
      #     be cleared. With triplicates this is a stricter, more conservative test.
      for (i in which(outlier_flags)) {
        conc_i     <- x_log_fit[i]
        same_conc  <- (x_log_fit == conc_i) & (seq_along(x_log_fit) != i)
        paired_idx <- which(same_conc)
        if (length(paired_idx) == 0L) next          # singleton at this conc  --  keep flag
        if (any(outlier_flags[paired_idx])) next     # another rep also flagged  --  keep flag
        if (all(abs(y_fit[i] - y_fit[paired_idx]) < 2 * rsdr_val))
          outlier_flags[i] <- FALSE                  # all other reps agree  --  clear flag
      }
      fit$outlier.adj <- outlier_flags
    }

    # ---- Systematic-residual filter ----
    # Runs AFTER the replicate-agreement filter on the surviving flags.
    # Clears flags that form a systematic same-sign consecutive pattern,
    # which indicates model misfit rather than true outliers.
    sys_result <- NULL
    if (any(fit$outlier.adj)) {
      sys_result      <- .detect_systematic_flags(fit$outlier.adj,
                                                  fit$sresiduals,
                                                  x_log_fit)
      fit$outlier.adj <- sys_result$flags
      if (sys_result$cleared && verbose)
        message(sprintf("%s: systematic residual filter cleared %d flag(s) -- %s",
                        cmpd, sys_result$n_cleared, sys_result$reason))
    }

    # ---- Extract curve parameters ----
    n_used  <- res_chosen$n_param
    f_model <- res_chosen$f_model

    df_out <- data.frame(
      compound          = cmpd,
      replicate         = rep_fit,
      bret_ratio        = y_fit,
      fitted            = f_model(fit$par[seq_len(n_used)], x_fit),
      residual          = fit$residuals,
      rsdr              = fit$rsdr,
      std_residual      = fit$sresiduals,
      outlier_raw       = fit$outlier,
      outlier_fdr       = fit$outlier.adj,
      converged         = fit$Converge,
      retried           = was_retried,
      model_used        = used_model,
      bottom            = fit$par[1L],
      top               = fit$par[2L],
      log10_EC50        = log10(exp(fit$par[3L])),  # par[3] is ln(EC50); convert to log10 for output
      hill_slope        = if (n_used == 4L) fit$par[4L] else hill_fixed,
      dynamic_range_pct = round(dyn_range, 1),
      .row_idx          = row_fit,
      .col_idx          = col_fit,
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
    # Insert log-concentration column at position 3 with the correct name
    # ("log10_conc" for log_base="log10", "ln_conc" for log_base="ln")
    df_out <- cbind(df_out[, 1:2, drop = FALSE],
                    setNames(data.frame(x_log_fit), conc_col_name),
                    df_out[, 3:ncol(df_out), drop = FALSE])

    # Return a list carrying both the results data.frame and any systematic-
    # filter record. The post-processing block below separates these.
    list(results = df_out,
         cleared = if (!is.null(sys_result) && sys_result$cleared)
           data.frame(compound        = cmpd,
                      n_flags_cleared = sys_result$n_cleared,
                      residual_sign   = sys_result$residual_sign,
                      reason          = sys_result$reason,
                      stringsAsFactors = FALSE)
         else NULL)
  })

  # DESIGN_1 fix (continued): separate skipped sentinels from real results.
  skipped_list <- Filter(Negate(is.null),
                         lapply(results_list, function(x)
                           if (is.list(x) && !is.null(x$skipped)) x$skipped else NULL))
  results_raw  <- Filter(Negate(is.null),
                         lapply(results_list, function(x)
                           if (is.list(x) && !is.null(x$results)) x$results else NULL))
  cleared_raw  <- Filter(Negate(is.null),
                         lapply(results_list, function(x)
                           if (is.list(x) && !is.null(x$cleared)) x$cleared else NULL))

  # Guard: if every compound was skipped, return an empty results data.frame
  # rather than NULL (which would crash nrow() and subsetting downstream).
  if (length(results_raw) == 0L) {
    results <- data.frame(
      compound = character(), replicate = integer(), bret_ratio = numeric(),
      fitted = numeric(), residual = numeric(), rsdr = numeric(),
      std_residual = numeric(), outlier_raw = logical(), outlier_fdr = logical(),
      converged = logical(), retried = logical(), model_used = character(),
      bottom = numeric(), top = numeric(), log10_EC50 = numeric(),
      hill_slope = numeric(), dynamic_range_pct = numeric(),
      stringsAsFactors = FALSE)
    results[[conc_col_name]] <- numeric(0)
  } else {
    results <- do.call(rbind, results_raw)
  }

  # ---- Build cleaned table (outliers replaced with NA) ----
  cleaned_table <- data
  outlier_rows  <- if (nrow(results) > 0L) results[results$outlier_fdr, ] else results[0L, ]

  replaced_log <- if (nrow(outlier_rows) > 0L) {
    rl <- data.frame(
      compound       = outlier_rows$compound,
      column         = colnames(data)[outlier_rows$.col_idx],
      row            = outlier_rows$.row_idx,
      original_value = round(outlier_rows$bret_ratio, 4),
      stringsAsFactors = FALSE,
      row.names      = NULL
    )
    # Insert concentration column with the correct name at position 4
    rl <- cbind(rl[, 1:3, drop = FALSE],
                setNames(data.frame(outlier_rows[[conc_col_name]]), conc_col_name),
                rl[, 4, drop = FALSE])
    rl
  } else {
    rl <- data.frame(compound=character(), column=character(), row=integer(),
                     original_value=numeric())
    rl[[conc_col_name]] <- numeric(0)
    rl
  }

  for (i in seq_len(nrow(outlier_rows)))
    cleaned_table[outlier_rows$.row_idx[i], outlier_rows$.col_idx[i]] <- NA

  # Attach replacement log as an attribute for self-documentation
  attr(cleaned_table, "outliers_replaced") <- replaced_log

  # ---- Outlier summary table ----
  outlier_table <- if (nrow(outlier_rows) > 0L) {
    ot <- data.frame(
      compound          = outlier_rows$compound,
      column            = colnames(data)[outlier_rows$.col_idx],
      replicate         = outlier_rows$replicate,
      row               = outlier_rows$.row_idx,
      bret_ratio        = round(outlier_rows$bret_ratio, 4),
      fitted            = round(outlier_rows$fitted, 4),
      std_residual      = round(outlier_rows$std_residual, 3),
      dynamic_range_pct = outlier_rows$dynamic_range_pct,
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
    # Insert concentration column with the correct name at position 5
    ot <- cbind(ot[, 1:4, drop = FALSE],
                setNames(data.frame(outlier_rows[[conc_col_name]]), conc_col_name),
                ot[, 5:ncol(ot), drop = FALSE])
    ot
  } else {
    empty_ot <- data.frame(
      compound          = character(),
      column            = character(),
      replicate         = integer(),
      row               = integer(),
      stringsAsFactors  = FALSE)
    empty_ot[[conc_col_name]] <- numeric(0)   # position 5  --  matches non-empty cbind layout
    empty_ot$bret_ratio        <- numeric(0)
    empty_ot$fitted            <- numeric(0)
    empty_ot$std_residual      <- numeric(0)
    empty_ot$dynamic_range_pct <- numeric(0)
    empty_ot
  }

  # ---- Skipped compounds table ----
  skipped_table <- if (length(skipped_list) > 0L) {
    do.call(rbind, skipped_list)
  } else {
    data.frame(compound=character(), reason=character(),
               n_valid=integer(), dynamic_range_pct=numeric())
  }

  # ---- Cleared-systematic table ----
  # Audit trail of flags cleared by the systematic-residual filter.
  # Always present (empty data frame when nothing was cleared).
  cleared_systematic <- if (length(cleared_raw) > 0L) {
    do.call(rbind, cleared_raw)
  } else {
    data.frame(compound        = character(),
               n_flags_cleared = integer(),
               residual_sign   = character(),
               reason          = character(),
               stringsAsFactors = FALSE)
  }
  rownames(cleared_systematic) <- NULL

  results$.row_idx <- NULL
  results$.col_idx <- NULL

  # ---- Verbose summary ----
  if (verbose) {
    # Convergence summary  --  distinguish three outcomes:
    #   (1) Converged on first try        --  no message (expected, silent)
    #   (2) Converged after retry         --  already reported per-compound above
    #   (3) Did not converge after retry  --  reported here as a batch warning
    if (nrow(results) > 0L) {
      non_conv <- unique(results$compound[!results$converged])
      if (length(non_conv) > 0L) {
        retry_attempted <- ntry_retry > 0L
        message(sprintf(
          "Warning: optimizer did not converge%s for: %s\n  Outlier calls for these compounds may be unreliable.",
          if (retry_attempted) " (even after retry)" else "",
          paste(non_conv, collapse = ", ")))
      }
    }

    if (nrow(skipped_table) > 0L) {
      cat("Skipped compounds:\n")
      print(skipped_table, row.names = FALSE)
      cat("\n")
    }

    if (nrow(cleared_systematic) > 0L) {
      cat("Systematic-residual filter (possible model misfit  --  flags cleared):\n")
      print(cleared_systematic[, c("compound", "n_flags_cleared", "residual_sign", "reason")],
            row.names = FALSE)
      cat("\n")
    }

    n_out   <- nrow(outlier_table)
    n_total <- nrow(results)
    cat(sprintf("Total dose points assessed : %d\n", n_total))
    cat(sprintf("Control rows excluded      : %d (never tested)\n", length(ctrl_rows)))
    pct_str <- if (n_total > 0L) sprintf("%.1f%%", 100 * n_out / n_total) else "N/A"
    cat(sprintf("Outliers (FDR Q=%.3f)      : %d (%s)\n", Q, n_out, pct_str))

    if (nrow(cleared_systematic) > 0L)
      cat(sprintf("Systematic flags cleared   : %d (see $cleared_systematic)\n",
                  sum(cleared_systematic$n_flags_cleared)))

    if (n_out > 0L) {
      cat("\nFlagged outliers:\n")
      print(outlier_table[, c("compound", "column", conc_col_name, "bret_ratio",
                              "fitted", "std_residual", "dynamic_range_pct")],
            row.names = FALSE, digits = 3)
    } else {
      cat("No outliers detected.\n")
    }
  }

  return(invisible(list(results            = results,
                        cleaned_table      = cleaned_table,
                        outlier_table      = outlier_table,
                        skipped_table      = skipped_table,
                        cleared_systematic = cleared_systematic,
                        params             = list(Q           = Q,
                                                  n_param     = n_param,
                                                  direction   = direction,
                                                  log_base    = log_base,
                                                  ntry_retry  = ntry_retry))))
}
