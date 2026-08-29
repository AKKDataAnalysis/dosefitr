# ---------------------------------------------------------------------------
# test-plausibility-limits.R -- user-configurable biological-plausibility limits
# ---------------------------------------------------------------------------
#
# Covers the ten limit arguments shared by fit_drc_3pl(), fit_drc_4pl() and
# batch_drc_analysis():
#   bottom/top_limits_{nanobret,viability}_{inhibition,activation},
#   logIC50_limits, hill_slope_limits.
#
# POLICY (raw-plateau policy): Bottom/Top plateaus are NEVER censored in the
# reported fit -- the reported value is always the raw fitted value, so the
# drawn curve tracks the data (GraphPad-consistent).  Degenerate fits are
# FLAGGED, not censored: a fit whose asymptote diverged (beyond the observed
# response range +/- 3x) or -- on the normalized 0-100 scale only -- whose
# plateau is wildly non-physical (outside [-50, 200]) is marked
# "Fit diverged - implausible plateau" in curve_quality (and the batch
# Exclusion column), with the raw value still reported.  The plateau limit
# arguments remain in force as (a) validation-checked inputs and (b) the port
# bounds of the constrained-fit recovery.  logIC50_limits and
# hill_slope_limits still censor (LogIC50 -> NA, HillSlope -> clamped
# default).
#
# All fixtures are synthetic (no inst/extdata dependency).  Response vectors
# are length-12 (one NA header row + 11 plate rows), matching the wide-format
# convention used by test-display-consistency.R.

# --- Fixture builders --------------------------------------------------------

# Wide-format data frame: first column log-concentrations, then one compound
# in duplicate.  `responses` must be length 12 (raw plate rows below the NA
# header row).
pl_build_data <- function(responses, compound_name = "Test.Cpd1") {
  stopifnot(length(responses) == 12L)
  df <- data.frame(
    log_conc = c(NA_real_, seq(-4.5, -10, length.out = 12)),
    r1       = c(NA_real_, responses),
    r2       = c(NA_real_, responses)
  )
  names(df) <- c("log_conc",
                 paste0(compound_name, "_rep1"),
                 paste0(compound_name, "_rep2"))
  df
}

# Wide-format data frame from per-replicate response vectors, with optional
# control anchor rows (first row response 0, last row response 100, both with
# NA log-conc) mirroring the batch_viability_analysis() output layout.  The
# anchors make normalize = TRUE an identity transform on a 0-100 fixture --
# the same mechanism that keeps the bundled plates unchanged under
# normalization.
pl_build_data2 <- function(r1, r2, compound_name = "Test.Cpd1", ctrl = FALSE) {
  stopifnot(length(r1) == 12L, length(r2) == 12L)
  lc <- c(NA_real_, seq(-4.5, -10, length.out = 12))
  a <- c(NA_real_, r1)
  b <- c(NA_real_, r2)
  if (ctrl) {
    lc <- c(lc, NA_real_)
    a  <- c(0, r1, 100)
    b  <- c(0, r2, 100)
  }
  df <- data.frame(log_conc = lc, a = a, b = b)
  names(df) <- c("log_conc",
                 paste0(compound_name, "_rep1"),
                 paste0(compound_name, "_rep2"))
  df
}

# Clean inhibition curve, Bottom ~10, Top ~100, LogIC50 ~ -7.  Classified
# deterministically as "inhibition"; passes all default plausibility limits.
pl_clean_inhibition <- function() {
  b <- 10; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# Clean activation curve, Bottom ~10, Top ~100, LogEC50 ~ -7.  Classified
# deterministically as "activation"; passes all default plausibility limits.
pl_clean_activation <- function() {
  b <- 10; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^-(x - li))
}

# Flat curve with one wild outlier: the unconstrained 3PL fit DIVERGES
# (Bottom ~ 1.05e6, far beyond the observed range +/- 3x), so the fit is
# flagged "Fit diverged - implausible plateau" and the raw value is reported.
# detect_curve_type() classifies the data as "flat".  Same fixture as
# test-display-consistency.R.
pl_wild_outlier <- function() {
  r <- rep(50, 12); r[1] <- 200; r
}

# Inhibition curve whose Bottom plateau sits at ~650 with Top ~700: outside
# the normalized-scale extreme band, but a genuine raw-scale NanoBRET ratio
# plateau.  Reported raw, uncensored, unflagged ("Good curve") in BOTH
# fitters.  (The wild-outlier fixture diverges only in 3PL, because the free
# 4PL Hill slope absorbs the outlier.)
pl_high_bottom_inhibition <- function() {
  b <- 650; t <- 700; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# Partial-plateau inhibition curve on the normalized scale: Bottom ~ 84,
# Top ~ 100, LogIC50 ~ -7.  Mirrors A375:LW456 (partial viability plateau):
# must be reported as fitted -- never clamped toward 60 -- and not flagged.
pl_partial_plateau <- function() {
  b <- 84; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# Inhibition curve with a NEGATIVE Bottom plateau (-80) on the normalized
# scale: outside the extreme band [-50, 200], so it is flagged as degenerate
# when normalize = TRUE.  On the raw (non-normalized) scale the band is
# inert and the same fit is a clean "Good curve".
pl_neg_plateau <- function() {
  b <- -80; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# LW579-like edge response (per-replicate): deep drop at the highest dose,
# partial at the second, plateau ~100 elsewhere.  The unconstrained 4PL fit
# diverges; the constrained-fit recovery then bounds the plateaus and is
# accepted ("IC50 from constrained fit (plateaus bounded)").
pl_edge_recovery <- function() {
  list(r1 = c(5, 72, 90, 98, 100, 99, 100, 101, 100, 99, 100, 101),
       r2 = c(6, 74, 88, 99, 100, 101, 99, 100, 101, 100, 99, 100))
}

# Minimal one-plate batch_results list for batch_drc_analysis(): a plate
# element is recognised via $result$modified_ratio_table.
pl_batch_results <- function(data_table, assay_source = "nanobret") {
  br <- list(plate_01 = list(result = list(modified_ratio_table = data_table)))
  attr(br, "assay_source") <- assay_source
  br
}

# --- Flag, don't censor: diverged degenerate fits -----------------------------

test_that("diverged degenerate fit is flagged, not censored (3PL)", {
  fit <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_wild_outlier()), normalize = FALSE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  # Flag-only: no parameter is censored ...
  expect_length(bpc$corrections_applied, 0L)
  expect_true(isTRUE(bpc$fit_diverged))
  # ... the raw diverged Bottom (~1.05e6) is still reported ...
  expect_gt(abs(unname(r$parameters$Value[1])), 1e5)
  # ... and the fit carries the divergence flag.
  expect_match(r$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
})

# --- Raw plateaus are reported as fitted --------------------------------------

test_that("raw NanoBRET-scale plateau is neither censored nor flagged (3PL and 4PL)", {
  # Bottom ~ 650 / Top ~ 700: outside the normalized-scale extreme band but a
  # genuine raw-scale NanoBRET ratio plateau.  Reported raw, no flag.
  for (fun in list(fit_drc_3pl, fit_drc_4pl)) {
    fit <- suppressWarnings(
      fun(pl_build_data(pl_high_bottom_inhibition()), normalize = FALSE, verbose = FALSE)
    )
    r <- fit$detailed_results[[1]]
    expect_false(isTRUE(r$biological_plausibility_check$needs_correction))
    expect_equal(unname(r$parameters$Value[1]), 650, tolerance = 0.05)
    expect_equal(unname(r$parameters$Value[2]), 700, tolerance = 0.05)
    expect_match(r$curve_quality, "Good curve", fixed = TRUE)
  }
})

test_that("partial plateau on the normalized scale is reported raw and unflagged (3PL and 4PL)", {
  # Regression test for the A375:LW456 GraphPad mismatch: a partial viability
  # plateau at ~84% must be reported as fitted (not clamped toward 60) so the
  # drawn curve tracks the points.
  for (fun in list(fit_drc_3pl, fit_drc_4pl)) {
    fit <- suppressWarnings(
      fun(pl_build_data2(pl_partial_plateau(), pl_partial_plateau(), ctrl = TRUE),
          normalize = TRUE, verbose = FALSE)
    )
    r <- fit$detailed_results[[1]]
    expect_false(isTRUE(r$biological_plausibility_check$needs_correction))
    expect_equal(unname(r$parameters$Value[1]), 84, tolerance = 0.05)
    expect_false(grepl("Fit diverged", r$curve_quality, fixed = TRUE))
  }
})

test_that("extreme plateau on the normalized scale is flagged but reported raw (3PL and 4PL)", {
  # Bottom = -80 on the normalized 0-100 scale is outside the extreme band
  # [-50, 200]: flagged as degenerate, but the raw value is still reported.
  for (fun in list(fit_drc_3pl, fit_drc_4pl)) {
    fit <- suppressWarnings(
      fun(pl_build_data2(pl_neg_plateau(), pl_neg_plateau(), ctrl = TRUE),
          normalize = TRUE, verbose = FALSE)
    )
    r <- fit$detailed_results[[1]]
    bpc <- r$biological_plausibility_check
    expect_true(isTRUE(bpc$needs_correction))
    expect_length(bpc$corrections_applied, 0L)
    # Band flag, not a divergence flag.
    expect_false(isTRUE(bpc$fit_diverged))
    expect_equal(unname(r$parameters$Value[1]), -80, tolerance = 0.05)
    expect_match(r$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
  }
})

test_that("extreme band is inert on the raw (non-normalized) scale (3PL and 4PL)", {
  # The fixed [-50, 200] band assumes the normalized 0-100 scale; on raw
  # NanoBRET ratios a plateau at -80 is unremarkable and must NOT be flagged.
  for (fun in list(fit_drc_3pl, fit_drc_4pl)) {
    fit <- suppressWarnings(
      fun(pl_build_data2(pl_neg_plateau(), pl_neg_plateau(), ctrl = TRUE),
          normalize = FALSE, verbose = FALSE)
    )
    r <- fit$detailed_results[[1]]
    expect_false(isTRUE(r$biological_plausibility_check$needs_correction))
    expect_equal(unname(r$parameters$Value[1]), -80, tolerance = 0.05)
    expect_match(r$curve_quality, "Good curve", fixed = TRUE)
  }
})

test_that("non-plateau corrections still apply (4PL)", {
  # The free 4PL Hill slope absorbs the wild outlier and lands outside the
  # hill-slope limits: the HillSlope correction still fires (main path), the
  # fit is labelled "Implausible fit", and the plateau flag does NOT apply.
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_wild_outlier()), normalize = FALSE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_true("HillSlope" %in% names(bpc$corrections_applied))
  expect_equal(bpc$fit_label, "Implausible fit")
  expect_false(isTRUE(bpc$fit_diverged))
  expect_false(grepl("Fit diverged - implausible plateau", r$curve_quality, fixed = TRUE))
})

# --- Plateau limit arguments never censor the reported fit --------------------

test_that("plateau limit arguments never censor the reported fit (3PL)", {
  # The divergence flag is scale-free: relaxing the bottom limits neither
  # censors nor unflags the diverged wild-outlier fit.
  fit <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_wild_outlier()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(-100, 1e7))
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_length(bpc$corrections_applied, 0L)
  expect_gt(abs(unname(r$parameters$Value[1])), 1e5)
  expect_match(r$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)

  # Tightening the top limits does not censor a clean fit either.
  fit_tight <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                top_limits_nanobret_inhibition = c(0, 90))
  )
  r_tight <- fit_tight$detailed_results[[1]]
  expect_false(isTRUE(r_tight$biological_plausibility_check$needs_correction))
  expect_equal(unname(r_tight$parameters$Value[2]), 100, tolerance = 0.05)
})

test_that("plateau limit arguments never censor the reported fit (4PL)", {
  # Relaxed bottom limits: raw-scale high plateau retained, unflagged.
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_high_bottom_inhibition()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(-100, 1e6))
  )
  r <- fit$detailed_results[[1]]
  expect_false(isTRUE(r$biological_plausibility_check$needs_correction))
  expect_equal(unname(r$parameters$Value[1]), 650, tolerance = 0.05)

  # Tightened top limits: clean fit untouched.
  fit_tight <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                top_limits_nanobret_inhibition = c(0, 90))
  )
  r_tight <- fit_tight$detailed_results[[1]]
  expect_false(isTRUE(r_tight$biological_plausibility_check$needs_correction))
  expect_equal(unname(r_tight$parameters$Value[2]), 100, tolerance = 0.05)
})

test_that("direction-specific plateau limits are validated but do not censor (3PL and 4PL)", {
  # Tightening either direction's bottom limits leaves both curve directions
  # untouched: the arguments remain validated (fail-fast tests below) and
  # still bound the constrained-fit recovery, but never censor a reported fit.
  for (fun in list(fit_drc_3pl, fit_drc_4pl)) {
    cases <- list(
      list(data = pl_clean_inhibition(), arg = "bottom_limits_nanobret_inhibition"),
      list(data = pl_clean_activation(), arg = "bottom_limits_nanobret_inhibition"),
      list(data = pl_clean_activation(), arg = "bottom_limits_nanobret_activation")
    )
    for (case in cases) {
      args <- list(data = pl_build_data(case$data), normalize = FALSE, verbose = FALSE)
      args[[case$arg]] <- c(20, 600)
      fit <- suppressWarnings(do.call(fun, args))
      r <- fit$detailed_results[[1]]
      expect_false(isTRUE(r$biological_plausibility_check$needs_correction))
      expect_equal(unname(r$parameters$Value[1]), 10, tolerance = 0.05)
    }
  }
})

# --- Constrained-fit recovery (unchanged by the raw-plateau policy) -----------

test_that("constrained-fit recovery still bounds plateaus (4PL)", {
  # LW579-like edge response: the unconstrained fit diverges, the data show a
  # genuine response at the highest doses, and the constrained refit is
  # accepted -- unchanged by the raw-plateau policy.
  rec <- pl_edge_recovery()
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data2(rec$r1, rec$r2), normalize = FALSE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  expect_match(r$curve_quality, "IC50 from constrained fit (plateaus bounded)", fixed = TRUE)
  # Bounded Bottom sits at the lower nanobret inhibition limit (-100).
  expect_gte(unname(r$parameters$Value[1]), -100 - 1e-6)
  expect_lte(unname(r$parameters$Value[1]), 600)
  # Recovered LogIC50 is finite and essentially in-window.
  li <- unname(r$parameters$Value[3])
  expect_true(is.finite(li))
  expect_gte(li, -11)
  expect_lte(li, -3.5)
})

test_that("recovery acceptance is not blocked by the plateau flag when normalized (4PL)", {
  # The accepted bounded fit lands at Bottom = -100, OUTSIDE the normalized
  # extreme band [-50, 200].  The band must not reject the constrained
  # recovery (the bounded fit is judged by the pre-flag criterion).
  rec <- pl_edge_recovery()
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data2(rec$r1, rec$r2, ctrl = TRUE), normalize = TRUE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  expect_match(r$curve_quality, "IC50 from constrained fit (plateaus bounded)", fixed = TRUE)
  expect_equal(unname(r$parameters$Value[1]), -100, tolerance = 1e-6)
})

test_that("plateau limits shape the constrained-fit recovery bounds (4PL)", {
  # The limit arguments remain the port bounds of the constrained recovery:
  # tightening the bottom limit moves the bounded Bottom accordingly.
  rec <- pl_edge_recovery()
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data2(rec$r1, rec$r2), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(-50, 600))
  )
  r <- fit$detailed_results[[1]]
  expect_match(r$curve_quality, "IC50 from constrained fit (plateaus bounded)", fixed = TRUE)
  expect_equal(unname(r$parameters$Value[1]), -50, tolerance = 1e-6)
})

# --- Shared logIC50 / hill_slope limits (still censor) ------------------------

test_that("tightened logIC50_limits sets LogIC50 to NA (3PL and 4PL)", {
  # Clean curve has LogIC50 ~ -7; exclude it with a window that omits -7.
  narrow <- c(-6, 5)
  f3 <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                logIC50_limits = narrow)
  )
  expect_true(is.na(f3$detailed_results[[1]]$parameters$Value[3]))
  f4 <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                logIC50_limits = narrow)
  )
  expect_true(is.na(f4$detailed_results[[1]]$parameters$Value[3]))
})

test_that("tightened hill_slope_limits resets HillSlope to the default (4PL)", {
  # Clean inhibition curve fits HillSlope ~ -1; forbid it.
  f4 <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                hill_slope_limits = c(-0.5, 0.5))
  )
  r <- f4$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true("HillSlope" %in% names(bpc$corrections_applied))
  # Reset to the inhibition default (-1) clamped into the new limits -> -0.5.
  expect_equal(unname(r$parameters$Value[4]), -0.5, tolerance = 1e-6)
})

test_that("hill_slope_limits is inert but validated in fit_drc_3pl", {
  # Valid but impossible range: no error, no effect (Hill fixed at -1).
  f3 <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                hill_slope_limits = c(-0.5, 0.5))
  )
  expect_true(isTRUE(f3$detailed_results[[1]]$success))
  # Malformed input still fails fast.
  expect_error(
    fit_drc_3pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                hill_slope_limits = c(1, -1)),
    "hill_slope_limits"
  )
})

# --- Validation (fail fast, informative) --------------------------------------

test_that("invalid limit arguments fail fast with the argument name (3PL)", {
  d <- pl_build_data(pl_clean_inhibition())
  expect_error(fit_drc_3pl(d, verbose = FALSE, logIC50_limits = c(5, -5)),
               "logIC50_limits")
  expect_error(fit_drc_3pl(d, verbose = FALSE, logIC50_limits = c(1, 2, 3)),
               "logIC50_limits")
  expect_error(fit_drc_3pl(d, verbose = FALSE, logIC50_limits = "a"),
               "logIC50_limits")
  expect_error(fit_drc_3pl(d, verbose = FALSE, logIC50_limits = c(NA, 5)),
               "logIC50_limits")
  expect_error(fit_drc_3pl(d, verbose = FALSE,
                           bottom_limits_viability_activation = c(60, -20)),
               "bottom_limits_viability_activation")
})

test_that("invalid limit arguments fail fast with the argument name (4PL)", {
  d <- pl_build_data(pl_clean_inhibition())
  expect_error(fit_drc_4pl(d, verbose = FALSE, hill_slope_limits = c(5, -5)),
               "hill_slope_limits")
  expect_error(fit_drc_4pl(d, verbose = FALSE, top_limits_nanobret_inhibition = 700),
               "top_limits_nanobret_inhibition")
  expect_error(fit_drc_4pl(d, verbose = FALSE,
                           bottom_limits_nanobret_activation = c(NA, Inf)),
               "bottom_limits_nanobret_activation")
})

test_that("invalid limit arguments fail fast in batch_drc_analysis", {
  br <- pl_batch_results(pl_build_data(pl_clean_inhibition()))
  expect_error(
    batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE,
                       logIC50_limits = c(2, -2)),
    "logIC50_limits"
  )
  expect_error(
    batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE,
                       top_limits_viability_inhibition = c(1, 2, 3)),
    "top_limits_viability_inhibition"
  )
})

# --- batch_drc_analysis pass-through ------------------------------------------

test_that("batch_drc_analysis forwards logIC50_limits to the fitter (3PL and 4PL)", {
  # logIC50_limits still censors, so it exercises the pass-through directly:
  # the clean curve's LogIC50 ~ -7 falls outside the narrowed window.
  br <- pl_batch_results(pl_build_data(pl_clean_inhibition()))
  get_fit <- function(res) res$drc_results[[1]]$drc_result$detailed_results[[1]]
  for (m in c("3pl", "4pl")) {
    res_default <- suppressWarnings(suppressMessages(
      batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE, model = m)
    ))
    res_tight <- suppressWarnings(suppressMessages(
      batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE, model = m,
                         logIC50_limits = c(-6, 5))
    ))
    expect_false(isTRUE(get_fit(res_default)$biological_plausibility_check$needs_correction))
    expect_true(isTRUE(get_fit(res_tight)$biological_plausibility_check$needs_correction))
    expect_true(is.na(get_fit(res_tight)$parameters$Value[3]))
  }
})

test_that("batch_drc_analysis no longer censors plateaus via top_limits (3PL and 4PL)", {
  br <- pl_batch_results(pl_build_data(pl_clean_inhibition()))
  get_fit <- function(res) res$drc_results[[1]]$drc_result$detailed_results[[1]]
  for (m in c("3pl", "4pl")) {
    res_tight <- suppressWarnings(suppressMessages(
      batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE, model = m,
                         top_limits_nanobret_inhibition = c(0, 90))
    ))
    fit <- get_fit(res_tight)
    expect_false(isTRUE(fit$biological_plausibility_check$needs_correction))
    expect_equal(unname(fit$parameters$Value[2]), 100, tolerance = 0.05)
  }
})
