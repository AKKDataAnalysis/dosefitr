# ---------------------------------------------------------------------------
# test-plausibility-limits.R -- user-configurable biological-plausibility limits
# ---------------------------------------------------------------------------
#
# Covers the ten limit arguments shared by fit_drc_3pl(), fit_drc_4pl() and
# batch_drc_analysis():
#   bottom/top_limits_{nanobret,viability}_{inhibition,activation},
#   logIC50_limits, hill_slope_limits.
#
# All fixtures are synthetic (no inst/extdata dependency). Response vectors
# are length-12 (one NA header row + 11 plate rows), matching the wide-format
# convention used by test-display-consistency.R.

# --- Fixture builders --------------------------------------------------------

# Wide-format data frame: first column log-concentrations, then one compound
# in duplicate. `responses` must be length 12 (raw plate rows below the NA
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

# Clean inhibition curve, Bottom ~10, Top ~100, LogIC50 ~ -7. Classified
# deterministically as "inhibition"; passes all default plausibility limits.
pl_clean_inhibition <- function() {
  b <- 10; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# Clean activation curve, Bottom ~10, Top ~100, LogEC50 ~ -7. Classified
# deterministically as "activation"; passes all default plausibility limits.
pl_clean_activation <- function() {
  b <- 10; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^-(x - li))
}

# Flat curve with one wild outlier: the unconstrained 3PL Bottom fits far
# outside the default nanobret inhibition/flat limits and is corrected to the
# experimental minimum (50). Same fixture as test-display-consistency.R.
pl_wild_outlier <- function() {
  r <- rep(50, 12); r[1] <- 200; r
}

# Inhibition curve whose Bottom plateau sits ABOVE the default nanobret
# inhibition Bottom limit (600) while Top stays inside the default Top limit
# (0-700): Bottom ~ 650, Top ~ 700. Under default limits Bottom is corrected
# (clamped to 600); relaxing bottom_limits_nanobret_inhibition removes the
# correction in BOTH fitters. (The wild-outlier fixture only triggers a Bottom
# correction in 3PL, because the free 4PL Hill slope absorbs the outlier.)
pl_high_bottom_inhibition <- function() {
  b <- 650; t <- 700; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# Minimal one-plate batch_results list for batch_drc_analysis(): a plate
# element is recognised via $result$modified_ratio_table.
pl_batch_results <- function(data_table, assay_source = "nanobret") {
  br <- list(plate_01 = list(result = list(modified_ratio_table = data_table)))
  attr(br, "assay_source") <- assay_source
  br
}

# --- Defaults reproduce previous hardcoded behaviour -------------------------

test_that("default limits reproduce the hardcoded Bottom correction (3PL)", {
  fit <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_wild_outlier()), normalize = FALSE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  expect_true(isTRUE(r$biological_plausibility_check$needs_correction))
  # Corrected Bottom is clamped to the experimental minimum (50).
  expect_equal(unname(r$parameters$Value[1]), 50, tolerance = 1e-6)
})

test_that("default limits reproduce the hardcoded Bottom correction (4PL)", {
  # The wild-outlier fixture does not trigger a Bottom correction in 4PL (the
  # free Hill slope absorbs the outlier), so use the high-Bottom inhibition
  # fixture: Bottom ~ 800 exceeds the default inhibition limit of 600.
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_high_bottom_inhibition()), normalize = FALSE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_true("Bottom" %in% names(bpc$corrections_applied))
  # Corrected Bottom is clamped into the default limits.
  expect_lte(unname(r$parameters$Value[1]), 600)
})

# --- Relaxing a limit removes the correction ---------------------------------

test_that("relaxed bottom_limits_nanobret_inhibition keeps the fitted Bottom (3PL)", {
  # The wild-outlier raw 3PL fit reaches Bottom ~ 1.05e6, so the relaxed
  # limit must exceed that to fully remove the correction.
  fit <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_wild_outlier()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(-100, 1e7))
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  # Bottom no longer corrected (any remaining correction must not involve Bottom).
  expect_false("Bottom" %in% names(bpc$corrections_applied))
})

test_that("relaxed bottom_limits_nanobret_inhibition keeps the fitted Bottom (4PL)", {
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_high_bottom_inhibition()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(-100, 1e6))
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_false("Bottom" %in% names(bpc$corrections_applied))
  # Fitted Bottom retained near the true plateau (~650).
  expect_equal(unname(r$parameters$Value[1]), 650, tolerance = 0.05)
})

# --- Tightening a limit triggers a correction --------------------------------

test_that("tightened top_limits_nanobret_inhibition triggers Top correction (3PL)", {
  fit <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                top_limits_nanobret_inhibition = c(0, 90))
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_true("Top" %in% names(bpc$corrections_applied))
  # Corrected Top is clamped to the new upper limit.
  expect_lte(unname(r$parameters$Value[2]), 90)
})

test_that("tightened top_limits_nanobret_inhibition triggers Top correction (4PL)", {
  fit <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                top_limits_nanobret_inhibition = c(0, 90))
  )
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_true("Top" %in% names(bpc$corrections_applied))
  expect_lte(unname(r$parameters$Value[2]), 90)
})

# --- Direction split ----------------------------------------------------------

test_that("activation curves use the *_activation limits, not *_inhibition (4PL)", {
  # Tighten ONLY the inhibition Bottom limits: an activation curve must be
  # unaffected (its Bottom ~10 is inside the default activation limits).
  fit_act <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_activation()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(20, 600))
  )
  bpc_act <- fit_act$detailed_results[[1]]$biological_plausibility_check
  expect_false("Bottom" %in% names(bpc_act$corrections_applied))

  # The same tightened inhibition limit DOES bite an inhibition curve
  # (Bottom ~10 < 20).
  fit_inh <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_inhibition()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(20, 600))
  )
  bpc_inh <- fit_inh$detailed_results[[1]]$biological_plausibility_check
  expect_true("Bottom" %in% names(bpc_inh$corrections_applied))

  # Tighten ONLY the activation Bottom limits: now the activation curve is
  # corrected while the inhibition curve (default limits) is not.
  fit_act2 <- suppressWarnings(
    fit_drc_4pl(pl_build_data(pl_clean_activation()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_activation = c(20, Inf))
  )
  bpc_act2 <- fit_act2$detailed_results[[1]]$biological_plausibility_check
  expect_true("Bottom" %in% names(bpc_act2$corrections_applied))
})

test_that("activation curves use the *_activation limits, not *_inhibition (3PL)", {
  fit_act <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_clean_activation()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_inhibition = c(20, 600))
  )
  bpc_act <- fit_act$detailed_results[[1]]$biological_plausibility_check
  expect_false("Bottom" %in% names(bpc_act$corrections_applied))

  fit_act2 <- suppressWarnings(
    fit_drc_3pl(pl_build_data(pl_clean_activation()), normalize = FALSE, verbose = FALSE,
                bottom_limits_nanobret_activation = c(20, Inf))
  )
  bpc_act2 <- fit_act2$detailed_results[[1]]$biological_plausibility_check
  expect_true("Bottom" %in% names(bpc_act2$corrections_applied))
})

# --- Shared logIC50 / hill_slope limits ---------------------------------------

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

test_that("batch_drc_analysis forwards tightened limits to the fitter (3PL)", {
  br <- pl_batch_results(pl_build_data(pl_clean_inhibition()))
  res_default <- suppressWarnings(suppressMessages(
    batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE, model = "3pl")
  ))
  res_tight <- suppressWarnings(suppressMessages(
    batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE, model = "3pl",
                       top_limits_nanobret_inhibition = c(0, 90))
  ))
  get_fit <- function(res) res$drc_results[[1]]$drc_result$detailed_results[[1]]
  expect_false(isTRUE(get_fit(res_default)$biological_plausibility_check$needs_correction))
  expect_true(isTRUE(get_fit(res_tight)$biological_plausibility_check$needs_correction))
  expect_lte(unname(get_fit(res_tight)$parameters$Value[2]), 90)
})

test_that("batch_drc_analysis forwards tightened limits to the fitter (4PL)", {
  br <- pl_batch_results(pl_build_data(pl_clean_inhibition()))
  res_tight <- suppressWarnings(suppressMessages(
    batch_drc_analysis(br, generate_reports = FALSE, verbose = FALSE, model = "4pl",
                       top_limits_nanobret_inhibition = c(0, 90))
  ))
  fit <- res_tight$drc_results[[1]]$drc_result$detailed_results[[1]]
  expect_true(isTRUE(fit$biological_plausibility_check$needs_correction))
  expect_lte(unname(fit$parameters$Value[2]), 90)
})
