# -----------------------------------------------------------------------------
# test-display-consistency.R -- IC50/plot consistency fixes + raw-plateau policy
# -----------------------------------------------------------------------------
#
# Regression tests for the three optional fixes documented in
# report_ic50_plot_consistency.md, plus the raw-plateau policy tests at the
# bottom of the file (reported Bottom/Top are ALWAYS the raw fitted values;
# degenerate fits are flagged "Fit diverged - implausible plateau", never
# censored):
#
#   Fix A  Vertical IC50 line source in plot_dose_response() reads the
#          *corrected* LogIC50 from result$parameters (not the raw drm coef),
#          matching plot_multiple_compounds().
#
#   Fix B  The fitted curve line is drawn from analytic_dose_response() using
#          the corrected parameters whenever biological_plausibility_check
#          fired.  Prevents the "curve overshoots plot" divergence when the
#          raw drm Bottom / Top are wildly out of range but the batch report
#          would show a corrected value.
#
#   Fix C  Optional propagation of the "N/D" and ">highest_conc" display
#          overrides used by batch_drc_analysis()'s Pharmacology_Summary.
#          Opt-in via show_display_overrides = TRUE + show_display_badge = TRUE.
#
# All three fixes are backward compatible: default behaviour on well-fit
# compounds is unchanged.  New args default to FALSE.
#
# Fixtures are tiny synthetic datasets (12 log-concentrations) rather than
# the bundled NanoBRET plates so behaviour is deterministic and fast across
# CRAN/Bioc check environments.  The fit result is accessed via
# fit$detailed_results[[i]] (the shape plot_dose_response() actually
# indexes; see R/plot_dose_response.R:365).

# --- Fixture builders --------------------------------------------------------

# Log-concentration column that mirrors the NanoBRET plate layout used by
# the bundled fixtures: NA header row + 12 concentrations, highest at
# 10^-4.5 (~31.6 uM), lowest at 10^-10 (~0.1 pM).
lc_seq <- function() {
  c(NA_real_, seq(-4.5, -10, length.out = 12))
}

# Build a wide-format data frame ready for fit_drc_3pl().  Response values
# must be a length-12 numeric vector (the raw plate rows below the NA
# header row).
build_data <- function(responses, compound_name = "Test.Cpd1") {
  stopifnot(length(responses) == 12L)
  df <- data.frame(
    log_conc = lc_seq(),
    r1       = c(NA_real_, responses),
    r2       = c(NA_real_, responses)
  )
  names(df) <- c("log_conc",
                 paste0(compound_name, "_rep1"),
                 paste0(compound_name, "_rep2"))
  df
}

# A "clean" inhibition curve, LogIC50 ~ -7, no correction needed.  Curve type
# is deterministically classified as "inhibition".
clean_responses <- function() {
  b <- 10;  t <- 100;  li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# A "wild" outlier response set whose unconstrained 3PL fit DIVERGES (raw
# Bottom ~ 1.05e6, far beyond the observed range +/- 3x).  Under the
# raw-plateau policy the fit is flagged "Fit diverged - implausible plateau"
# and the raw value is reported (flag, don't censor).
# detect_curve_type() classifies this as "flat" because the noise-aware
# threshold dominates the head-vs-tail-3 mean difference.
wild_outlier_responses <- function() {
  r <- rep(50, 12); r[1] <- 200; r
}

# Fit a single-compound batch and return the fit object.
fit_one <- function(responses, compound_name = "Test.Cpd1") {
  suppressWarnings(
    fit_drc_3pl(data      = build_data(responses, compound_name),
                normalize = FALSE, verbose = FALSE)
  )
}

# --- ggplot layer extractors -------------------------------------------------

# Every vertical line (geom_vline layer) in a ggplot object.  Returns the
# xintercept values as an *unnamed* numeric vector (some ggplot layers
# attach a `name` attribute we don't want to compare against).
vline_xintercepts <- function(gg) {
  out <- vapply(gg$layers, function(L) {
    geom_cls <- class(L$geom)[1]
    if (identical(geom_cls, "GeomVline")) {
      xi <- L$data$xintercept
      if (length(xi) > 0) xi[1] else NA_real_
    } else NA_real_
  }, numeric(1))
  unname(out[!is.na(out)])
}

# All non-empty labels from geom_text layers.  ggplot2::annotate() stores
# the label inside L$aes_params$label rather than L$data$label, so we check
# both locations.
annotation_labels <- function(gg) {
  out <- vapply(gg$layers, function(L) {
    geom_cls <- class(L$geom)[1]
    if (!identical(geom_cls, "GeomText")) return(NA_character_)
    if (!is.null(L$data) && "label" %in% names(L$data)) {
      return(as.character(L$data$label[1]))
    }
    if (!is.null(L$aes_params) && "label" %in% names(L$aes_params)) {
      return(as.character(L$aes_params$label))
    }
    NA_character_
  }, character(1))
  out[!is.na(out)]
}

# The y values from the first geom_line layer.  Used to verify Fix B (the
# curve line should sit within the corrected Bottom/Top band, not spike out).
line_y <- function(gg) {
  for (L in gg$layers) {
    geom_cls <- class(L$geom)[1]
    if (identical(geom_cls, "GeomLine")) {
      d <- L$data
      if (!is.null(d) && "y" %in% names(d)) return(d$y)
      if (!is.null(d) && "response" %in% names(d)) return(d$response)
    }
  }
  numeric(0)
}

# --- Fix A -------------------------------------------------------------------

test_that("Fix A: vline matches parameters$Value LogIC50 for a clean fit", {
  # No correction fires here, so raw coef == corrected parameters$Value.
  # This is the backward-compatibility baseline.
  fit <- fit_one(clean_responses())
  expect_gte(length(fit$detailed_results), 1L)
  r  <- fit$detailed_results[[1]]
  li_from_params <- r$parameters$Value[r$parameters$Parameter == "LogIC50"]
  expect_true(is.finite(li_from_params[1]))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  xi <- vline_xintercepts(gg)
  expect_length(xi, 1L)
  expect_equal(xi[1], li_from_params[1], tolerance = 1e-10)
})

test_that("Fix A: vline is suppressed by default when the fit is classified flat", {
  # Wild outlier -> Bottom correction fires and detect_curve_type() classifies
  # the fit as "flat".  A flat classification means the IC50 is not defined,
  # so the default plot must not draw an IC50 reference line at the
  # (meaningless) corrected LogIC50.  The corrected value itself is still
  # reported in parameters / the batch summary.
  fit <- fit_one(wild_outlier_responses())
  r <- fit$detailed_results[[1]]
  expect_true(isTRUE(r$biological_plausibility_check$needs_correction))
  expect_true(isTRUE(r$curve_type == "flat"))
  li_corrected <- r$parameters$Value[r$parameters$Parameter == "LogIC50"][1]
  expect_true(is.finite(li_corrected))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  expect_length(vline_xintercepts(gg), 0L)
})

# --- Fix B -------------------------------------------------------------------

test_that("Fix B: flat-classified curve is drawn flat at the median response", {
  # Wild outlier: the unconstrained fit diverges (raw Bottom ~ 1.05e6) and is
  # flagged "Fit diverged - implausible plateau" -- the raw value is reported
  # (flag, don't censor).  For flat-classified fits the plot draws a
  # horizontal line at the median of the per-dose mean responses -- the
  # median is robust to the very outliers that triggered the flat
  # classification.
  fit <- fit_one(wild_outlier_responses())
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_length(bpc$corrections_applied, 0L)
  expect_true(isTRUE(bpc$fit_diverged))
  expect_match(r$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
  expect_gt(abs(unname(r$parameters$Value[1])), 1e5)
  expect_true(isTRUE(r$curve_type == "flat"))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  ys <- line_y(gg)
  expect_gt(length(ys), 0L)
  expect_true(all(is.finite(ys)))
  expect_true(all(ys == ys[1]),
              info = sprintf("curve not flat: y range = [%.4f, %.4f]",
                             min(ys), max(ys)))
  expect_equal(ys[1], stats::median(wild_outlier_responses()), tolerance = 1e-6)
})

test_that("Fix B: clean fit still uses predict() and stays finite", {
  # Guardrail: an uncorrected fit must not silently switch curve sources.
  # We don't assert an exact match with predict() -- just finite + within
  # the sigmoid response range.
  fit <- fit_one(clean_responses())
  r <- fit$detailed_results[[1]]
  expect_false(isTRUE(r$biological_plausibility_check$needs_correction))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  ys <- line_y(gg)
  expect_gt(length(ys), 10L)
  expect_true(all(is.finite(ys)))
  # Sigmoid range: min slightly below Bottom (10), max slightly above Top
  # (100); a generous [-5, 120] envelope is enough here.
  expect_true(all(ys >= -5 & ys <= 120),
              info = sprintf("clean y range = [%.3f, %.3f]", min(ys), max(ys)))
})

test_that("Fix B: diverged out-of-window fit is flagged and draws a raw-model tendency line", {
  # The fixture rises sharply at the low-concentration edge (steeper than the
  # fixed unit Hill slope can fit in-window), so the unconstrained 3PL fit
  # diverges (raw Top ~ 3.5e5) with LogIC50 below the tested range.  Under
  # the raw-plateau policy the fit is flagged "Fit diverged - implausible
  # plateau" and reported raw -- NOT reclassified flat and NOT censored.
  # (Pre-policy this fixture was relabelled "IC50 below tested range (N/D)"
  # via the plateau clamp; the N/D labels now require a LogIC50/HillSlope
  # correction to fire, which clean synthetic shapes no longer trigger.)
  # The plot draws the RAW model as a tendency line (following the points)
  # instead of the flat data level; the IC50 vline stays suppressed because
  # the fitted LogIC50 is outside the dose window.
  nd_below <- c(20, 20, 21, 20, 21, 20, 21, 20, 21, 22, 40, 95)
  fit <- fit_one(nd_below)
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  expect_true(isTRUE(bpc$needs_correction))
  expect_length(bpc$corrections_applied, 0L)
  expect_true(isTRUE(bpc$fit_diverged))
  expect_match(r$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
  # Flag, don't censor: the raw diverged Top is reported, and the curve type
  # is left alone (no flat reclassification).
  expect_gt(abs(unname(r$parameters$Value[2])), 1e5)
  expect_false(isTRUE(r$curve_type == "flat"))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  ys <- line_y(gg)
  expect_gt(length(ys), 10L)
  expect_true(all(is.finite(ys)))
  # Not the flat data-level line: the tendency line must vary with dose and
  # track the raw model predictions at the tested concentrations.
  expect_gt(stats::sd(ys), 1)
  L <- NULL
  for (ly in gg$layers) {
    if (identical(class(ly$geom)[1], "GeomLine")) { L <- ly; break }
  }
  cx <- if ("x" %in% names(L$data)) L$data$x else L$data$log_inhibitor
  cy <- if ("y" %in% names(L$data)) L$data$y else L$data$response
  xs <- sort(unique(r$data$log_inhibitor))
  drawn_at <- stats::approx(cx, cy, xout = xs, ties = "ordered")$y
  raw_at <- as.numeric(stats::predict(r$model,
                                      newdata = data.frame(log_inhibitor = xs)))
  expect_lt(max(abs(drawn_at - raw_at)), 2)
  # IC50 vline stays suppressed (the fitted LogIC50 is out of window).
  expect_length(vline_xintercepts(gg), 0L)
})

test_that("plot_multiple_compounds Fix B: diverged out-of-window compound draws tendency line", {
  # Same rule in the batch plotter: a diverged out-of-window compound (raw
  # Top exploded, LogIC50 below the tested range, flagged "Fit diverged -
  # implausible plateau") gets the raw-model tendency line (excluded from the
  # y-range, clipped at the panel edge), while a genuine no-response flat
  # compound keeps the flat data-level line.
  nd_below <- c(20, 20, 21, 20, 21, 20, 21, 20, 21, 22, 40, 95)
  df <- data.frame(
    log_conc = lc_seq(),
    Test.nd_rep1   = c(NA_real_, nd_below),
    Test.nd_rep2   = c(NA_real_, nd_below),
    Test.wild_rep1 = c(NA_real_, wild_outlier_responses()),
    Test.wild_rep2 = c(NA_real_, wild_outlier_responses())
  )
  fit <- suppressWarnings(
    fit_drc_3pl(data = df, normalize = FALSE, verbose = FALSE)
  )
  cmpds <- vapply(fit$detailed_results, function(x) x$compound, character(1))
  nd_idx <- which(grepl("nd", cmpds))
  expect_gt(length(nd_idx), 0L)
  r_nd <- fit$detailed_results[[nd_idx[1]]]
  expect_match(r_nd$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
  expect_false(isTRUE(r_nd$curve_type == "flat"))

  out_dir <- tempfile("plmc_nd_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  gg <- suppressWarnings(plot_multiple_compounds(
    results = fit, compound_indices = seq_along(cmpds),
    save_plot = file.path(out_dir, "multi_nd.png"),
    verbose = FALSE
  ))

  line_dat <- NULL
  for (L in gg$layers) {
    if (identical(class(L$geom)[1], "GeomLine")) { line_dat <- L$data; break }
  }
  expect_true(!is.null(line_dat) && nrow(line_dat) > 0)
  col <- intersect(names(line_dat), c("compound", "Compound", "group", "colour"))[1]
  y_col <- intersect(names(line_dat), c("y", "response", "value"))[1]
  nd_rows <- line_dat[grepl("nd", line_dat[[col]]), ]
  expect_gt(nrow(nd_rows), 10L)
  ys <- nd_rows[[y_col]]
  expect_true(all(is.finite(ys)))
  # Tendency line: varies with dose and stays near the point envelope
  # (data range 20..95, clipped to +/-10%).
  expect_gt(stats::sd(ys), 1)
  expect_true(all(ys >= 20 - 0.10 * 75 - 1 & ys <= 95 + 0.10 * 75 + 1),
              info = sprintf("nd tendency y range = [%.3f, %.3f]", min(ys), max(ys)))
  # Contrast: the genuine-flat wild compound still gets the flat level.
  wild_rows <- line_dat[grepl("wild", line_dat[[col]]), ]
  expect_gt(nrow(wild_rows), 0L)
  expect_true(all(wild_rows[[y_col]] == wild_rows[[y_col]][1]),
              info = "wild compound curve should remain the flat data level")
})

# --- Fix C -------------------------------------------------------------------

test_that("Fix C: default (overrides OFF) suppresses the flat-fit vline, no badge", {
  # Approved behaviour change: the IC50 line is suppressed by default for
  # flat-classified fits (an IC50 is not defined for "no dose response").
  # No override badge is shown unless explicitly requested.
  fit <- fit_one(wild_outlier_responses())
  r <- fit$detailed_results[[1]]
  expect_true(isTRUE(r$curve_type == "flat"))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  expect_length(vline_xintercepts(gg), 0L)
  labs <- annotation_labels(gg)
  expect_false(any(grepl("N/D|IC50: >", labs)))
})

test_that("Fix C: vline is suppressed by default when LogIC50 is out of window", {
  # A non-flat fit whose fitted LogIC50 lands outside the tested dose window
  # must not draw the IC50 line by default -- the marker would sit off the
  # data range and stretch the x-axis.  Simulated by shifting the stored
  # LogIC50 of a clean in-window fit below the tested range.
  fit <- fit_one(clean_responses())
  r <- fit$detailed_results[[1]]
  expect_false(isTRUE(r$curve_type == "flat"))
  li_row <- which(r$parameters$Parameter == "LogIC50")[1]
  fit$detailed_results[[1]]$parameters$Value[li_row] <-
    min(fit$detailed_results[[1]]$data$log_inhibitor, na.rm = TRUE) - 5

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  expect_length(vline_xintercepts(gg), 0L)
})

test_that("Fix C: show_display_overrides=TRUE drops the vline for flat curve", {
  fit <- fit_one(wild_outlier_responses())
  r <- fit$detailed_results[[1]]
  expect_true(isTRUE(r$curve_type == "flat"))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            show_display_overrides = TRUE,
                                            verbose = FALSE))
  expect_length(vline_xintercepts(gg), 0L)
  labs <- annotation_labels(gg)
  # Badge is separate opt-in; must not appear yet.
  expect_false(any(grepl("^IC50:", labs)))
})

test_that("Fix C: show_display_badge=TRUE adds an 'IC50: N/D' text annotation", {
  fit <- fit_one(wild_outlier_responses())
  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            show_display_overrides = TRUE,
                                            show_display_badge = TRUE,
                                            verbose = FALSE))
  labs <- annotation_labels(gg)
  expect_true(any(grepl("^IC50: N/D$", labs)),
              info = paste("annotations:", paste(labs, collapse = " | ")))
})

# --- plot_multiple_compounds Fix B / Fix C -----------------------------------

test_that("plot_multiple_compounds Fix B: wild curve stays inside corrected band", {
  # Two-compound batch: one clean, one wild.  The wild compound's fit
  # diverges and is flagged (raw Bottom ~ 1.05e6 reported, not censored);
  # classified flat, its drawn curve must be the flat data level, not blow
  # up like raw drm predict() did before Fix B.
  df <- data.frame(
    log_conc = lc_seq(),
    Test.clean_rep1 = c(NA_real_, clean_responses()),
    Test.clean_rep2 = c(NA_real_, clean_responses()),
    Test.wild_rep1  = c(NA_real_, wild_outlier_responses()),
    Test.wild_rep2  = c(NA_real_, wild_outlier_responses())
  )
  fit <- suppressWarnings(
    fit_drc_3pl(data = df, normalize = FALSE, verbose = FALSE)
  )
  cmpds <- vapply(fit$detailed_results, function(x) x$compound, character(1))
  wild_idx <- which(grepl("wild", cmpds))
  expect_gt(length(wild_idx), 0L)
  r_wild <- fit$detailed_results[[wild_idx[1]]]
  expect_true(isTRUE(r_wild$biological_plausibility_check$needs_correction))

  out_dir <- tempfile("plmc_disp_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  gg <- suppressWarnings(plot_multiple_compounds(
    results = fit, compound_indices = seq_along(cmpds),
    save_plot = file.path(out_dir, "multi.png"),
    verbose = FALSE
  ))

  # Pull the curve-line data and check the wild compound's y values.
  line_dat <- NULL
  for (L in gg$layers) {
    if (identical(class(L$geom)[1], "GeomLine")) {
      line_dat <- L$data; break
    }
  }
  expect_true(!is.null(line_dat) && nrow(line_dat) > 0)
  cand <- intersect(names(line_dat), c("compound", "Compound", "group", "colour"))
  expect_gt(length(cand), 0L)
  col <- cand[1]
  wild_rows <- line_dat[grepl("wild", line_dat[[col]]), ]
  expect_gt(nrow(wild_rows), 0L)
  y_cand <- intersect(names(wild_rows), c("y", "response", "value"))
  expect_gt(length(y_cand), 0L)
  ys <- wild_rows[[y_cand[1]]]
  expect_true(all(is.finite(ys)))
  # Flat-classified diverged fit -> flat data-level line near 50.
  # Envelope of [0, 120] is generous vs the pre-fix runaway to y > 180.
  expect_true(all(ys >= 0 & ys <= 120),
              info = sprintf("wild curve y range = [%.3f, %.3f]",
                             min(ys), max(ys)))
})

test_that("plot_multiple_compounds Fix C: overrides drop vlines and stack badges", {
  # Two wild-outlier compounds so both get classified as flat -> both are
  # overridden, and the badge stacks vertically with '\n'.
  df <- data.frame(
    log_conc = lc_seq(),
    Test.wildA_rep1 = c(NA_real_, wild_outlier_responses()),
    Test.wildA_rep2 = c(NA_real_, wild_outlier_responses()),
    Test.wildB_rep1 = c(NA_real_, wild_outlier_responses()),
    Test.wildB_rep2 = c(NA_real_, wild_outlier_responses())
  )
  fit <- suppressWarnings(
    fit_drc_3pl(data = df, normalize = FALSE, verbose = FALSE)
  )
  cmpds <- vapply(fit$detailed_results, function(x) x$compound, character(1))
  ctypes <- vapply(fit$detailed_results, function(x) {
    if (is.null(x$curve_type)) NA_character_ else as.character(x$curve_type)
  }, character(1))
  expect_true(all(ctypes == "flat"))

  out_dir <- tempfile("plmc_badge_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  gg <- suppressWarnings(plot_multiple_compounds(
    results = fit, compound_indices = seq_along(cmpds),
    save_plot = file.path(out_dir, "badge.png"),
    show_display_overrides = TRUE,
    show_display_badge = TRUE,
    verbose = FALSE
  ))
  # All vlines suppressed
  expect_length(vline_xintercepts(gg), 0L)
  # Badge annotation with N/D for both compounds
  labs <- annotation_labels(gg)
  combined <- paste(labs, collapse = " ")
  expect_true(grepl("N/D", combined))
  expect_true(grepl("wildA", combined) && grepl("wildB", combined),
              info = paste("badge text:", combined))
})

# --- Regression against the summary path ------------------------------------

test_that("Fix A: plot vline agrees with summary_table LogIC50", {
  # Cross-check that the vline (post-Fix A) equals the value that would be
  # written into fit$summary_table -- both should come from
  # parameters$Value[LogIC50].
  fit <- fit_one(clean_responses())
  st <- fit$summary_table
  expect_true("LogIC50" %in% names(st))
  li_summary <- as.numeric(st$LogIC50[1])
  expect_true(is.finite(li_summary))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  xi <- vline_xintercepts(gg)
  expect_length(xi, 1L)
  expect_equal(xi[1], li_summary, tolerance = 1e-4)
})

# --- Raw-plateau policy: reported value is always the raw fitted value -------
#
# The biological-plausibility check no longer clamps Bottom/Top plateaus: the
# reported parameters ARE the raw fitted parameters, so the curve drawn from
# them is the raw model curve by construction (the GraphPad-consistent
# behaviour, and the Fix D redraw branch is gone).  Degenerate fits are
# flagged ("Fit diverged - implausible plateau") but still reported raw.
# These tests pin the policy end-to-end:
#   (a) a partial-plateau fit (Bottom ~ 84 on the normalized scale, the
#       A375:LW456 case) reports the raw plateau and both plotters draw the
#       raw model curve -- no dive toward the old clamp value;
#   (b) a diverged degenerate fit is flagged, reports the raw (huge) Bottom,
#       and keeps the flat data-level fallback line.

# Wide-format data frame with control anchor rows (first row response 0, last
# row response 100, both NA log-conc), mirroring the batch_viability_analysis()
# output layout: normalize = TRUE is then an identity transform on a 0-100
# fixture.
build_data_ctrl <- function(responses, compound_name = "Test.Cpd1") {
  stopifnot(length(responses) == 12L)
  df <- data.frame(
    log_conc = c(NA_real_, seq(-4.5, -10, length.out = 12), NA_real_),
    r1       = c(0, responses, 100),
    r2       = c(0, responses, 100)
  )
  names(df) <- c("log_conc",
                 paste0(compound_name, "_rep1"),
                 paste0(compound_name, "_rep2"))
  df
}

# Partial-plateau inhibition curve on the normalized scale (Bottom ~ 84,
# Top ~ 100, LogIC50 ~ -7) -- mirrors A375:LW456.
partial_plateau_responses <- function() {
  b <- 84; t <- 100; li <- -7
  x <- seq(-4.5, -10, length.out = 12)
  b + (t - b) / (1 + 10^(x - li))
}

# First geom_line layer's x/y (the fitted curve), sorted by x.
curve_xy <- function(gg) {
  L <- NULL
  for (ly in gg$layers) {
    if (identical(class(ly$geom)[1], "GeomLine")) { L <- ly; break }
  }
  stopifnot(!is.null(L))
  cx <- if ("x" %in% names(L$data)) L$data$x else L$data$log_inhibitor
  cy <- if ("y" %in% names(L$data)) L$data$y else L$data$response
  ord <- order(cx)
  list(x = cx[ord], y = cy[ord])
}

test_that("raw-plateau policy: partial-plateau fit reports raw Bottom and draws the raw curve", {
  fit <- suppressWarnings(
    fit_drc_3pl(build_data_ctrl(partial_plateau_responses()),
                normalize = TRUE, verbose = FALSE)
  )
  r <- fit$detailed_results[[1]]
  # No correction, no flag: the raw partial plateau (~84) is reported as
  # fitted -- never clamped toward the old viability limit of 60.
  expect_false(isTRUE(r$biological_plausibility_check$needs_correction))
  expect_equal(unname(r$parameters$Value[1]), 84, tolerance = 0.05)
  expect_true(inherits(r$model, "nls"))

  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  cc <- curve_xy(gg)
  raw_at <- as.numeric(stats::predict(
    r$model, newdata = data.frame(log_inhibitor = cc$x)))
  # Drawn curve IS the raw model prediction ...
  expect_lt(max(abs(cc$y - raw_at)), 1e-6)
  # ... following the raw Bottom (~84), never diving toward 60.
  expect_gt(min(cc$y), 75)
})

test_that("raw-plateau policy: diverged fit is flagged, reported raw, keeps flat fallback", {
  fit <- fit_one(wild_outlier_responses())
  r <- fit$detailed_results[[1]]
  bpc <- r$biological_plausibility_check
  # Flag-only: nothing censored, raw diverged Bottom (~1.05e6) reported.
  expect_true(isTRUE(bpc$needs_correction))
  expect_length(bpc$corrections_applied, 0L)
  expect_true(isTRUE(bpc$fit_diverged))
  expect_match(r$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
  expect_gt(abs(unname(r$parameters$Value[1])), 1e5)

  # The degenerate fit keeps the flat data-level fallback line (the raw
  # runaway curve is never drawn for a flat-classified fit).
  gg <- suppressWarnings(plot_dose_response(fit, compound_index = 1,
                                            verbose = FALSE))
  ys <- line_y(gg)
  expect_gt(length(ys), 0L)
  expect_true(all(ys == ys[1]))
  expect_equal(ys[1], stats::median(wild_outlier_responses()), tolerance = 1e-6)
})

test_that("plot_multiple_compounds raw-plateau policy: raw partial plateau, flagged diverged", {
  df <- data.frame(
    log_conc       = c(NA_real_, seq(-4.5, -10, length.out = 12), NA_real_),
    Test.part_rep1 = c(0, partial_plateau_responses(), 100),
    Test.part_rep2 = c(0, partial_plateau_responses(), 100),
    Test.wild_rep1 = c(0, wild_outlier_responses(), 100),
    Test.wild_rep2 = c(0, wild_outlier_responses(), 100)
  )
  fit <- suppressWarnings(fit_drc_3pl(data = df, normalize = TRUE, verbose = FALSE))
  cmpds <- vapply(fit$detailed_results, function(x) x$compound, character(1))
  part_i <- which(grepl("part", cmpds))
  wild_i <- which(grepl("wild", cmpds))
  expect_gt(length(part_i), 0L)
  expect_gt(length(wild_i), 0L)
  # Partial-plateau compound: raw Bottom reported, unflagged.
  r_part <- fit$detailed_results[[part_i[1]]]
  expect_false(isTRUE(r_part$biological_plausibility_check$needs_correction))
  expect_equal(unname(r_part$parameters$Value[1]), 84, tolerance = 0.05)
  # Wild compound: flagged diverged, raw Bottom reported.
  r_wild <- fit$detailed_results[[wild_i[1]]]
  expect_match(r_wild$curve_quality, "Fit diverged - implausible plateau", fixed = TRUE)
  expect_gt(abs(unname(r_wild$parameters$Value[1])), 1e5)

  out_dir <- tempfile("plmc_rp_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  gg <- suppressWarnings(plot_multiple_compounds(
    results = fit, compound_indices = seq_along(cmpds),
    save_plot = file.path(out_dir, "multi_rp.png"), verbose = FALSE
  ))

  line_dat <- NULL
  for (L in gg$layers) {
    if (identical(class(L$geom)[1], "GeomLine")) { line_dat <- L$data; break }
  }
  expect_true(!is.null(line_dat) && nrow(line_dat) > 0)
  col   <- intersect(names(line_dat),
                     c("compound", "compound_name", "Compound", "group"))[1]
  y_col <- intersect(names(line_dat), c("y", "response", "value"))[1]
  part_rows <- line_dat[grepl("part", line_dat[[col]]), ]
  wild_rows <- line_dat[grepl("wild", line_dat[[col]]), ]
  expect_gt(nrow(part_rows), 10L)
  expect_gt(nrow(wild_rows), 10L)
  # Partial-plateau compound draws the raw sigmoid (bottom ~84) ...
  expect_gt(min(part_rows[[y_col]]), 75)
  expect_gt(stats::sd(part_rows[[y_col]]), 1)
  # ... the diverged flat compound keeps the flat data-level line.
  expect_true(all(wild_rows[[y_col]] == wild_rows[[y_col]][1]))
  expect_equal(wild_rows[[y_col]][1], stats::median(wild_outlier_responses()),
               tolerance = 1e-6)
})
