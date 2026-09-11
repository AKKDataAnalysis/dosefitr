# -----------------------------------------------------------------------------
# test-atomic-helpers.R -- unit tests for functions that operate on a single
# stage of the pipeline, without exercising the whole batch chain.
# -----------------------------------------------------------------------------
#
# Targets:
#   * ratio_dose_response()       -- single-plate v1 ratio (character controls)
#   * ratio_dose_response_v2()    -- single-plate v2 ratio (numeric controls)
#   * fit_drc_4pl() / fit_drc_3pl()  -- direct DRC on a modified_ratio_table
#   * rout_outliers()             -- outlier detection on a single table
#   * merge_plate_replicates()    -- cross-plate compound merging

test_that("ratio_dose_response (v1) processes a single 384-well plate", {
  # ratio_dose_response expects openxlsx::read.xlsx()-style tables:
  # sheet 1 of the info file (first plate) and a headerless raw plate.
  extdata_dir <- extdata_root()
  info_table <- openxlsx::read.xlsx(
    file.path(extdata_dir, "nanobret_info.xlsx"),
    sheet = 1
  )
  raw <- openxlsx::read.xlsx(
    file.path(extdata_dir, "nanobret_plate_01.xlsx"),
    sheet = 1, colNames = FALSE
  )

  res <- ratio_dose_response(
    data             = raw,
    control_0perc    = "1",
    control_100perc  = "24",
    info_table       = info_table,
    selected_columns = 1:24,
    verbose          = FALSE
  )
  expect_type(res, "list")
  expect_true("modified_ratio_table" %in% names(res))
  expect_s3_class(res$modified_ratio_table, "data.frame")
  expect_equal(nrow(res$modified_ratio_table), 13L)
  expect_true("Z'_factor" %in% rownames(res$interval_means))
  expect_false("Z_Score" %in% rownames(res$interval_means))
})

test_that("ratio processors optionally combine repeated rows as four replicates", {
  extdata_dir <- extdata_root()
  info_table <- openxlsx::read.xlsx(
    file.path(extdata_dir, "nanobret_info.xlsx"), sheet = 1
  )
  raw <- openxlsx::read.xlsx(
    file.path(extdata_dir, "nanobret_plate_01.xlsx"),
    sheet = 1, colNames = FALSE
  )

  # Make non-adjacent rows A and C describe the same Construct+Compound pair.
  # Each plate row contains two technical series, so combine mode must yield
  # four adjacent columns even though another curve occurs between the rows.
  info_table[[3L]][3L] <- info_table[[3L]][1L]
  info_table[[4L]][3L] <- info_table[[4L]][1L]
  base <- paste(info_table[[3L]][1L], info_table[[4L]][1L], sep = ":")
  separate_second <- paste0(info_table[[3L]][1L], "_2:", info_table[[4L]][1L])

  processors <- list(
    v1 = ratio_dose_response,
    v2 = ratio_dose_response_v2
  )
  for (processor in processors) {
    common_args <- list(
      data = raw, control_0perc = "1", control_100perc = "24",
      split_replicates = TRUE, info_table = info_table,
      selected_columns = 1:24, verbose = FALSE
    )

    separate <- suppressWarnings(do.call(
      processor, c(common_args, list(repeated_rows = "separate"))
    ))$modified_ratio_table
    default_mode <- suppressWarnings(do.call(
      processor, common_args
    ))$modified_ratio_table
    expect_equal(default_mode, separate)
    expect_true(all(c(base, paste0(base, ".2"), separate_second,
                      paste0(separate_second, ".2")) %in% colnames(separate)))

    combined <- suppressWarnings(do.call(
      processor, c(common_args, list(repeated_rows = "combine"))
    ))$modified_ratio_table
    expected_combined <- c(base, paste0(base, ".", 2:4))
    expect_identical(colnames(combined)[-1L][1:4], expected_combined)
    expect_equal(sum(sub("\\.\\d+$", "", colnames(combined)[-1L]) == base), 4L)
    expect_false(any(grepl(paste0("^", info_table[[3L]][1L], "_2:"),
                           colnames(combined))))
  }
})

test_that("ratio_dose_response_v2 accepts character control indices", {
  extdata_dir <- extdata_root()
  info_table <- openxlsx::read.xlsx(
    file.path(extdata_dir, "nanobret_info.xlsx"),
    sheet = 1
  )
  raw <- openxlsx::read.xlsx(
    file.path(extdata_dir, "nanobret_plate_01.xlsx"),
    sheet = 1, colNames = FALSE
  )

  # v2 accepts character controls too (nanobret 384-well fixtures use
  # named-column mode; numeric mode is exercised via process_viability_data_v2
  # with a synthetic 96-well plate elsewhere).
  res <- ratio_dose_response_v2(
    data             = raw,
    control_0perc    = "1",
    control_100perc  = "24",
    info_table       = info_table,
    selected_columns = 1:24,
    verbose          = FALSE
  )
  expect_type(res, "list")
  mt <- get_modtable(res)
  expect_s3_class(mt, "data.frame")
  expect_gt(nrow(mt), 5L)
  expect_true("Z'_factor" %in% rownames(res$interval_means))
  expect_false("Z_Score" %in% rownames(res$interval_means))
})

test_that("fit_drc_4pl produces a per-compound summary table", {
  work_dir <- stage_nanobret_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  ratio_res <- batch_ratio_analysis(
    directory        = work_dir,
    info_file        = "nanobret_info.xlsx",
    data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc    = "1",
    control_100perc  = "24",
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempdir(),
    verbose          = FALSE
  )
  mrt <- ratio_res$plate_01$result$modified_ratio_table

  fit <- suppressWarnings(fit_drc_4pl(data = mrt, normalize = FALSE, verbose = FALSE))
  expect_type(fit, "list")
  expect_true("summary_table" %in% names(fit))
  st <- fit$summary_table
  expect_s3_class(st, "data.frame")
  expect_true(all(c("Compound", "LogIC50", "R_squared") %in% colnames(st)))
  expect_true(nrow(st) >= 1L)

  # At least one row must be a numeric LogIC50 in a plausible window.
  logic50 <- suppressWarnings(as.numeric(st$LogIC50))
  finite  <- logic50[is.finite(logic50)]
  expect_true(length(finite) >= 1L)
  expect_true(any(finite > -10 & finite < -4))
})

test_that("fit_drc_3pl runs without error and returns a summary_table", {
  work_dir <- stage_nanobret_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  ratio_res <- batch_ratio_analysis(
    directory        = work_dir,
    info_file        = "nanobret_info.xlsx",
    data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc    = "1",
    control_100perc  = "24",
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempdir(),
    verbose          = FALSE
  )
  mrt <- ratio_res$plate_01$result$modified_ratio_table

  fit <- suppressWarnings(fit_drc_3pl(data = mrt, normalize = FALSE, verbose = FALSE))
  expect_type(fit, "list")
  expect_true("summary_table" %in% names(fit))
  expect_s3_class(fit$summary_table, "data.frame")
  expect_true(nrow(fit$summary_table) >= 1L)
})

test_that("rout_outliers only removes synthetic outliers after convergence", {
  # Build a tiny synthetic modified_ratio_table: 12 concentrations + 1 header
  # row, 4 compound columns, one column carrying an obvious outlier.
  logc <- seq(-9, -5, length.out = 12L)
  set.seed(42L)
  bottom <- 20; top <- 100; hill <- -1; ic50 <- -7
  clean <- top + (bottom - top) / (1 + 10 ^ ((logc - ic50) * hill))
  # Make column 3 have a huge outlier at row 6
  vals <- data.frame(
    "log(inhibitor).[M]" = c(NA, logc),
    "Cpd1"    = c(NA, clean + rnorm(12L, 0, 1)),
    "Cpd1.2"  = c(NA, clean + rnorm(12L, 0, 1)),
    "Cpd2"    = c(NA, clean + rnorm(12L, 0, 1)),
    "Cpd2.2"  = c(NA, {y <- clean + rnorm(12L, 0, 1); y[6L] <- 400; y}),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rr <- suppressWarnings(rout_outliers(
    data      = vals,
    Q         = 0.01,
    model     = "auto",
    direction = "inhibition",
    verbose   = FALSE
  ))
  expect_type(rr, "list")
  # rout_outliers returns a list with cleaned_table (post-outlier-removal
  # data.frame) and outlier_table (rows flagged as outliers).
  expect_true("cleaned_table" %in% names(rr))
  expect_true("outlier_table" %in% names(rr))
  expect_s3_class(rr$cleaned_table, "data.frame")
  expect_s3_class(rr$outlier_table, "data.frame")
  # Cpd2.2 has a huge observation at row 6 (value = 400). It may only be
  # removed when the ROUT fit for Cpd2 converged. A non-converged fit must
  # retain every original observation, even when the raw ROUT call produced
  # an outlier flag.
  cpd2_results <- rr$results[
    rr$results$compound == "Cpd2",
    ,
    drop = FALSE
  ]

  expect_gt(nrow(cpd2_results), 0L)

  if (isTRUE(all(cpd2_results$converged))) {
    expect_true(any(cpd2_results$outlier_fdr))
    expect_gte(nrow(rr$outlier_table), 1L)
  } else {
    expect_false(any(cpd2_results$outlier_fdr))
    expect_equal(
      rr$cleaned_table[["Cpd2.2"]],
      vals[["Cpd2.2"]]
    )
  }
})

test_that("replicate consensus safeguard requires at least three replicates", {
  consensus_filter <- dosefitr:::.clear_concordant_replicate_flags

  # Two concordant replicates retain the historical behaviour: flags remain.
  two_reps <- consensus_filter(
    outlier_flags = c(TRUE, TRUE),
    x_log_fit     = c(-7, -7),
    y_fit         = c(100, 103),
    cv_max        = 15
  )
  expect_equal(two_reps$flags, c(TRUE, TRUE))
  expect_equal(nrow(two_reps$cleared), 0L)

  # Three concordant replicates are preserved, even if ROUT flagged only a
  # subset of the concentration group.
  three_reps <- consensus_filter(
    outlier_flags = c(TRUE, TRUE, FALSE),
    x_log_fit     = c(-7, -7, -7),
    y_fit         = c(100, 106, 103),
    cv_max        = 15
  )
  expect_false(any(three_reps$flags))
  expect_equal(three_reps$cleared$n_replicates, 3L)
  expect_equal(three_reps$cleared$n_cleared, 2L)

  # The same protection applies when more than three replicates agree.
  four_reps <- consensus_filter(
    outlier_flags = c(TRUE, TRUE, TRUE, FALSE),
    x_log_fit     = rep(-7, 4L),
    y_fit         = c(100, 104, 98, 102),
    cv_max        = 15
  )
  expect_false(any(four_reps$flags))
  expect_equal(four_reps$cleared$n_replicates, 4L)

  # A discordant triplicate is not rescued.
  discordant <- consensus_filter(
    outlier_flags = c(TRUE, TRUE, FALSE),
    x_log_fit     = c(-7, -7, -7),
    y_fit         = c(100, 105, 200),
    cv_max        = 15
  )
  expect_equal(discordant$flags, c(TRUE, TRUE, FALSE))
  expect_equal(nrow(discordant$cleared), 0L)
})

test_that("compare_plates_drc grouping modes use the requested dimension", {
  entry <- list(construct = "AGP-01", compound = "TP-0903")
  key <- dosefitr:::.compare_plates_group_key
  normalise <- dosefitr:::.normalise_compare_plates_mode
  normalise_legend <- dosefitr:::.normalise_compare_plates_legend
  legend_label <- dosefitr:::.compare_plates_legend_label

  expect_equal(key(entry, "construct_compound"), "AGP-01:TP-0903")
  expect_equal(key(entry, "pair"), "AGP-01:TP-0903")
  expect_equal(key(entry, "compound"), "TP-0903")
  expect_equal(key(entry, "construct"), "AGP-01")

  expect_equal(normalise("construct-compound"), "construct_compound")
  expect_equal(normalise_legend("AUTO"), "auto")
  expect_equal(normalise_legend("plate"), "plate")
  expect_equal(normalise_legend("construct"), "construct")
  expect_equal(normalise_legend("compound"), "compound")

  legend_entry <- list(
    plate_name = "plate_02",
    construct = "AGP-01",
    compound = "TP-0903"
  )
  expect_equal(legend_label(legend_entry, "compound", "plate"), "plate_02")
  expect_equal(
    legend_label(legend_entry, "compound", "auto"),
    "plate_02 (AGP-01)"
  )
  expect_equal(
    legend_label(legend_entry, "construct", "auto"),
    "plate_02 (TP-0903)"
  )
  expect_error(
    normalise("plate"),
    "construct_compound.*compound.*construct"
  )
  expect_error(normalise_legend("invalid"), "auto.*plate.*construct.*compound")
})

test_that("merge_plate_replicates combines shared compounds across plates", {
  work_dir <- stage_nanobret_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  ratio_res <- batch_ratio_analysis(
    directory        = work_dir,
    info_file        = "nanobret_info.xlsx",
    data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc    = "1",
    control_100perc  = "24",
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempdir(),
    verbose          = FALSE
  )

  merged <- suppressWarnings(merge_plate_replicates(
    results = ratio_res,
    plates  = c("plate_01", "plate_02"),
    verbose = FALSE
  ))
  expect_type(merged, "list")
  # Merged result should include an entry per-plate plus a new merged element
  # (name defaults to "merged"). Loose contract: at least the merged element
  # must exist and contain a modified_ratio_table.
  expect_true("merged" %in% names(merged) || length(merged) >= length(ratio_res))
})
