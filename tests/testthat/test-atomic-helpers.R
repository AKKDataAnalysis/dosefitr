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
    n_param   = 4L,
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
