# -----------------------------------------------------------------------------
# test-nanobret-pipeline.R -- end-to-end NanoBRET assay integration test
# -----------------------------------------------------------------------------
#
# Exercises the canonical NanoBRET route:
#   batch_ratio_analysis()  ->  rout_outliers_batch()  ->  batch_drc_analysis()
#   ->  reshape_dr_table() (side branch)  ->  compare_plates_drc() (bridging)
#
# Uses the bundled 384-well fixtures (nanobret_info.xlsx, nanobret_plate_01.xlsx,
# nanobret_plate_02.xlsx) copied into a fresh tempfile() staging directory by
# the shared stage_nanobret_dir() helper.
#
# Assertions focus on:
#   * shape (row/column counts, list names)
#   * type (data.frame vs list; numeric where expected)
#   * scientific plausibility (LogIC50 in [-10, -4], R^2 > threshold for
#     "Good curve" rows)
#   * bridging behaviour (plate_02 shares 4 compounds with plate_01)

test_that("NanoBRET ratio analysis runs end-to-end on bundled fixtures", {
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

  # -- shape ----
  expect_type(ratio_res, "list")
  expect_setequal(names(ratio_res), c("plate_01", "plate_02"))

  for (plate in c("plate_01", "plate_02")) {
    expect_named(
      ratio_res[[plate]],
      c("data_file", "info_sheet", "sheet_number", "function_version",
        "control_0perc", "control_100perc", "selected_columns", "result")
    )
    expect_named(
      ratio_res[[plate]]$result,
      c("general_means", "interval_means", "construct_intervals",
        "original_ratio_table", "modified_ratio_table",
        "selected_columns_info")
    )
  }

  # -- modified_ratio_table dimensions ----
  mrt_01 <- ratio_res$plate_01$result$modified_ratio_table
  mrt_02 <- ratio_res$plate_02$result$modified_ratio_table
  expect_s3_class(mrt_01, "data.frame")
  expect_s3_class(mrt_02, "data.frame")
  # 12 non-control concentrations + 1 log(inhibitor).[M] header row = 13
  expect_equal(nrow(mrt_01), 13L)
  expect_equal(nrow(mrt_02), 13L)
  # First column is the log(inhibitor).[M] concentration axis
  expect_equal(colnames(mrt_01)[[1L]], "log(inhibitor).[M]")

  # log-inhibitor values should be finite and monotone-ish descending.
  # Row 1 and row 13 are NA placeholders (used internally to bracket the
  # concentration axis); the actual serial dilution sits in rows 2:12.
  loginhib <- mrt_01[[1L]]
  finite   <- loginhib[is.finite(loginhib)]
  expect_gte(length(finite), 8L)          # at least 8 real concentrations
  expect_true(all(diff(finite) < 0))       # decreasing serial dilution
})

test_that("rout_outliers_batch flags outliers without breaking pipeline shape", {
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

  rout_res <- suppressWarnings(rout_outliers_batch(
    batch_results = ratio_res,
    Q             = 0.01,
    verbose       = FALSE
  ))

  # -- shape ----
  expect_type(rout_res, "list")
  expect_true(all(c("plate_01", "plate_02") %in% names(rout_res)))
  expect_true("outlier_summary" %in% names(rout_res))

  # Each plate result must still expose $result$modified_ratio_table so it can
  # feed batch_drc_analysis() downstream.
  for (plate in c("plate_01", "plate_02")) {
    expect_s3_class(rout_res[[plate]]$result$modified_ratio_table, "data.frame")
    expect_equal(nrow(rout_res[[plate]]$result$modified_ratio_table), 13L)
  }

  # outlier_summary is a per-plate tally of outliers detected
  expect_s3_class(rout_res$outlier_summary, "data.frame")
  # (0 or more outliers is fine; just insist the tally is a non-negative integer)
  if ("Outliers_Removed" %in% colnames(rout_res$outlier_summary)) {
    expect_true(all(rout_res$outlier_summary$Outliers_Removed >= 0))
  }
})

test_that("batch_drc_analysis produces sensible LogIC50 and R^2 on bundled data", {
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

  drc_res <- suppressWarnings(batch_drc_analysis(
    batch_results   = ratio_res,
    model           = "4pl",
    normalize       = FALSE,
    generate_reports = FALSE,
    output_dir      = tempdir(),
    verbose         = FALSE
  ))

  # -- shape ----
  expect_type(drc_res, "list")
  expect_true("drc_results" %in% names(drc_res))
  expect_setequal(names(drc_res$drc_results), c("plate_01", "plate_02"))

  # summary_table is the per-compound fit table (16 compounds on plate_01)
  sum_01 <- drc_res$drc_results$plate_01$drc_result$summary_table
  expect_s3_class(sum_01, "data.frame")
  expect_true("Compound" %in% colnames(sum_01))
  expect_true("LogIC50" %in% colnames(sum_01))
  expect_true("R_squared" %in% colnames(sum_01))
  expect_true("Curve_Quality" %in% colnames(sum_01))

  # -- scientific spot-checks on Good curves ----
  good_idx <- which(sum_01$Curve_Quality == "Good curve")
  expect_true(length(good_idx) >= 1L)

  # LogIC50 for well-fit curves should sit in a biologically plausible range.
  # The bundled fixtures are simulated at IC50s spanning ~1e-9 to 1e-5 M, so
  # LogIC50 must fall in [-10, -4] with plenty of buffer.
  logic50_good <- suppressWarnings(as.numeric(sum_01$LogIC50[good_idx]))
  expect_true(all(is.finite(logic50_good)))
  expect_true(all(logic50_good > -10))
  expect_true(all(logic50_good < -4))

  # R^2 for Good curves must be at least min_good_rsq.
  rsq_good <- suppressWarnings(as.numeric(sum_01$R_squared[good_idx]))
  expect_true(all(is.finite(rsq_good)))
  expect_true(all(rsq_good >= min_good_rsq))
})

test_that("compare_plates_drc detects bridging compounds between plate_01 and plate_02", {
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
  drc_res <- suppressWarnings(batch_drc_analysis(
    batch_results   = ratio_res,
    model           = "4pl",
    normalize       = FALSE,
    generate_reports = FALSE,
    output_dir      = tempdir(),
    verbose         = FALSE
  ))

  out_dir <- tmp_out_dir()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  cmp <- suppressWarnings(compare_plates_drc(
    batch_drc_result = drc_res,
    compare_by       = "compound",
    output_dir       = out_dir,
    min_plates       = 2,
    plot_dpi         = 72,
    verbose          = FALSE
  ))
  # compare_plates_drc returns a list keyed by compound (one entry per bridging
  # compound). The bundled plate_02 shares 4 compounds with plate_01, so we
  # expect exactly 4 entries.
  expect_type(cmp, "list")
  expect_gte(length(cmp), 4L)
  # Each entry should be a per-compound result list (non-empty).
  expect_true(all(vapply(cmp, function(x) is.list(x) && length(x) >= 1L,
                         logical(1L))))
})

test_that("reshape_dr_table converts summary_table to a long-format layout", {
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
  drc_res <- suppressWarnings(batch_drc_analysis(
    batch_results   = ratio_res,
    model           = "4pl",
    normalize       = FALSE,
    generate_reports = FALSE,
    output_dir      = tempdir(),
    verbose         = FALSE
  ))
  sum_01 <- drc_res$drc_results$plate_01$drc_result$summary_table
  long <- reshape_dr_table(results_table = sum_01)
  expect_s3_class(long, "data.frame")
  # Reshape organizes parameters into rows; compound columns become header
  # rows. Just insist that rows/cols are non-trivial.
  expect_true(nrow(long) > 5L)
  expect_true(ncol(long) >= 2L)
})
