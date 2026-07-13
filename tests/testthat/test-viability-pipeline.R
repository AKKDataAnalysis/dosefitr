# -----------------------------------------------------------------------------
# test-viability-pipeline.R -- end-to-end viability assay integration test
# -----------------------------------------------------------------------------
#
# Exercises the canonical viability route:
#   batch_viability_analysis()  ->  batch_drc_analysis(normalize = TRUE)
#     ->  scarab_viability() (side branch)
#
# Uses the bundled viability_info.xlsx + viability_plate_01/02.xlsx fixtures.
# Assertion strategy mirrors test-nanobret-pipeline.R: shape + type + numeric
# spot-checks against min_good_rsq / min plausible LogIC50.

test_that("batch_viability_analysis runs end-to-end on bundled fixtures", {
  work_dir <- stage_viability_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  via_res <- batch_viability_analysis(
    directory        = work_dir,
    info_file        = "viability_info.xlsx",
    data_pattern     = "viability_plate_\\d+\\.xlsx$",
    control_0perc    = 1,
    control_100perc  = 24,
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempdir(),
    verbose          = FALSE
  )

  # -- shape ----
  expect_type(via_res, "list")
  expect_setequal(names(via_res), c("plate_01", "plate_02"))
  for (plate in names(via_res)) {
    expect_true("result" %in% names(via_res[[plate]]))
    # viability wraps the processed data in $result$modified_ratio_table
    expect_s3_class(via_res[[plate]]$result$modified_ratio_table, "data.frame")
    expect_equal(nrow(via_res[[plate]]$result$modified_ratio_table), 13L)
    # Log-inhibitor axis in first column
    expect_equal(
      colnames(via_res[[plate]]$result$modified_ratio_table)[[1L]],
      "log(inhibitor).[M]"
    )
  }
})

test_that("viability DRC with normalize=TRUE produces plausible bottoms/tops", {
  work_dir <- stage_viability_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  via_res <- batch_viability_analysis(
    directory        = work_dir,
    info_file        = "viability_info.xlsx",
    data_pattern     = "viability_plate_\\d+\\.xlsx$",
    control_0perc    = 1,
    control_100perc  = 24,
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempdir(),
    verbose          = FALSE
  )

  drc_via <- suppressWarnings(batch_drc_analysis(
    batch_results   = via_res,
    model           = "4pl",
    normalize       = TRUE,
    generate_reports = FALSE,
    output_dir      = tempdir(),
    verbose         = FALSE
  ))

  # -- shape ----
  expect_true("drc_results" %in% names(drc_via))
  expect_setequal(names(drc_via$drc_results), c("plate_01", "plate_02"))

  sum_01 <- drc_via$drc_results$plate_01$drc_result$summary_table
  expect_s3_class(sum_01, "data.frame")
  expect_true(all(c("Compound", "Bottom", "Top", "LogIC50", "R_squared",
                    "Curve_Quality") %in% colnames(sum_01)))

  # -- Good-curve rate: >= 80% of the 16 fits should be "Good curve" ----
  good_mask <- sum_01$Curve_Quality == "Good curve"
  good_frac <- sum(good_mask) / nrow(sum_01)
  expect_gte(good_frac, 0.80)

  good_idx <- which(good_mask)
  bottoms  <- suppressWarnings(as.numeric(sum_01$Bottom[good_idx]))
  tops     <- suppressWarnings(as.numeric(sum_01$Top[good_idx]))
  rsq      <- suppressWarnings(as.numeric(sum_01$R_squared[good_idx]))

  # Bottom (saturated inhibition) after normalize=TRUE should sit near 0%.
  # Top (control baseline) should sit near 100%. Give generous windows to
  # avoid brittle regressions.
  expect_true(all(is.finite(bottoms)))
  expect_true(all(is.finite(tops)))
  expect_true(all(bottoms >= -20 & bottoms <= 50),
              info = paste("bottoms:", paste(round(bottoms, 1), collapse = ", ")))
  expect_true(all(tops    >=  60 & tops    <= 130),
              info = paste("tops:", paste(round(tops, 1), collapse = ", ")))

  # R^2 for Good curves must clear min_good_rsq.
  expect_true(all(is.finite(rsq)))
  expect_true(all(rsq >= min_good_rsq))
})

test_that("scarab_viability writes a SCARAB-formatted xlsx", {
  work_dir <- stage_viability_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  via_res <- batch_viability_analysis(
    directory        = work_dir,
    info_file        = "viability_info.xlsx",
    data_pattern     = "viability_plate_\\d+\\.xlsx$",
    control_0perc    = 1,
    control_100perc  = 24,
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempdir(),
    verbose          = FALSE
  )
  drc_via <- suppressWarnings(batch_drc_analysis(
    batch_results   = via_res,
    model           = "4pl",
    normalize       = TRUE,
    generate_reports = FALSE,
    output_dir      = tempdir(),
    verbose         = FALSE
  ))

  outfile <- tempfile(fileext = ".xlsx")
  on.exit(unlink(outfile), add = TRUE)

  scb <- scarab_viability(
    results_list      = via_res,
    drc_results_list  = drc_via,
    plate_name        = "plate_01",
    save              = TRUE,
    file_name         = outfile,
    decimal_separator = "."
  )
  # Confirm the write side effect and a non-null return object.
  expect_true(file.exists(outfile))
  expect_gt(file.info(outfile)$size, 1024L)  # non-trivial xlsx (>1kB)
  expect_true(is.list(scb) || is.data.frame(scb))
})
