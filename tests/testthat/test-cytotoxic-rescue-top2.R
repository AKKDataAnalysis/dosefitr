# -----------------------------------------------------------------------------
# test-cytotoxic-rescue-top2.R
# -----------------------------------------------------------------------------
# Regression test for the keep_cytotoxic rescue rule in rout_outliers_batch().
#
# Documented contract (@param keep_cytotoxic): a ROUT-flagged point is rescued
# ONLY when BOTH conditions hold simultaneously:
#   (1) replicate agreement : both replicates at that concentration are flagged
#   (2) top-2 concentration : the concentration is one of the two highest tested
#
# A prior version enforced only (1), so reproducible drops at ANY concentration
# (including low/mid doses, which are not cytotoxicity) were silently rescued.
# On the bundled fixtures that bug rescued 8 mid-dose points (concs ~ -8 to -9)
# even though the two highest concentrations are ~ -4.6 and -5.1. This test
# locks in that (2) is a HARD filter by driving the public API and asserting
# every rescued point sits at a top-2 concentration.
# -----------------------------------------------------------------------------

# NB: cleanup uses base on.exit(add = TRUE), matching the rest of the suite --
# no external cleanup dependency (e.g. withr) is introduced here.

run_ratio_analysis <- function(work_dir) {
  batch_ratio_analysis(
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
}

top2_concs_of <- function(mrt) {
  concs <- sort(unique(stats::na.omit(suppressWarnings(as.numeric(mrt[[1L]])))),
                decreasing = TRUE)
  concs[seq_len(min(2L, length(concs)))]
}

test_that("every rescued point lies at one of the plate's top-2 concentrations", {
  work_dir <- stage_nanobret_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  ratio_res <- run_ratio_analysis(work_dir)

  # Use a permissive Q so ROUT flags plenty of points to stress the filter.
  res <- suppressWarnings(rout_outliers_batch(
    batch_results  = ratio_res,
    Q              = 0.05,
    keep_cytotoxic = TRUE,
    verbose        = FALSE
  ))

  rs <- res$rescued_summary
  expect_s3_class(rs, "data.frame")

  # Per-plate top-2 concentration sets.
  top2 <- lapply(c("plate_01", "plate_02"), function(pl)
    top2_concs_of(res[[pl]]$result$modified_ratio_table))
  names(top2) <- c("plate_01", "plate_02")

  if (nrow(rs) > 0L) {
    conc_col <- intersect(c("log10_conc", "ln_conc"), names(rs))[1L]
    ok <- vapply(seq_len(nrow(rs)), function(i) {
      any(abs(rs[[conc_col]][i] - top2[[rs$plate[i]]]) < 1e-6)
    }, logical(1))
    # The core invariant: NO rescued point may sit below the top-2 concs.
    expect_true(all(ok),
      info = "A rescued point was NOT at a top-2 concentration (top-2 filter broken).")
  } else {
    succeed("No points met both rescue conditions on this fixture (acceptable).")
  }
})

test_that("keep_cytotoxic = FALSE never rescues (feature is opt-in)", {
  work_dir <- stage_nanobret_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  ratio_res <- run_ratio_analysis(work_dir)

  res <- suppressWarnings(rout_outliers_batch(
    batch_results  = ratio_res,
    Q              = 0.05,
    keep_cytotoxic = FALSE,
    verbose        = FALSE
  ))
  expect_equal(nrow(res$rescued_summary), 0L)
})

test_that("rescued rows carry the cytotoxicity reason label", {
  work_dir <- stage_nanobret_dir()
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  ratio_res <- run_ratio_analysis(work_dir)

  res <- suppressWarnings(rout_outliers_batch(
    batch_results  = ratio_res,
    Q              = 0.05,
    keep_cytotoxic = TRUE,
    verbose        = FALSE
  ))
  rs <- res$rescued_summary
  if (nrow(rs) > 0L) {
    expect_true(all(grepl("top-2 concentration \\(cytotoxicity\\)",
                          rs$rescue_reason)))
  } else {
    succeed("No rescues on this fixture; label check not applicable.")
  }
})
