# -----------------------------------------------------------------------------
# test-batch_read_tables.R -- round-trip tests for the pre-computed table
# importer.
# -----------------------------------------------------------------------------
#
# batch_read_tables() lets a user skip the raw-signal ratio step (typical when
# the normalisation was done in Prism or in a legacy script) and jump straight
# into batch_drc_analysis().  This file exercises:
#
#   1. Round-trip nanobret: batch_ratio_analysis() -> write xlsx ->
#      batch_read_tables() -> identical modified_ratio_table.
#   2. Object attributes (assay_source, qc_available).
#   3. End-to-end equivalence: LogIC50 from a batch_drc_analysis() on the
#      re-imported plate matches the LogIC50 from the original ratio result
#      within 1e-6.
#   4. Multi-plate: two xlsx files -> two plates named plate_01, plate_02.
#   5. Bad header (no colon) -> informative error naming both file and
#      offending column.
#   6. Invalid assay_source string ("elisa") -> error.
#   7. Default assay_source (no arg) -> "nanobret".
#   8. Default log_conc_col_name = NULL -> read column 1 whatever it is
#      called (numeric check still enforced).
#   9. Explicit log_conc_col_name pointing to a mid-sheet column ->
#      reorder so the log-conc column becomes column 1, keeping the
#      remaining columns in their original order.
#  10. Explicit log_conc_col_name not found in the file -> error naming
#      the missing column.
#  11. Viability round-trip: batch_viability_analysis() -> xlsx ->
#      batch_read_tables(assay_source = "viability") -> batch_drc_analysis()
#      still succeeds.
#  12. file_map = character vector: reads exactly the listed files, in
#      the given order, and IGNORES file_pattern completely (even when a
#      pattern that would match nothing is passed).
#  13. file_map = named list (batch_ratio_analysis style): list names
#      are silently ignored; plate identifiers still come from the
#      default plate_01, plate_02 sequence.
#  14. file_map referencing a missing file -> informative error naming
#      the missing path.
#  15. Invalid file_map shapes (non-character, empty, duplicate) -> error.
#
# Uses the bundled fixtures via the shared stage_*_dir() helpers in
# helper-fixtures.R.

# --- Setup helpers -----------------------------------------------------------

# Write a data.frame to a single-sheet xlsx, preserving column names verbatim.
write_ratio_table <- function(df, path) {
  openxlsx::write.xlsx(df, file = path, colNames = TRUE, rowNames = FALSE)
}

# Build a NanoBRET ratio result on the bundled fixtures, silently.
build_nanobret_ratio <- function() {
  work_dir <- stage_nanobret_dir("brt_nb_")
  ratio_res <- batch_ratio_analysis(
    directory        = work_dir,
    info_file        = "nanobret_info.xlsx",
    data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc    = "1",
    control_100perc  = "24",
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempfile("brt_nb_out_"),
    verbose          = FALSE
  )
  list(work_dir = work_dir, ratio_res = ratio_res)
}

# Build a viability ratio result on the bundled fixtures, silently.
build_viability_ratio <- function() {
  via_dir <- stage_viability_dir("brt_vb_")
  via_res <- batch_viability_analysis(
    directory        = via_dir,
    info_file        = "viability_info.xlsx",
    data_pattern     = "viability_plate_\\d+\\.xlsx$",
    control_0perc    = 13,
    control_100perc  = 12,
    selected_columns = 2:23,
    generate_reports = FALSE,
    output_dir       = tempfile("brt_vb_out_"),
    verbose          = FALSE
  )
  list(work_dir = via_dir, via_res = via_res)
}

# --- Tests -------------------------------------------------------------------

test_that("nanobret round-trip preserves modified_ratio_table", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)

  # Export plate_01 modified_ratio_table to xlsx
  shipped <- tempfile("brt_shipped_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  src_tbl <- bundle$ratio_res$plate_01$result$modified_ratio_table
  write_ratio_table(src_tbl, file.path(shipped, "plate_01.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(shipped, assay_source = "nanobret", verbose = FALSE)
  )

  # Basic shape
  expect_type(imported, "list")
  expect_named(imported, "plate_01")

  # Column names identical
  ret_tbl <- imported$plate_01$result$modified_ratio_table
  expect_equal(colnames(ret_tbl), colnames(src_tbl))

  # Values equal within numeric tolerance (xlsx round-trip is not bit-exact)
  expect_equal(ret_tbl, src_tbl, tolerance = 1e-8,
               ignore_attr = TRUE)

  # interval_means MUST be NULL (QC gone by construction)
  expect_null(imported$plate_01$result$interval_means)
})

test_that("nanobret result carries assay_source and qc_available attributes", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_attrs_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)
  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(shipped, assay_source = "nanobret", verbose = FALSE)
  )
  expect_identical(attr(imported, "assay_source"), "nanobret")
  expect_identical(attr(imported, "qc_available"), FALSE)
  expect_s3_class(imported, "dosefitr_batch_result")
})

test_that("batch_drc_analysis produces matching LogIC50s from imported plate", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_drc_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)
  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(shipped, assay_source = "nanobret", verbose = FALSE)
  )

  # Restrict the original ratio_res to plate_01 so both DRC runs act on the
  # exact same plate.
  ratio_only_01 <- bundle$ratio_res["plate_01"]
  attr(ratio_only_01, "assay_source") <- "nanobret"

  drc_orig <- batch_drc_analysis(ratio_only_01, model = "3pl",
                                 generate_reports = FALSE, verbose = FALSE)
  drc_imp  <- batch_drc_analysis(imported,      model = "3pl",
                                 generate_reports = FALSE, verbose = FALSE)

  # The tidy per-compound IC50 table lives at
  #   drc_results$<plate>$drc_result$summary_table
  # with columns Compound, LogIC50, R_squared, Curve_Quality, ...
  st_orig <- drc_orig$drc_results$plate_01$drc_result$summary_table
  st_imp  <- drc_imp $drc_results$plate_01$drc_result$summary_table

  expect_s3_class(st_orig, "data.frame")
  expect_s3_class(st_imp,  "data.frame")

  # Same compounds, same set (order-independent).
  expect_setequal(st_orig$Compound, st_imp$Compound)

  # LogIC50 must agree to numerical precision for compounds that fit
  # successfully on both branches (the imported plate carries the identical
  # ratio table, so every fit should reproduce).
  ord_imp <- st_imp[match(st_orig$Compound, st_imp$Compound), , drop = FALSE]
  finite_both <- is.finite(st_orig$LogIC50) & is.finite(ord_imp$LogIC50)
  expect_true(any(finite_both))
  expect_equal(st_orig$LogIC50[finite_both],
               ord_imp$LogIC50[finite_both],
               tolerance = 1e-6)
})

test_that("multi-plate import produces plate_01, plate_02 in sorted order", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_multi_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))
  write_ratio_table(bundle$ratio_res$plate_02$result$modified_ratio_table,
                    file.path(shipped, "plate_02.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(shipped, assay_source = "nanobret", verbose = FALSE)
  )
  expect_equal(names(imported), c("plate_01", "plate_02"))

  # Each plate carries a data.frame at the expected slot
  for (p in names(imported)) {
    expect_s3_class(imported[[p]]$result$modified_ratio_table, "data.frame")
  }
})

test_that("bad column header raises an informative error", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_bad_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  src_tbl <- bundle$ratio_res$plate_01$result$modified_ratio_table
  # Sabotage: rename column 2 to one without a colon.
  colnames(src_tbl)[[2L]] <- "LRRK2_Cpd1"

  bad_path <- file.path(shipped, "plate_01.xlsx")
  write_ratio_table(src_tbl, bad_path)

  # Error must name the file and the offending column.
  expect_error(
    batch_read_tables(shipped, assay_source = "nanobret", verbose = FALSE),
    regexp = "plate_01\\.xlsx.*LRRK2_Cpd1"
  )
})

test_that("invalid assay_source string is rejected", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_asy_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))
  expect_error(
    batch_read_tables(shipped, assay_source = "elisa", verbose = FALSE),
    regexp = "assay_source"
  )
})

test_that("assay_source defaults to 'nanobret' when omitted", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_def_asy_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))
  # No assay_source arg
  imported <- suppressMessages(
    batch_read_tables(shipped, verbose = FALSE)
  )
  expect_identical(attr(imported, "assay_source"), "nanobret")
})

test_that("log_conc_col_name = NULL accepts a first column of any name", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_lcn_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  src_tbl <- bundle$ratio_res$plate_01$result$modified_ratio_table
  # Rename column 1 to something non-canonical.
  colnames(src_tbl)[[1L]] <- "log_M"
  write_ratio_table(src_tbl, file.path(shipped, "plate_01.xlsx"))

  # Default log_conc_col_name = NULL: take column 1 regardless of name.
  imported <- suppressMessages(
    batch_read_tables(shipped, verbose = FALSE)
  )
  ret_tbl <- imported$plate_01$result$modified_ratio_table
  expect_identical(colnames(ret_tbl)[[1L]], "log_M")
  # Column 1 must still be numeric with same values (round-trip via xlsx).
  expect_type(ret_tbl[[1L]], "double")
  expect_equal(ret_tbl[[1L]], src_tbl[[1L]], tolerance = 1e-8)
})

test_that("explicit log_conc_col_name finds a mid-sheet column and moves it to column 1", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_reord_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  src_tbl <- bundle$ratio_res$plate_01$result$modified_ratio_table
  # Put a compound column first and shuffle the log-conc column to
  # position 5 (a mid-sheet location).
  target_idx <- 5L
  reorder    <- c(2:target_idx, 1L,
                  seq_len(ncol(src_tbl))[-c(1L, 2:target_idx)])
  shuffled   <- src_tbl[, reorder, drop = FALSE]

  # Sanity: log(inhibitor).[M] is at column `target_idx` in `shuffled`.
  stopifnot(colnames(shuffled)[[target_idx]] == "log(inhibitor).[M]")

  write_ratio_table(shuffled, file.path(shipped, "plate_01.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(shipped,
                      log_conc_col_name = "log(inhibitor).[M]",
                      verbose           = FALSE)
  )
  ret_tbl <- imported$plate_01$result$modified_ratio_table

  # After reorder: column 1 is the log-conc column, and the remaining
  # columns keep the shuffled order (with the moved-out column removed).
  expect_identical(colnames(ret_tbl)[[1L]], "log(inhibitor).[M]")
  expect_type(ret_tbl[[1L]], "double")

  # Numerical equality of the log-conc values with the source.
  expect_equal(ret_tbl[[1L]], src_tbl[[1L]], tolerance = 1e-8)

  # Compound columns must be the same set, just reordered.
  expect_setequal(colnames(ret_tbl)[-1L], colnames(src_tbl)[-1L])
})

test_that("explicit log_conc_col_name not found in file is an error", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_notfound_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))

  expect_error(
    batch_read_tables(shipped,
                      log_conc_col_name = "not_a_real_column",
                      verbose           = FALSE),
    regexp = "not_a_real_column"
  )
})

test_that("viability round-trip still feeds batch_drc_analysis", {
  bundle <- build_viability_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_vb_ship_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  src_tbl <- bundle$via_res$plate_01$result$modified_ratio_table
  write_ratio_table(src_tbl, file.path(shipped, "plate_01.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(shipped, assay_source = "viability", verbose = FALSE)
  )
  expect_identical(attr(imported, "assay_source"), "viability")

  drc_imp <- suppressWarnings(
    batch_drc_analysis(imported, model = "3pl",
                       generate_reports = FALSE, verbose = FALSE)
  )
  expect_type(drc_imp, "list")
  expect_true("drc_results" %in% names(drc_imp))
  expect_setequal(names(drc_imp$drc_results), "plate_01")
})

test_that("file_map = character vector reads exactly the listed files and ignores file_pattern", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_fmap_char_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  # Write TWO plates but under NON-standard filenames that would not be
  # matched by the default "\\.xlsx$" pattern in a decoy directory, and
  # ALSO write a decoy xlsx that would match the pattern but is NOT in
  # the map.  Reading with an explicit file_map must ignore the decoy.
  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "run_A.xlsx"))
  write_ratio_table(bundle$ratio_res$plate_02$result$modified_ratio_table,
                    file.path(shipped, "run_B.xlsx"))
  # Decoy: also a .xlsx, would be picked up by auto-discovery
  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "decoy_ignore_me.xlsx"))

  # Pass a file_pattern that MATCHES NOTHING to prove it is ignored
  imported <- suppressMessages(
    batch_read_tables(
      directory    = shipped,
      assay_source = "nanobret",
      file_pattern = "^__nomatch__$",
      file_map     = c("run_A.xlsx", "run_B.xlsx"),
      verbose      = FALSE
    )
  )
  expect_equal(names(imported), c("plate_01", "plate_02"))
  # Order preserved from file_map (NOT alphabetical): plate_01 = run_A
  expect_true(grepl("run_A\\.xlsx$", imported$plate_01$data_file))
  expect_true(grepl("run_B\\.xlsx$", imported$plate_02$data_file))
  # Decoy file was NOT read
  data_files <- vapply(imported, function(x) basename(x$data_file), character(1))
  expect_false(any(data_files == "decoy_ignore_me.xlsx"))
})

test_that("file_map = named list ignores names, uses default plate_01, plate_02 order", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_fmap_list_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "alpha.xlsx"))
  write_ratio_table(bundle$ratio_res$plate_02$result$modified_ratio_table,
                    file.path(shipped, "beta.xlsx"))

  imported <- suppressMessages(
    batch_read_tables(
      directory    = shipped,
      assay_source = "nanobret",
      file_map     = list("weird_name_1" = "alpha.xlsx",
                          "weird_name_2" = "beta.xlsx"),
      verbose      = FALSE
    )
  )
  # Names from file_map are dropped
  expect_equal(names(imported), c("plate_01", "plate_02"))
  # Order preserved from list
  expect_true(grepl("alpha\\.xlsx$", imported$plate_01$data_file))
  expect_true(grepl("beta\\.xlsx$",  imported$plate_02$data_file))
})

test_that("file_map referencing a non-existent file raises an error", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_fmap_miss_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)

  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "real.xlsx"))

  expect_error(
    batch_read_tables(
      directory    = shipped,
      assay_source = "nanobret",
      file_map     = c("real.xlsx", "does_not_exist.xlsx"),
      verbose      = FALSE
    ),
    regexp = "does_not_exist\\.xlsx"
  )
})

test_that("invalid file_map shapes are rejected", {
  bundle <- build_nanobret_ratio()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  shipped <- tempfile("brt_fmap_bad_"); dir.create(shipped)
  on.exit(unlink(shipped, recursive = TRUE), add = TRUE)
  write_ratio_table(bundle$ratio_res$plate_01$result$modified_ratio_table,
                    file.path(shipped, "plate_01.xlsx"))

  # Numeric file_map
  expect_error(
    batch_read_tables(shipped, file_map = 1:3, verbose = FALSE),
    regexp = "file_map"
  )
  # Empty file_map
  expect_error(
    batch_read_tables(shipped, file_map = character(0), verbose = FALSE),
    regexp = "file_map"
  )
  # Duplicate entries
  expect_error(
    batch_read_tables(shipped,
                      file_map = c("plate_01.xlsx", "plate_01.xlsx"),
                      verbose = FALSE),
    regexp = "duplicate"
  )
  # List with non-scalar entries
  expect_error(
    batch_read_tables(shipped,
                      file_map = list(c("a.xlsx", "b.xlsx")),
                      verbose = FALSE),
    regexp = "file_map"
  )
})
