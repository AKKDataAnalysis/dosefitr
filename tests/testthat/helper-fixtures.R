# ---------------------------------------------------------------------------
# helper-fixtures.R -- shared staging + tolerance helpers for the test suite
# ---------------------------------------------------------------------------
#
# testthat sources every helper-*.R file BEFORE running any test-*.R file, so
# functions and constants defined here are available to all tests without
# being re-declared per-file.  See ?testthat::test_dir.
#
# All fixtures follow the "temp-dir sandbox" pattern: they copy the required
# subset of inst/extdata files into a fresh tempfile()-derived directory,
# leaving the installed package untouched.  Directories are torn down
# automatically when R exits (tempdir() cleanup) but tests also `unlink()`
# aggressively to keep parallel-check runs tidy.

# --- Tolerance constants ---------------------------------------------------
# LogIC50 tolerance (log10 M): 0.5 corresponds to ~3x IC50 wiggle, which is
# tighter than typical inter-plate variance for a well-behaved NanoBRET curve
# and loose enough to survive package-side algorithm nudges.
tol_logic50 <- 0.5

# R^2 tolerance for spot-checking curve quality on the bundled data.
tol_rsq <- 0.05

# Minimum acceptable R^2 for a "Good curve" row.
min_good_rsq <- 0.80

# --- Directory staging -----------------------------------------------------

# Return the inst/extdata directory of the *installed* dosefitr package.
extdata_root <- function() {
  d <- system.file("extdata", package = "dosefitr")
  if (!nzchar(d)) {
    stop("dosefitr is not installed or has no inst/extdata")
  }
  d
}

# Stage a fresh scratch directory containing only the nanobret_* fixtures
# (info file + two 384-well raw data files).  Returns the absolute path to
# the staging directory.
stage_nanobret_dir <- function(prefix = "test_nanobret_") {
  extdata_dir <- extdata_root()
  work_dir <- tempfile(prefix)
  dir.create(work_dir, recursive = TRUE)
  files <- list.files(
    extdata_dir,
    pattern = "^nanobret_",
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("no nanobret_ fixtures found in inst/extdata")
  }
  ok <- file.copy(files, work_dir, overwrite = TRUE)
  if (!all(ok)) {
    stop("failed to stage nanobret fixtures into ", work_dir)
  }
  work_dir
}

# Stage a fresh scratch directory containing only the viability_* fixtures.
stage_viability_dir <- function(prefix = "test_viability_") {
  extdata_dir <- extdata_root()
  work_dir <- tempfile(prefix)
  dir.create(work_dir, recursive = TRUE)
  files <- list.files(
    extdata_dir,
    pattern = "^viability_",
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("no viability_ fixtures found in inst/extdata")
  }
  ok <- file.copy(files, work_dir, overwrite = TRUE)
  if (!all(ok)) {
    stop("failed to stage viability fixtures into ", work_dir)
  }
  work_dir
}

# Fresh scratch subdirectory (for outputs, plots, sheets, etc.).
tmp_out_dir <- function(prefix = "dosefitr_out_") {
  d <- tempfile(prefix)
  dir.create(d, recursive = TRUE)
  d
}

# --- Small helpers ---------------------------------------------------------

# Safely read the first numeric column of a modified_ratio_table (v1) or
# modified_table (v2); used for shape assertions across pipeline stages.
get_modtable <- function(plate_res) {
  # v1 puts the table in $modified_ratio_table; v2 in $modified_table.
  if (!is.null(plate_res$modified_ratio_table)) {
    return(plate_res$modified_ratio_table)
  }
  if (!is.null(plate_res$modified_table)) {
    return(plate_res$modified_table)
  }
  stop("plate result has neither $modified_ratio_table nor $modified_table")
}

# Given a DRC result list (one element per compound), return a data.frame
# with LogIC50 + R_squared + Fit_type columns for easy assertion.
drc_summary_frame <- function(drc_res) {
  if (is.null(drc_res) || length(drc_res) == 0L) {
    return(data.frame(
      Compound = character(0),
      LogIC50 = numeric(0),
      R_squared = numeric(0),
      Fit_type = character(0),
      stringsAsFactors = FALSE
    ))
  }
  cmpd_names <- names(drc_res)
  # DRC output shape varies across the package -- try common slot names.
  extract_scalar <- function(x, key) {
    if (is.null(x)) return(NA_real_)
    if (!is.null(x[[key]])) return(x[[key]][[1L]])
    # Some rows nest scalars in $fit
    if (!is.null(x$fit) && !is.null(x$fit[[key]])) return(x$fit[[key]][[1L]])
    NA_real_
  }
  extract_str <- function(x, key) {
    if (is.null(x)) return(NA_character_)
    v <- x[[key]]
    if (is.null(v)) return(NA_character_)
    as.character(v)[[1L]]
  }
  data.frame(
    Compound  = cmpd_names,
    LogIC50   = vapply(drc_res, extract_scalar, numeric(1L), key = "LogIC50"),
    R_squared = vapply(drc_res, extract_scalar, numeric(1L), key = "R_squared"),
    Fit_type  = vapply(drc_res, extract_str,    character(1L), key = "Fit_type"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
