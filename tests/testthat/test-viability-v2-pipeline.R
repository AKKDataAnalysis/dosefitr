# -----------------------------------------------------------------------------
# test-viability-v2-pipeline.R -- unit + integration tests for the v2 viability
# processor (process_viability_data_v2) and its batch wiring.
# -----------------------------------------------------------------------------
#
# v2 ports the NanoBRET control style (ratio_dose_response_v2) into viability:
#   * 0% control  = a fixed scalar -> constant first row (Fixed_0perc)
#   * 100% control = plate column(s), averaged into the last row (Mean_100perc)
#                    and then consumed (removed) from the experimental block
#   * control_mean_scope in {row, construct, global} controls how the 100%
#     columns are averaged.
#
# Unlike the other viability tests, these are FULLY SELF-CONTAINED: they build
# synthetic 96-well and 384-well plates in-line and never touch inst/extdata,
# so they run in any environment (including ones without the bundled fixtures).

# ---------------------------------------------------------------------------
# Local synthetic-data builders (kept inside the test file on purpose so they
# do not collide with the fixture-based helpers in helper-fixtures.R).
# ---------------------------------------------------------------------------

# Pad a short descriptor vector with NA up to length n (info_table columns must
# share a common length; the concentration column is usually the longest).
.v2_pad <- function(x, n) c(x, rep(NA, n - length(x)))

# Build a synthetic BMG-style raw plate: `n_hdr` blank header rows, a column-
# label row ("", 1..n_cols), then one data row per plate-row letter. The LAST
# data column is a saturated (~100%) control; the remaining data columns follow
# a sigmoid in log-concentration space.
.v2_make_plate <- function(n_rows, n_cols, logc,
                           ec50 = -7, top = 800, bot = 40, seed = 1L,
                           noise = 4) {
  set.seed(seed)
  n_exp <- n_cols - 1L
  stopifnot(length(logc) == n_exp)
  n_hdr <- 12L
  header_rows   <- as.data.frame(matrix("", nrow = n_hdr, ncol = n_cols + 1L),
                                 stringsAsFactors = FALSE)
  col_label_row <- as.data.frame(matrix(c("", as.character(seq_len(n_cols))),
                                        nrow = 1L), stringsAsFactors = FALSE)
  data_rows     <- as.data.frame(matrix("", nrow = n_rows, ncol = n_cols + 1L),
                                 stringsAsFactors = FALSE)
  for (i in seq_len(n_rows)) {
    frac <- 1 / (1 + 10^((logc - ec50)))
    vals <- bot + (top - bot) * frac + rnorm(n_exp, 0, noise)
    data_rows[i, 1L]                <- LETTERS[i]
    data_rows[i, 2L:(n_exp + 1L)]   <- as.character(round(vals))
    data_rows[i, n_cols + 1L]       <- as.character(round(top + rnorm(1L, 0, 3)))
  }
  raw <- rbind(header_rows, col_label_row, data_rows)
  colnames(raw) <- paste0("V", seq_len(ncol(raw)))
  raw
}

# Build a matching info_table. `constructs`/`compounds` are per-plate-row
# vectors (length n_rows); the first column holds the log-concentration axis
# c(NA, logc, NA) and shorter columns are NA-padded to that length.
.v2_make_info <- function(n_rows, logc, constructs, compounds,
                          first_col_name = "log(inhibitor).[M]") {
  conc <- c(NA, logc, NA)
  L <- length(conc)
  info <- data.frame(
    conc      = conc,
    Plate_Row = .v2_pad(LETTERS[seq_len(n_rows)], L),
    Construct = .v2_pad(constructs, L),
    Compound  = .v2_pad(compounds, L),
    stringsAsFactors = FALSE
  )
  colnames(info)[1] <- first_col_name
  info
}

# Convenience: a canonical 96-well plate (cols 1-12, control col 12, 11 exp
# wells) with 2 constructs (A-D, E-H) + matching info_table.
.v2_fixture_96 <- function(seed = 1L) {
  logc <- seq(-9, -5, length.out = 11L)
  list(
    raw  = .v2_make_plate(8L, 12L, logc, seed = seed),
    info = .v2_make_info(8L, logc,
                         constructs = c(rep("KinaseA", 4L), rep("KinaseB", 4L)),
                         compounds  = c(rep("Cpd1", 4L),   rep("Cpd2", 4L))),
    logc = logc
  )
}

# Convenience: a canonical 384-well plate (cols 1-24, control col 24, 23 exp
# wells) with 4 constructs of 4 rows each + matching info_table.
.v2_fixture_384 <- function(seed = 2L) {
  logc <- seq(-9, -5, length.out = 23L)
  list(
    raw  = .v2_make_plate(16L, 24L, logc, seed = seed),
    info = .v2_make_info(16L, logc,
                         constructs = rep(paste0("K", 1:4), each = 4L),
                         compounds  = rep(paste0("C", 1:4), each = 4L)),
    logc = logc
  )
}

# Write a synthetic multi-plate directory for batch tests. Returns the dir path.
.v2_write_batch_dir <- function(n_plates = 2L, n_rows = 8L, n_cols = 12L) {
  logc <- seq(-9, -5, length.out = n_cols - 1L)
  info <- .v2_make_info(
    n_rows, logc,
    constructs = if (n_rows == 8L) c(rep("KinaseA", 4L), rep("KinaseB", 4L))
                 else rep(paste0("K", seq_len(n_rows %/% 4L)), each = 4L),
    compounds  = if (n_rows == 8L) c(rep("Cpd1", 4L), rep("Cpd2", 4L))
                 else rep(paste0("C", seq_len(n_rows %/% 4L)), each = 4L)
  )
  dir <- tempfile("v2_batch_")
  dir.create(dir, recursive = TRUE)
  info_wb <- openxlsx::createWorkbook()
  for (p in seq_len(n_plates)) {
    sheet <- paste0("info_", p)
    openxlsx::addWorksheet(info_wb, sheet)
    openxlsx::writeData(info_wb, sheet, info)
  }
  openxlsx::saveWorkbook(info_wb, file.path(dir, "info_tables.xlsx"), overwrite = TRUE)
  for (p in seq_len(n_plates)) {
    openxlsx::write.xlsx(
      .v2_make_plate(n_rows, n_cols, logc, seed = p),
      file.path(dir, sprintf("viability_plate_%d.xlsx", p)),
      colNames = FALSE
    )
  }
  list(dir = dir, info = info, logc = logc)
}

# ===========================================================================
# UNIT TESTS -- process_viability_data_v2()
# ===========================================================================

test_that("v2 returns the documented list shape and metadata", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  expect_type(out, "list")
  expect_true(all(c("original_table", "modified_table", "processing_info",
                    "selected_columns_info", "version", "data_type",
                    "auto_detect", "behavior_mode") %in% names(out)))
  expect_identical(out$version, "v2")
  expect_identical(out$data_type, "viability")
  expect_identical(out$behavior_mode, "v2-fixed0-colmean100")
  expect_s3_class(out$modified_table, "data.frame")
})

test_that("v2 96-well shape: 7 rows, Fixed_0perc first, Mean_100perc last", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  # k=12, m=1 -> exp=11 -> floor(11/2)=5 -> nrow = 5 + 2 = 7
  expect_equal(nrow(out$modified_table), 7L)
  expect_identical(rownames(out$modified_table)[1], "Fixed_0perc")
  expect_identical(rownames(out$modified_table)[nrow(out$modified_table)], "Mean_100perc")
  # First column is the (verbatim) info_table concentration column name.
  expect_identical(colnames(out$modified_table)[1], "log(inhibitor).[M]")
  # First row (excluding the leading NA conc cell) is the fixed 0% scalar.
  first_row <- unlist(out$modified_table[1, -1], use.names = FALSE)
  expect_true(all(first_row == 0))
  # Leading conc cell of the two control rows is NA.
  expect_true(is.na(out$modified_table[1, 1]))
  expect_true(is.na(out$modified_table[nrow(out$modified_table), 1]))
})

test_that("v2 splits technical replicates with .2-suffixed columns", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  data_cols <- colnames(out$modified_table)[-1]
  # Every biological construct column should have a paired .2 replicate.
  base <- sub("\\.2$", "", data_cols)
  rep2 <- grepl("\\.2$", data_cols)
  expect_true(any(rep2))
  # For each base id, both the rep1 and rep2 columns exist.
  for (b in unique(base)) {
    expect_true(b %in% data_cols)
    expect_true(paste0(b, ".2") %in% data_cols)
  }
})

test_that("v2 no-split keeps a single set of construct columns and 13 rows", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", split_replicates = FALSE, verbose = FALSE
  )
  # No split: nrow = (k - m) + 2 = 11 + 2 = 13
  expect_equal(nrow(out$modified_table), 13L)
  expect_false(any(grepl("\\.2$", colnames(out$modified_table)[-1])))
})

test_that("v2 384-well plate is detected and yields 13 rows", {
  fx <- .v2_fixture_384()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 24,
    info_table = fx$info, selected_columns = 1:24,
    control_mean_scope = "row", verbose = FALSE
  )
  # k=24, m=1 -> exp=23 -> floor(23/2)=11 -> nrow = 11 + 2 = 13
  expect_equal(nrow(out$modified_table), 13L)
  expect_identical(out$version, "v2")
})

test_that("v2 accepts numeric and character control_100perc identically", {
  fx <- .v2_fixture_96()
  out_num <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  out_chr <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = "12",
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  expect_equal(out_num$modified_table, out_chr$modified_table)
})

test_that("v2 consumes the 100% control column(s) from the experimental block", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  # original_table keeps all 12 selected columns...
  expect_equal(ncol(out$original_table), 12L)
  # ...but the experimental rows in modified_table exclude the control column:
  # 7 total rows - 2 control rows = 5 experimental rows (floor(11/2)).
  expect_equal(nrow(out$modified_table) - 2L, 5L)
  # control_100_info records which column(s) were consumed.
  expect_identical(out$processing_info$control_100_info$name, "12")
  expect_identical(out$processing_info$control_100_info$scope, "row")
})

test_that("v2 scope='global' gives one constant Mean_100perc across constructs", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "global", verbose = FALSE
  )
  last_row <- unlist(out$modified_table[nrow(out$modified_table), -1], use.names = FALSE)
  last_row <- last_row[!is.na(last_row)]
  expect_gt(length(last_row), 0L)
  # All construct columns share a single global mean.
  expect_equal(length(unique(round(last_row, 6))), 1L)
})

test_that("v2 scope='row' varies Mean_100perc per plate row", {
  # Use a plate where each row's 100% control differs substantially so the
  # per-row means are distinguishable.
  logc <- seq(-9, -5, length.out = 11L)
  raw <- .v2_make_plate(8L, 12L, logc, seed = 5L)
  # Overwrite the control column (col 12 -> data col 13) with row-specific values.
  ctrl_vals <- seq(700, 910, length.out = 8L)
  for (i in seq_len(8L)) raw[13L + i - 1L, 13L] <- as.character(round(ctrl_vals[i]))
  info <- .v2_make_info(8L, logc,
                        constructs = paste0("K", 1:8),   # 8 distinct constructs
                        compounds  = paste0("C", 1:8))
  out <- process_viability_data_v2(
    data = raw, control_0perc = 0, control_100perc = 12,
    info_table = info, selected_columns = 1:12,
    control_mean_scope = "row", verbose = FALSE
  )
  last_row <- unlist(out$modified_table[nrow(out$modified_table), -1], use.names = FALSE)
  last_row <- last_row[!is.na(last_row)]
  # Row scope should produce more than one distinct 100% value here.
  expect_gt(length(unique(round(last_row, 3))), 1L)
})

test_that("v2 scope='construct' is constant within a construct's columns", {
  fx <- .v2_fixture_96()  # KinaseA = A-D, KinaseB = E-H
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 0, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "construct", verbose = FALSE
  )
  last <- out$modified_table[nrow(out$modified_table), -1, drop = FALSE]
  cols <- colnames(last)
  # Group columns by their construct prefix (strip compound + .2 suffix).
  construct_of <- sub(":.*$", "", cols)
  vals <- unlist(last, use.names = FALSE)
  for (cn in unique(construct_of)) {
    v <- vals[construct_of == cn]
    v <- v[!is.na(v)]
    if (length(v) > 1L)
      expect_equal(length(unique(round(v, 6))), 1L,
                   info = paste("construct", cn, "should share one value"))
  }
})

test_that("v2 validation errors fire for bad inputs", {
  fx <- .v2_fixture_96()
  # 0% must be a single finite scalar
  expect_error(
    process_viability_data_v2(fx$raw, control_0perc = c(0, 1),
      control_100perc = 12, info_table = fx$info, verbose = FALSE),
    "single numeric value"
  )
  # 100% column index out of range
  expect_error(
    process_viability_data_v2(fx$raw, control_0perc = 0,
      control_100perc = 99, info_table = fx$info, verbose = FALSE),
    "out of range"
  )
  # info_table required
  expect_error(
    process_viability_data_v2(fx$raw, control_0perc = 0,
      control_100perc = 12, info_table = NULL, verbose = FALSE),
    "requires an info_table"
  )
  # every 100% column must be in selected_columns
  expect_error(
    process_viability_data_v2(fx$raw, control_0perc = 0,
      control_100perc = 12, info_table = fx$info,
      selected_columns = 1:11, verbose = FALSE),
    "must be included in selected_columns"
  )
  # bad control_mean_scope
  expect_error(
    process_viability_data_v2(fx$raw, control_0perc = 0,
      control_100perc = 12, info_table = fx$info,
      control_mean_scope = "nonsense", verbose = FALSE)
  )
})

test_that("v2 processing_info records the fixed 0% value and scope", {
  fx <- .v2_fixture_96()
  out <- process_viability_data_v2(
    data = fx$raw, control_0perc = 7.5, control_100perc = 12,
    info_table = fx$info, selected_columns = 1:12,
    control_mean_scope = "construct", verbose = FALSE
  )
  c0 <- out$processing_info$control_0_info
  expect_identical(c0$name, "Fixed_0perc")
  expect_false(c0$is_relative)
  expect_equal(c0$fixed_value, 7.5)
  expect_identical(out$processing_info$control_mean_scope, "construct")
  # The fixed value propagates to the first row.
  first_row <- unlist(out$modified_table[1, -1], use.names = FALSE)
  expect_true(all(first_row == 7.5))
})

# ===========================================================================
# UNIT TESTS -- control_mean_scope on the ROW-BASED processors
#   process_viability_data_v3 (row-based; both controls are plate columns) and
#   process_viability_data   (v1; 384-well only). Here "row" is an IDENTITY
#   (each plate row keeps its own control values); "construct" replaces each
#   Target group with its mean; "global" uses one plate-wide mean per control
#   column. The scope is applied JOINTLY to the 0% and 100% control columns and
#   is gated by apply_control_means.
# ===========================================================================

# Build a row-based plate with KNOWN 0% (col 1) and 100% (last col) control
# values so the per-scope aggregation can be checked exactly. Experimental
# wells sit in the interior columns.
.rb_make_plate <- function(n_rows, n_cols, logc, ctrl0_vals, ctrl100_vals,
                           ec50 = -7, seed = 11L) {
  set.seed(seed)
  n_exp <- n_cols - 2L
  stopifnot(length(logc) == n_exp,
            length(ctrl0_vals) == n_rows, length(ctrl100_vals) == n_rows)
  n_hdr <- 12L
  header_rows   <- as.data.frame(matrix("", nrow = n_hdr, ncol = n_cols + 1L),
                                 stringsAsFactors = FALSE)
  col_label_row <- as.data.frame(matrix(c("", as.character(seq_len(n_cols))),
                                        nrow = 1L), stringsAsFactors = FALSE)
  data_rows     <- as.data.frame(matrix("", nrow = n_rows, ncol = n_cols + 1L),
                                 stringsAsFactors = FALSE)
  for (i in seq_len(n_rows)) {
    frac <- 1 / (1 + 10^((logc - ec50)))
    expv <- 40 + (800 - 40) * frac + rnorm(n_exp, 0, 2)
    data_rows[i, 1L]                 <- LETTERS[i]           # plate-row label
    data_rows[i, 2L]                 <- as.character(ctrl0_vals[i])   # col 1 = 0%
    data_rows[i, 3L:(n_exp + 2L)]    <- as.character(round(expv))     # interior
    data_rows[i, n_cols + 1L]        <- as.character(ctrl100_vals[i]) # last col = 100%
  }
  raw <- rbind(header_rows, col_label_row, data_rows)
  colnames(raw) <- paste0("V", seq_len(ncol(raw)))
  raw
}

.rb_make_info <- function(n_rows, logc, targets, compounds,
                          first_col_name = "log(inhibitor).[M]") {
  conc <- c(NA, logc, NA)   # length = 1 ctrl0 + n_exp + 1 ctrl100
  L <- length(conc)
  info <- data.frame(
    conc      = conc,
    Plate_Row = .v2_pad(LETTERS[seq_len(n_rows)], L),
    Target    = .v2_pad(targets, L),
    Compound  = .v2_pad(compounds, L),
    stringsAsFactors = FALSE
  )
  colnames(info)[1] <- first_col_name
  info
}

# Pull the two control columns out of viability_modified_with_means, indexed by
# their (character) column names, ordered by plate-row letter.
.rb_controls <- function(out, c0_name, c100_name) {
  vm <- out$processing_info$viability_modified_with_means
  ord <- order(rownames(vm))
  list(c0 = vm[ord, c0_name], c100 = vm[ord, c100_name],
       rows = rownames(vm)[ord])
}

test_that("v3 (row-based) scope=construct replaces each Target with its mean", {
  logc <- seq(-9, -5, length.out = 10L)              # 10 interior wells
  c0   <- c(10, 20, 30, 40, 50, 60, 70, 80)          # A..H
  c100 <- c(900, 910, 920, 930, 940, 950, 960, 970)
  raw  <- .rb_make_plate(8L, 12L, logc, c0, c100)
  info <- .rb_make_info(8L, logc,
                        targets   = c(rep("KinaseA", 4L), rep("KinaseB", 4L)),
                        compounds = c(rep("Cpd1", 4L),    rep("Cpd2", 4L)))
  out <- process_viability_data_v3(
    data = raw, control_0perc = 1, control_100perc = 12,
    info_table = info, selected_columns = 1:12,
    apply_control_means = TRUE, control_mean_scope = "construct", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "12")
  expect_equal(unname(ctl$c0),   c(rep(25, 4), rep(65, 4)))      # per-Target 0% mean
  expect_equal(unname(ctl$c100), c(rep(915, 4), rep(955, 4)))    # per-Target 100% mean
})

test_that("v3 (row-based) scope=global uses one plate-wide mean per control column", {
  logc <- seq(-9, -5, length.out = 10L)
  c0   <- c(10, 20, 30, 40, 50, 60, 70, 80)
  c100 <- c(900, 910, 920, 930, 940, 950, 960, 970)
  raw  <- .rb_make_plate(8L, 12L, logc, c0, c100)
  info <- .rb_make_info(8L, logc,
                        targets   = c(rep("KinaseA", 4L), rep("KinaseB", 4L)),
                        compounds = c(rep("Cpd1", 4L),    rep("Cpd2", 4L)))
  out <- process_viability_data_v3(
    data = raw, control_0perc = 1, control_100perc = 12,
    info_table = info, selected_columns = 1:12,
    apply_control_means = TRUE, control_mean_scope = "global", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "12")
  expect_equal(unname(ctl$c0),   rep(45, 8))      # mean(10..80)
  expect_equal(unname(ctl$c100), rep(935, 8))     # mean(900..970)
})

test_that("v3 (row-based) scope=row is identity (each row keeps its own controls)", {
  logc <- seq(-9, -5, length.out = 10L)
  c0   <- c(10, 20, 30, 40, 50, 60, 70, 80)
  c100 <- c(900, 910, 920, 930, 940, 950, 960, 970)
  raw  <- .rb_make_plate(8L, 12L, logc, c0, c100)
  info <- .rb_make_info(8L, logc,
                        targets   = c(rep("KinaseA", 4L), rep("KinaseB", 4L)),
                        compounds = c(rep("Cpd1", 4L),    rep("Cpd2", 4L)))
  out <- process_viability_data_v3(
    data = raw, control_0perc = 1, control_100perc = 12,
    info_table = info, selected_columns = 1:12,
    apply_control_means = TRUE, control_mean_scope = "row", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "12")
  expect_equal(unname(ctl$c0),   c0)      # raw per-row values, untouched
  expect_equal(unname(ctl$c100), c100)
})

test_that("v3 (row-based) apply_control_means=FALSE leaves controls raw (scope ignored)", {
  logc <- seq(-9, -5, length.out = 10L)
  c0   <- c(10, 20, 30, 40, 50, 60, 70, 80)
  c100 <- c(900, 910, 920, 930, 940, 950, 960, 970)
  raw  <- .rb_make_plate(8L, 12L, logc, c0, c100)
  info <- .rb_make_info(8L, logc,
                        targets   = c(rep("KinaseA", 4L), rep("KinaseB", 4L)),
                        compounds = c(rep("Cpd1", 4L),    rep("Cpd2", 4L)))
  out <- process_viability_data_v3(
    data = raw, control_0perc = 1, control_100perc = 12,
    info_table = info, selected_columns = 1:12,
    apply_control_means = FALSE, control_mean_scope = "construct", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "12")
  expect_equal(unname(ctl$c0),   c0)      # gate off -> no replacement
  expect_equal(unname(ctl$c100), c100)
})

test_that("v3 (row-based) records control_mean_scope in processing_info", {
  logc <- seq(-9, -5, length.out = 10L)
  c0   <- c(10, 20, 30, 40, 50, 60, 70, 80)
  c100 <- c(900, 910, 920, 930, 940, 950, 960, 970)
  raw  <- .rb_make_plate(8L, 12L, logc, c0, c100)
  info <- .rb_make_info(8L, logc,
                        targets   = c(rep("KinaseA", 4L), rep("KinaseB", 4L)),
                        compounds = c(rep("Cpd1", 4L),    rep("Cpd2", 4L)))
  for (sc in c("construct", "row", "global")) {
    out <- process_viability_data_v3(
      data = raw, control_0perc = 1, control_100perc = 12,
      info_table = info, selected_columns = 1:12,
      apply_control_means = TRUE, control_mean_scope = sc, verbose = FALSE)
    expect_identical(out$processing_info$control_mean_scope, sc)
    expect_identical(out$version, "v3")
  }
})


# ===========================================================================
# UNIT TESTS -- control_mean_scope on v1 (process_viability_data, 384-well)
#   v1 accepts both apply_control_means (outer gate) and control_mean_scope.
#   Default scope is "construct" (legacy per-construct means). "row" is an
#   identity (one control well per row -> nothing to average within a row).
# ===========================================================================

test_that("v1 scope=construct replaces each construct with its per-group mean", {
  logc <- seq(-9, -5, length.out = 22L)              # 22 interior wells (cols 2-23)
  c0   <- seq(10, by = 5, length.out = 16L)          # A..P, distinct per row
  c100 <- seq(900, by = 3, length.out = 16L)
  raw  <- .rb_make_plate(16L, 24L, logc, c0, c100, seed = 21L)
  # 4 constructs of 4 rows each: rows A-D, E-H, I-L, M-P.
  info <- .rb_make_info(16L, logc,
                        targets   = rep(paste0("K", 1:4), each = 4L),
                        compounds = paste0("C", 1:16))  # distinct compound per row -> Target defines the construct group
  out <- process_viability_data(
    data = raw, control_0perc = 1, control_100perc = 24,
    info_table = info, selected_columns = 1:24,
    apply_control_means = TRUE, control_mean_scope = "construct", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "24")
  # Expected per-construct means (groups of 4 consecutive rows).
  exp0   <- rep(tapply(c0,   rep(1:4, each = 4L), mean), each = 4L)
  exp100 <- rep(tapply(c100, rep(1:4, each = 4L), mean), each = 4L)
  expect_equal(unname(ctl$c0),   unname(exp0))
  expect_equal(unname(ctl$c100), unname(exp100))
})

test_that("v1 scope=global uses one plate-wide mean per control column", {
  logc <- seq(-9, -5, length.out = 22L)
  c0   <- seq(10, by = 5, length.out = 16L)
  c100 <- seq(900, by = 3, length.out = 16L)
  raw  <- .rb_make_plate(16L, 24L, logc, c0, c100, seed = 22L)
  info <- .rb_make_info(16L, logc,
                        targets   = rep(paste0("K", 1:4), each = 4L),
                        compounds = paste0("C", 1:16))  # distinct compound per row -> Target defines the construct group
  out <- process_viability_data(
    data = raw, control_0perc = 1, control_100perc = 24,
    info_table = info, selected_columns = 1:24,
    apply_control_means = TRUE, control_mean_scope = "global", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "24")
  expect_equal(unname(ctl$c0),   rep(mean(c0),   16L))
  expect_equal(unname(ctl$c100), rep(mean(c100), 16L))
})

test_that("v1 scope=row is identity (each row keeps its own controls)", {
  logc <- seq(-9, -5, length.out = 22L)
  c0   <- seq(10, by = 5, length.out = 16L)
  c100 <- seq(900, by = 3, length.out = 16L)
  raw  <- .rb_make_plate(16L, 24L, logc, c0, c100, seed = 23L)
  info <- .rb_make_info(16L, logc,
                        targets   = rep(paste0("K", 1:4), each = 4L),
                        compounds = paste0("C", 1:16))  # distinct compound per row -> Target defines the construct group
  out <- process_viability_data(
    data = raw, control_0perc = 1, control_100perc = 24,
    info_table = info, selected_columns = 1:24,
    apply_control_means = TRUE, control_mean_scope = "row", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "24")
  expect_equal(unname(ctl$c0),   c0)
  expect_equal(unname(ctl$c100), c100)
})

test_that("v1 apply_control_means=FALSE leaves controls raw (scope ignored)", {
  logc <- seq(-9, -5, length.out = 22L)
  c0   <- seq(10, by = 5, length.out = 16L)
  c100 <- seq(900, by = 3, length.out = 16L)
  raw  <- .rb_make_plate(16L, 24L, logc, c0, c100, seed = 24L)
  info <- .rb_make_info(16L, logc,
                        targets   = rep(paste0("K", 1:4), each = 4L),
                        compounds = paste0("C", 1:16))  # distinct compound per row -> Target defines the construct group
  out <- process_viability_data(
    data = raw, control_0perc = 1, control_100perc = 24,
    info_table = info, selected_columns = 1:24,
    apply_control_means = FALSE, control_mean_scope = "construct", verbose = FALSE)
  ctl <- .rb_controls(out, "1", "24")
  expect_equal(unname(ctl$c0),   c0)
  expect_equal(unname(ctl$c100), c100)
})

test_that("v1 records control_mean_scope in processing_info", {
  logc <- seq(-9, -5, length.out = 22L)
  c0   <- seq(10, by = 5, length.out = 16L)
  c100 <- seq(900, by = 3, length.out = 16L)
  raw  <- .rb_make_plate(16L, 24L, logc, c0, c100, seed = 25L)
  info <- .rb_make_info(16L, logc,
                        targets   = rep(paste0("K", 1:4), each = 4L),
                        compounds = paste0("C", 1:16))  # distinct compound per row -> Target defines the construct group
  for (sc in c("construct", "row", "global")) {
    out <- process_viability_data(
      data = raw, control_0perc = 1, control_100perc = 24,
      info_table = info, selected_columns = 1:24,
      apply_control_means = TRUE, control_mean_scope = sc, verbose = FALSE)
    expect_identical(out$processing_info$control_mean_scope, sc)
  }
})

# ===========================================================================
# INTEGRATION TESTS -- batch_viability_analysis(version = "v2")
# ===========================================================================

test_that("batch v2 runs end-to-end on a synthetic 2-plate directory", {
  bd <- .v2_write_batch_dir(n_plates = 2L, n_rows = 8L, n_cols = 12L)
  on.exit(unlink(bd$dir, recursive = TRUE), add = TRUE)
  out_dir <- tmp_out_dir()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res <- suppressWarnings(batch_viability_analysis(
    directory        = bd$dir,
    control_0perc    = 0,
    control_100perc  = 12,
    info_file        = "info_tables.xlsx",
    data_pattern     = "_\\d+\\.xlsx$",
    output_dir       = out_dir,
    version          = "v2",
    generate_reports = TRUE,
    verbose          = FALSE
  ))

  expect_type(res, "list")
  expect_equal(length(res), 2L)
  expect_identical(attr(res, "assay_source"), "viability")
  for (nm in names(res)) {
    r <- res[[nm]]$result
    expect_identical(r$version, "v2")
    # batch renames modified_table -> modified_ratio_table
    expect_s3_class(r$modified_ratio_table, "data.frame")
    expect_null(r$modified_table)
    expect_equal(nrow(r$modified_ratio_table), 7L)
  }
})

test_that("batch v2 stores control_mean_scope and defaults it to 'row'", {
  bd <- .v2_write_batch_dir(2L, 8L, 12L)
  on.exit(unlink(bd$dir, recursive = TRUE), add = TRUE)
  out_dir <- tmp_out_dir(); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # No control_mean_scope supplied -> v2 default 'row'.
  res_default <- suppressWarnings(batch_viability_analysis(
    directory = bd$dir, control_0perc = 0, control_100perc = 12,
    info_file = "info_tables.xlsx", data_pattern = "_\\d+\\.xlsx$",
    output_dir = out_dir, version = "v2", generate_reports = FALSE, verbose = FALSE
  ))
  expect_identical(
    res_default[[1]]$result$processing_info$control_mean_scope, "row"
  )

  # Explicit 'construct' is honoured.
  res_con <- suppressWarnings(batch_viability_analysis(
    directory = bd$dir, control_0perc = 0, control_100perc = 12,
    info_file = "info_tables.xlsx", data_pattern = "_\\d+\\.xlsx$",
    output_dir = out_dir, version = "v2", control_mean_scope = "construct",
    generate_reports = FALSE, verbose = FALSE
  ))
  expect_identical(
    res_con[[1]]$result$processing_info$control_mean_scope, "construct"
  )
})

test_that("batch v2 QC: fixed 0% -> Mean_Background scalar, SD=0, CV=NA", {
  bd <- .v2_write_batch_dir(1L, 8L, 12L)
  on.exit(unlink(bd$dir, recursive = TRUE), add = TRUE)
  out_dir <- tmp_out_dir(); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  suppressWarnings(batch_viability_analysis(
    directory = bd$dir, control_0perc = 0, control_100perc = 12,
    info_file = "info_tables.xlsx", data_pattern = "_\\d+\\.xlsx$",
    output_dir = out_dir, version = "v2", generate_reports = TRUE, verbose = FALSE
  ))

  qc_file <- file.path(out_dir, "drc_quality", "viability_results_1.xlsx")
  expect_true(file.exists(qc_file))
  sheets <- openxlsx::getSheetNames(qc_file)
  expect_true("Quality_Metrics" %in% sheets)

  qm <- openxlsx::read.xlsx(qc_file, sheet = "Quality_Metrics")
  # Metrics are laid out with a 'Metric' column and one column per construct.
  get_metric <- function(metric) {
    row <- qm[qm$Metric == metric, , drop = FALSE]
    unlist(row[, setdiff(colnames(qm), "Metric")], use.names = FALSE)
  }
  mean_bg <- suppressWarnings(as.numeric(get_metric("Mean_Background")))
  sd_bg   <- suppressWarnings(as.numeric(get_metric("SD_Background")))
  cv_bg   <- get_metric("CV_Background_pct")

  expect_true(all(mean_bg == 0))              # fixed 0% scalar
  expect_true(all(sd_bg == 0))                # constant -> SD 0
  # CV should be NA (blank / non-numeric) for every construct.
  expect_true(all(is.na(suppressWarnings(as.numeric(cv_bg)))))
})

test_that("batch v2 output is compatible with batch_drc_analysis(normalize=TRUE)", {
  skip_if_not_installed("OptimModel")
  bd <- .v2_write_batch_dir(1L, 8L, 12L)
  on.exit(unlink(bd$dir, recursive = TRUE), add = TRUE)
  out_dir <- tmp_out_dir(); on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  via_res <- suppressWarnings(batch_viability_analysis(
    directory = bd$dir, control_0perc = 0, control_100perc = 12,
    info_file = "info_tables.xlsx", data_pattern = "_\\d+\\.xlsx$",
    output_dir = out_dir, version = "v2", generate_reports = FALSE, verbose = FALSE
  ))

  drc <- suppressWarnings(batch_drc_analysis(
    batch_results    = via_res,
    model            = "4pl",
    normalize        = TRUE,
    generate_reports = FALSE,
    output_dir       = out_dir,
    verbose          = FALSE
  ))

  expect_true("drc_results" %in% names(drc))
  expect_equal(length(drc$drc_results), 1L)
  sum_tbl <- drc$drc_results[[1]]$drc_result$summary_table
  expect_s3_class(sum_tbl, "data.frame")
  expect_true(all(c("Compound", "Bottom", "Top", "LogIC50", "R_squared") %in%
                    colnames(sum_tbl)))
  # At least one curve should fit (finite LogIC50) on this clean synthetic data.
  logic50 <- suppressWarnings(as.numeric(sum_tbl$LogIC50))
  expect_true(any(is.finite(logic50)))
})

test_that("batch: v1/v3 scope='row' is a valid identity choice (no fallback warning)", {
  # B1: on the row-based processors a single control well per row means "row"
  # scope is an identity (no cross-row averaging) -- it must NOT warn or fall
  # back. The scope bridge runs before directory validation, so point at a
  # nonexistent directory and capture any messages emitted before the error.
  capture_msgs <- function(version) {
    msgs <- character(0)
    withCallingHandlers(
      tryCatch(
        batch_viability_analysis(
          directory = tempfile("does_not_exist_"),
          control_0perc = 1, control_100perc = 12,
          version = version, control_mean_scope = "row",
          apply_control_means = TRUE, verbose = TRUE
        ),
        error = function(e) invisible(NULL)
      ),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
    )
    msgs
  }
  for (v in c("v1", "v3")) {
    msgs <- capture_msgs(v)
    # No legacy "only supported by version" fallback message any more.
    expect_false(any(grepl("only supported by version", msgs)))
    # And no "scope is ignored" message, because apply_control_means = TRUE.
    expect_false(any(grepl("scope is ignored", msgs)))
  }
})

test_that("batch: B2 message fires when apply_control_means=FALSE + explicit non-construct scope", {
  # B2: if the user turns the control-mean gate off but still sets a non-default
  # scope explicitly, batch emits a one-line message that the scope is ignored.
  capture_msgs <- function(version, scope, apply_means, explicit = TRUE) {
    msgs <- character(0)
    args <- list(
      directory = tempfile("does_not_exist_"),
      control_0perc = 1, control_100perc = 12,
      version = version, apply_control_means = apply_means, verbose = TRUE
    )
    if (explicit) args$control_mean_scope <- scope
    withCallingHandlers(
      tryCatch(do.call(batch_viability_analysis, args),
               error = function(e) invisible(NULL)),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
    )
    msgs
  }
  # v1 + explicit "global" + gate OFF -> message fires.
  msgs <- capture_msgs("v1", "global", apply_means = FALSE)
  expect_true(any(grepl("scope is ignored", msgs)))
  # v3 + explicit "row" + gate OFF -> message fires.
  msgs <- capture_msgs("v3", "row", apply_means = FALSE)
  expect_true(any(grepl("scope is ignored", msgs)))
  # Gate OFF but scope left at default (construct) -> NO message.
  msgs <- capture_msgs("v1", "construct", apply_means = FALSE)
  expect_false(any(grepl("scope is ignored", msgs)))
  # Gate OFF, scope not supplied at all -> NO message.
  msgs <- capture_msgs("v3", NULL, apply_means = FALSE, explicit = FALSE)
  expect_false(any(grepl("scope is ignored", msgs)))
})

test_that("batch v2 rejects a non-finite fixed 0% value", {
  bd <- .v2_write_batch_dir(1L, 8L, 12L)
  on.exit(unlink(bd$dir, recursive = TRUE), add = TRUE)
  expect_error(
    batch_viability_analysis(
      directory = bd$dir, control_0perc = Inf, control_100perc = 12,
      info_file = "info_tables.xlsx", data_pattern = "_\\d+\\.xlsx$",
      version = "v2", generate_reports = FALSE, verbose = FALSE
    ),
    "finite numeric value"
  )
})

test_that("process_viability_data_v2 is exported in NAMESPACE", {
  # Guards against the export being dropped; the function must be visible.
  expect_true(exists("process_viability_data_v2", mode = "function"))
})
