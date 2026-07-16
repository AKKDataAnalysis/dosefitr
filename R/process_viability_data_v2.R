#' Process viability data from a plate with row-based replicates (v2 layout)
#'
#' Variant of \code{\link{process_viability_data}} for a different plate-layout
#' convention. In v2:
#' \itemize{
#'   \item Replicates are encoded as repeated rows in \code{info_table} that
#'         share the same \code{Target}/\code{Compound} pair. The first
#'         occurrence keeps the name \code{Target:Compound}; subsequent rows
#'         get suffixes \code{.2}, \code{.3}, \code{.4}, ... appended to the
#'         column name in the final output (the construct name itself is not
#'         modified).
#'   \item Concentrations live in the first column of \code{info_table} and
#'         map to wells (columns of the plate), not to plate rows. Only the
#'         non-NA numeric entries are used; the function builds the final
#'         concentration column as \code{c(NA, concs, NA)}, with NA at the
#'         positions corresponding to \code{control_0perc} (first row of the
#'         final table) and \code{control_100perc} (last row).
#'   \item Optimised for a 96-well plate (rows A-H, columns 1-12) by default,
#'         but also accepts a 384-well plate (rows A-P, columns 1-24) if the
#'         auto-detected data block matches that shape.
#' }
#'
#' @param data           Raw plate data.frame / matrix, including any header
#'                       rows/columns. Auto-detection will locate the data
#'                       block.
#' @param control_0perc  Single integer in \code{[1, n_cols]} -- well used as
#'                       the 0\% (background) control. Must be one of
#'                       \code{selected_columns}.
#' @param control_100perc Single integer in \code{[1, n_cols]} -- well used as
#'                       the 100\% (positive) control. Must be one of
#'                       \code{selected_columns}.
#' @param split_replicates Accepted for API parity with v1. v2 never splits a
#'                       row (replicates are defined by repeated rows in
#'                       \code{info_table}), so if \code{TRUE} the function
#'                       forces it to \code{FALSE} and emits a one-line
#'                       message.
#' @param info_table     data.frame with at least four columns, in this order:
#'                       \code{log(inhibitor)}, \code{Plate_Row}, \code{Target},
#'                       \code{Compound}. The first column may contain NAs
#'                       (they are dropped); \code{Plate_Row} may have trailing
#'                       NA rows (they exist only to extend the concentration
#'                       vector and are ignored).
#' @param selected_columns Integer vector of well/column indices in
#'                       \code{[1, n_cols]} to keep. Both controls must be
#'                       included. Default \code{NULL} means \code{1:n_cols}.
#' @param low_value_threshold Numeric. Values strictly below this threshold
#'                       are replaced with NA. Default \code{0}.
#' @param verbose        Logical. If \code{TRUE}, print progress messages.
#' @param apply_control_means Logical. If \code{TRUE}, replace each row's
#'                       \code{control_0perc} and \code{control_100perc}
#'                       wells with the per-Target mean across all plate rows
#'                       that share the same Target (across all compounds for
#'                       that Target).
#' @param auto_detect    Logical. If \code{TRUE} (default), auto-detect the
#'                       plate block. If \code{FALSE}, attempt the same fixed
#'                       fallback positions as v1.
#'
#' @return A list with the same shape as v1's \code{\link{process_viability_data}}:
#' \describe{
#'   \item{original_table}{Numeric data.frame of extracted viability values
#'         (rownames = A-H or A-P; colnames = "1".."12" or "1".."24").}
#'   \item{modified_table}{Transposed, control-reordered, replicate-renamed
#'         table with a leading log(inhibitor) column.}
#'   \item{processing_info}{Intermediate objects used by
#'         \code{\link{batch_viability_analysis}} for QC computation.}
#'   \item{selected_columns_info}{Description of column selection.}
#'   \item{version}{\code{"v2"}.}
#'   \item{data_type}{\code{"viability"}.}
#'   \item{auto_detect}{Echo of the input.}
#'   \item{behavior_mode}{\code{"v2-row-replicates"}.}
#' }
#'
#' @seealso \code{\link{process_viability_data}} (v1, two-replicates-per-row
#' layout) and \code{\link{batch_viability_analysis}} (use
#' \code{version = "v2"} to dispatch to this function).
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' # v2 targets a 96-well plate (rows A-H, cols 1-12) with row-based replicates.
#' # The bundled fixtures are 384-well, so we build a small synthetic plate here.
#' n_hdr <- 12L
#' header_rows <- as.data.frame(
#'   matrix("", nrow = n_hdr, ncol = 13L), stringsAsFactors = FALSE
#' )
#' col_label_row <- as.data.frame(
#'   matrix(c("", as.character(1:12)), nrow = 1L), stringsAsFactors = FALSE
#' )
#' set.seed(1L)
#' data_rows <- as.data.frame(
#'   matrix("", nrow = 8L, ncol = 13L), stringsAsFactors = FALSE
#' )
#' logc <- seq(-9, -5, length.out = 10L)
#' for (i in seq_len(8L)) {
#'   data_rows[i, 1L] <- LETTERS[i]
#'   frac <- 1 / (1 + 10^((logc - (-7))))
#'   vals <- 40 + (800 - 40) * frac + rnorm(10L, 0, 5)
#'   data_rows[i, 2L]     <- as.character(round(800 + rnorm(1L, 0, 3)))
#'   data_rows[i, 3L:12L] <- as.character(round(vals))
#'   data_rows[i, 13L]    <- as.character(round(40 + rnorm(1L, 0, 3)))
#' }
#' raw_plate <- rbind(header_rows, col_label_row, data_rows)
#' colnames(raw_plate) <- paste0("V", seq_len(ncol(raw_plate)))
#'
#' info_table <- data.frame(
#'   log_conc  = c(NA, logc, NA),
#'   Plate_Row = c(LETTERS[1:8], NA, NA, NA, NA),
#'   Target    = c(rep("KinaseZ", 8L), NA, NA, NA, NA),
#'   Compound  = c(paste0("Cpd", 1:8), NA, NA, NA, NA),
#'   stringsAsFactors = FALSE
#' )
#'
#' out <- process_viability_data_v2(
#'   data             = raw_plate,
#'   control_0perc    = 1,
#'   control_100perc  = 12,
#'   info_table       = info_table,
#'   selected_columns = 1:12,
#'   verbose          = FALSE
#' )
#' dim(out$modified_table)
#' }
#' @importFrom stats ave
#' @export

process_viability_data_v2 <- function(data,
                                      control_0perc       = NULL,
                                      control_100perc     = NULL,
                                      split_replicates    = TRUE,
                                      info_table          = NULL,
                                      selected_columns    = NULL,
                                      low_value_threshold = 0,
                                      verbose             = TRUE,
                                      apply_control_means = TRUE,
                                      auto_detect         = TRUE) {

  # -- 1. Argument sanity -----------------------------------------------------
  if (is.null(info_table))
    stop("v2 requires an info_table.")
  if (ncol(info_table) < 4L)
    stop("info_table must have at least 4 columns: log(inhibitor), Plate_Row, Target, Compound.")

  if (isTRUE(split_replicates)) {
    if (verbose)
      message("split_replicates is ignored in v2; replicates are defined by repeated rows in info_table.")
    split_replicates <- FALSE
  }

  # -- 2. Plate-block auto-detection (96- or 384-well) ------------------------
  # Parameterised version of v1's detect_specific_pattern: scans for either an
  # 8-letter (A-H) or 16-letter (A-P) label column, with a matching numeric
  # header (1..12 or 1..24).
  detect_block <- function(data) {
    if (verbose) message("Looking for 96-well (A-H, 1-12) or 384-well (A-P, 1-24) block...")

    if (!is.data.frame(data))
      data <- as.data.frame(data, stringsAsFactors = FALSE)

    plate_specs <- list(
      list(n_rows = 8L,  n_cols = 12L, row_letters = LETTERS[1:8]),
      list(n_rows = 16L, n_cols = 24L, row_letters = LETTERS[1:16])
    )

    for (spec in plate_specs) {
      n_rows <- spec$n_rows; n_cols <- spec$n_cols
      expected_letters <- spec$row_letters

      header_row <- NULL; header_start_col <- NULL; found_by_colnames <- FALSE

      # 2a. Look in current column names first.
      current_colnames <- colnames(data)
      if (length(current_colnames) >= n_cols) {
        for (start_col in seq_len(length(current_colnames) - n_cols + 1L)) {
          potential <- current_colnames[start_col:(start_col + n_cols - 1L)]
          nums <- suppressWarnings(as.numeric(potential))
          valid <- nums[!is.na(nums)]
          if (length(valid) >= n_cols - 4L &&
              all(valid >= 1 & valid <= n_cols)) {
            header_row <- NA
            header_start_col <- start_col
            found_by_colnames <- TRUE
            break
          }
        }
      }

      # 2b. Otherwise scan up to 30 leading rows for a numeric 1..n_cols run.
      if (is.null(header_row)) {
        for (i in seq_len(min(30L, nrow(data)))) {
          row_values <- as.character(unlist(data[i, ], use.names = FALSE))
          if (length(row_values) < n_cols) next
          for (start_col in seq_len(length(row_values) - n_cols + 1L)) {
            potential <- row_values[start_col:(start_col + n_cols - 1L)]
            nums <- suppressWarnings(as.numeric(potential))
            valid <- nums[!is.na(nums)]
            if (length(valid) >= n_cols - 4L &&
                all(valid >= 1 & valid <= n_cols)) {
              header_row <- i
              header_start_col <- start_col
              break
            }
          }
          if (!is.null(header_row)) break
        }
      }
      if (is.null(header_row)) next

      # 2c. Look for the expected row letters in column header_start_col,
      #     starting one row after the header (or row 1 if found_by_colnames).
      search_start <- if (found_by_colnames) 1L else header_row + 1L
      data_start_row <- NULL; data_end_row <- NULL
      end_scan <- min(search_start + 50L, nrow(data))
      for (i in search_start:end_scan) {
        if (header_start_col > ncol(data)) break
        first_val <- as.character(data[i, header_start_col])
        if (first_val %in% expected_letters) {
          if (is.null(data_start_row)) data_start_row <- i
          data_end_row <- i
        } else if (!is.null(data_start_row)) {
          break
        }
      }
      if (is.null(data_start_row)) next

      found_letters <- as.character(unlist(data[data_start_row:data_end_row, header_start_col]))
      if (!all(expected_letters %in% found_letters)) next

      if (verbose) {
        message("  Detected ", n_rows, "x", n_cols, " plate (",
                if (n_rows == 8L) "96-well" else "384-well", ").")
      }

      return(list(
        n_rows = n_rows, n_cols = n_cols,
        header_row = header_row,
        data_start_row = data_start_row,
        data_end_row = data_end_row,
        label_col = header_start_col,
        first_data_col = header_start_col + 1L,
        last_data_col = header_start_col + n_cols,
        found_by_colnames = found_by_colnames
      ))
    }
    return(NULL)
  }

  # -- 3. Run detection (or fall back to fixed positions) ---------------------
  if (auto_detect) {
    if (verbose) message("\n=== v2 AUTO-DETECTION ===")
    pos <- detect_block(data)
    if (is.null(pos))
      stop("v2 expects a 96-well (8x12) or 384-well (16x24) plate; no such block was detected in the input.")
    n_rows <- pos$n_rows; n_cols <- pos$n_cols
  } else {
    # Fallback: assume 96-well at fixed positions analogous to v1's fallback.
    pos <- list(
      n_rows = 8L, n_cols = 12L,
      header_row = 1L, data_start_row = 2L, data_end_row = 9L,
      label_col = 1L, first_data_col = 2L, last_data_col = 13L,
      found_by_colnames = FALSE
    )
    n_rows <- pos$n_rows; n_cols <- pos$n_cols
    if (verbose) message("Auto-detection disabled; assuming 96-well fixed positions.")
  }

  # Hard guard on plate shape (the plan promises 8x12 or 16x24).
  if (!(n_rows == 8L && n_cols == 12L) && !(n_rows == 16L && n_cols == 24L))
    stop(sprintf("v2 expects a 96-well (8x12) or 384-well (16x24) plate; detected %dx%d.",
                 n_rows, n_cols))

  # -- 4. Extract the numeric plate block -------------------------------------
  block <- data[pos$data_start_row:pos$data_end_row, pos$label_col:pos$last_data_col]
  row_letters <- as.character(unlist(block[, 1L]))
  via <- as.data.frame(lapply(block[, -1L, drop = FALSE], function(x) suppressWarnings(as.numeric(as.character(x)))),
                       stringsAsFactors = FALSE)
  rownames(via) <- row_letters
  colnames(via) <- as.character(seq_len(n_cols))

  # -- 5. Validate controls and selected_columns ------------------------------
  if (is.null(control_0perc) || is.null(control_100perc))
    stop("v2 requires both control_0perc and control_100perc.")
  if (!(is.numeric(control_0perc) && length(control_0perc) == 1L &&
        control_0perc >= 1 && control_0perc <= n_cols))
    stop(sprintf("control_0perc must be a single integer in [1, %d].", n_cols))
  if (!(is.numeric(control_100perc) && length(control_100perc) == 1L &&
        control_100perc >= 1 && control_100perc <= n_cols))
    stop(sprintf("control_100perc must be a single integer in [1, %d].", n_cols))
  if (control_0perc == control_100perc)
    stop("control_0perc and control_100perc must be different wells.")

  if (is.null(selected_columns))
    selected_columns <- seq_len(n_cols)
  if (!is.numeric(selected_columns))
    stop("selected_columns must be a numeric vector of well indices.")
  if (any(selected_columns < 1) || any(selected_columns > n_cols))
    stop(sprintf("selected_columns must be in [1, %d].", n_cols))
  selected_columns <- as.integer(selected_columns)
  if (anyDuplicated(selected_columns) > 0)
    stop("selected_columns must not contain duplicates.")
  if (!(control_0perc   %in% selected_columns))
    stop("control_0perc must be included in selected_columns.")
  if (!(control_100perc %in% selected_columns))
    stop("control_100perc must be included in selected_columns.")

  # Keep only the requested wells, preserving input order.
  via_selected <- via[, as.character(selected_columns), drop = FALSE]

  # -- 6. Low-value filtering -------------------------------------------------
  for (col in colnames(via_selected)) {
    low <- !is.na(via_selected[[col]]) & via_selected[[col]] < low_value_threshold
    if (any(low)) {
      if (verbose)
        warning(sprintf("Replaced %d value(s) < %s with NA in column '%s'.",
                        sum(low), format(low_value_threshold), col))
      via_selected[[col]][low] <- NA
    }
  }
  if (any(via_selected == 0, na.rm = TRUE)) {
    if (verbose) warning("Zero values detected - replacing with NA.")
    via_selected[via_selected == 0] <- NA
  }

  # -- 7. Build replicate-aware info_table ------------------------------------
  info_raw <- info_table
  plate_row_col <- info_raw[[2L]]
  target_col    <- info_raw[[3L]]
  compound_col  <- info_raw[[4L]]

  keep_idx <- !is.na(plate_row_col) & !is.na(target_col) & !is.na(compound_col)
  info_kept <- info_raw[keep_idx, , drop = FALSE]

  if (nrow(info_kept) == 0L)
    stop("info_table has no rows with a non-NA Plate_Row/Target/Compound.")

  # Validate plate-row letters
  expected_letters <- LETTERS[seq_len(n_rows)]
  bad_rows <- !(info_kept[[2L]] %in% expected_letters)
  if (any(bad_rows))
    stop(sprintf("info_table contains Plate_Row entries not in A-%s: %s",
                 expected_letters[length(expected_letters)],
                 paste(unique(info_kept[[2L]][bad_rows]), collapse = ", ")))

  if (anyDuplicated(info_kept[[2L]]) > 0L)
    stop("info_table contains duplicated Plate_Row letters; each plate row must appear at most once.")

  # Suffix assignment: first occurrence keeps base name, subsequent get .2/.3/...
  base_id <- paste(info_kept[[3L]], info_kept[[4L]], sep = ":")
  occurrence <- ave(seq_along(base_id), base_id, FUN = seq_along)
  output_colname <- ifelse(occurrence == 1L,
                           base_id,
                           paste0(base_id, ".", occurrence))

  info_kept$Base_ID            <- base_id
  info_kept$Construct_Modified <- info_kept[[3L]]   # construct name NOT suffixed
  info_kept$ID                 <- output_colname    # output column name (with suffix)

  if (verbose) {
    n_reps <- table(base_id)
    multi <- names(n_reps)[n_reps > 1L]
    if (length(multi) > 0L)
      message("Found ", length(multi), " replicate group(s): ",
              paste(sprintf("%s (n=%d)", multi, n_reps[multi]), collapse = ", "))
  }

  # -- 8. Apply per-Target control means (across all compounds of that Target) -
  via_modified <- via_selected
  if (isTRUE(apply_control_means)) {
    ctrl0_name <- as.character(control_0perc)
    ctrl1_name <- as.character(control_100perc)

    # Map plate-row letter -> integer index in via_modified rownames
    row_idx_of_letter <- setNames(seq_len(nrow(via_modified)), rownames(via_modified))

    for (tg in unique(info_kept[[3L]])) {
      letters_for_tg <- info_kept[[2L]][info_kept[[3L]] == tg]
      idx <- row_idx_of_letter[letters_for_tg]
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0L) next

      m0 <- mean(via_modified[idx, ctrl0_name], na.rm = TRUE)
      m1 <- mean(via_modified[idx, ctrl1_name], na.rm = TRUE)
      via_modified[idx, ctrl0_name] <- m0
      via_modified[idx, ctrl1_name] <- m1

      if (verbose)
        message(sprintf("  Target '%s' (rows %s): ctrl_0=%.3g, ctrl_100=%.3g",
                        tg, paste(letters_for_tg, collapse = ","), m0, m1))
    }
  }

  # -- 9. Column reorder (controls at first/last positions) -------------------
  ctrl0_name <- as.character(control_0perc)
  ctrl1_name <- as.character(control_100perc)
  other_cols <- setdiff(colnames(via_modified), c(ctrl0_name, ctrl1_name))
  via_modified <- via_modified[, c(ctrl0_name, other_cols, ctrl1_name), drop = FALSE]

  # -- 10. Transpose: wells become rows, plate rows become columns -----------
  via_t <- as.data.frame(t(via_modified), stringsAsFactors = FALSE)
  colnames(via_t) <- rownames(via_modified)   # plate-row letters

  # -- 11. Rename plate-row columns using the info_table mapping --------------
  letter_to_id <- setNames(info_kept$ID, info_kept[[2L]])
  new_colnames <- letter_to_id[colnames(via_t)]
  unmapped <- colnames(via_t)[is.na(new_colnames)]
  if (length(unmapped) > 0L && verbose)
    warning("Plate rows without info_table entry kept their letter names: ",
            paste(unmapped, collapse = ", "))
  colnames(via_t) <- ifelse(is.na(new_colnames), colnames(via_t), new_colnames)

  # -- 12. Build the leading log(inhibitor) column ---------------------------
  raw_log <- info_raw[[1L]]
  concs <- suppressWarnings(as.numeric(raw_log))
  concs <- concs[!is.na(concs)]
  n_expected <- length(selected_columns)
  if (length(concs) != n_expected - 2L)
    stop(sprintf(paste0("Found %d numeric concentration(s) in info_table[[1]] ",
                        "but selected_columns has %d well(s); expected %d-2 = %d."),
                 length(concs), n_expected, n_expected, n_expected - 2L))

  log_col <- c(NA_real_, concs, NA_real_)
  if (nrow(via_t) != length(log_col))
    stop(sprintf("Internal: transposed table has %d rows but expected %d.",
                 nrow(via_t), length(log_col)))

  final_table <- cbind(log_col, via_t)
  colnames(final_table)[1L] <- colnames(info_raw)[1L]
  rownames(final_table) <- rownames(via_t)

  # -- 13. Assemble result ----------------------------------------------------
  control_0_info <- list(
    name = ctrl0_name, user_index = control_0perc, is_relative = TRUE
  )
  control_100_info <- list(
    name = ctrl1_name, user_index = control_100perc, is_relative = TRUE
  )

  result <- list(
    original_table = via_selected,
    modified_table = final_table,
    processing_info = list(
      viability_data                = via_selected,
      viability_modified_with_means = via_modified,
      control_0_info                = control_0_info,
      control_100_info              = control_100_info,
      info_table                    = info_kept,
      final_rownames                = row_letters,
      selected_columns              = selected_columns,
      apply_control_means           = apply_control_means,
      auto_detected                 = auto_detect,
      detection_method              = if (auto_detect) {
        if (isTRUE(pos$found_by_colnames)) "column_names" else "row_header"
      } else "fixed_positions",
      table_positions               = pos
    ),
    selected_columns_info = list(
      user_indices       = selected_columns,
      description        = paste0("Selected data columns: ",
                                  paste(selected_columns, collapse = ", "),
                                  " (well indices in [1, ", n_cols, "])"),
      all_data_columns   = seq_len(n_cols),
      behavior           = "v2 (well indices map to plate columns; replicates from rows)"
    ),
    version       = "v2",
    data_type     = "viability",
    auto_detect   = auto_detect,
    behavior_mode = "v2-row-replicates"
  )

  return(result)
}
