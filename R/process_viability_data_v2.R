#' Process viability data with a fixed 0\% baseline and column-averaged 100\%
#' control (v2 layout)
#'
#' @description
#' Third viability processor for \pkg{dosefitr}. It ports the control style of
#' \code{\link{ratio_dose_response_v2}} (the NanoBRET processor) into the
#' cell-viability world:
#'
#' \itemize{
#'   \item The \strong{0\% control} is a \strong{fixed scalar} supplied by the
#'         user (\code{control_0perc}). It becomes a constant first row
#'         (\code{Fixed_0perc}) of \code{modified_table}. It is \emph{never} a
#'         plate well, and \code{control_mean_scope} does not affect it.
#'   \item The \strong{100\% control} is one or more \strong{plate columns}
#'         (\code{control_100perc}). Their values are averaged into a single
#'         \code{Mean_100perc} row (the \strong{last row} of
#'         \code{modified_table}) and the columns are then \strong{removed} from
#'         the experimental wells (consumed as controls), exactly like
#'         \code{\link{ratio_dose_response_v2}} drops its 100\% columns.
#' }
#'
#' Like \code{\link{process_viability_data}} (v1), the plate uses a layout in
#' which each construct spans a set of plate rows and, when
#' \code{split_replicates = TRUE}, the middle (experimental) rows are halved
#' into two technical replicates that share the bracketing control rows.
#'
#' @details
#' \strong{How the 100\% mean is computed --- \code{control_mean_scope}.}
#' The last row (\code{Mean_100perc}) is built from the columns named in
#' \code{control_100perc}. The scope argument selects one of three behaviours;
#' the 0\% first row is always the fixed scalar regardless of scope.
#'
#' Suppose \code{control_100perc = c(23, 24)} (two columns) and the plate has
#' constructs A (plate rows A-D) and B (plate rows E-H):
#'
#' \describe{
#'   \item{\code{"row"} (default for v2)}{For \emph{each plate row}, average that
#'         row's values across columns 23 and 24. After the table is transposed
#'         (wells become rows, plate rows become columns), the \code{Mean_100perc}
#'         row therefore carries a \emph{different} value for every construct
#'         column --- each equal to that plate row's own mean of the 100\%
#'         columns. This is the direct analogue of the NanoBRET
#'         \code{rowMeans(ratio[, control_100_cols])} behaviour.}
#'   \item{\code{"construct"}}{For \emph{each construct}, average columns 23 and
#'         24 over \emph{all plate rows belonging to that construct}. Every column
#'         of construct A shares one value (its mean over rows A-D); every column
#'         of construct B shares another (its mean over rows E-H). This mirrors
#'         the per-construct control averaging used by v1/v2, applied to the
#'         100\% control only.}
#'   \item{\code{"global"}}{Average \emph{all cells} in columns 23 and 24 across
#'         the whole plate into a single number. The \code{Mean_100perc} row is
#'         then identical for every construct column.}
#' }
#'
#' \strong{Concentration axis.} \code{info_table[[1]]} must provide one
#' log-concentration per \emph{experimental} well (i.e. after the 100\% control
#' columns are removed), per replicate. The leading column of the output is
#' assembled as \code{c(NA, exp_concs, NA)}, with \code{NA} at the
#' \code{Fixed_0perc} (row 1) and \code{Mean_100perc} (last row) positions.
#'
#' \strong{Output row count.} With \code{k = length(selected_columns)} and
#' \code{m = } number of distinct 100\% control columns, the retained
#' experimental wells number \code{k - m}. When \code{split_replicates = TRUE},
#' those are halved (an odd count drops the last well, as in v1), giving
#' \code{e = floor((k - m) / 2)} rows per replicate. The final table has
#' \code{e + 2} rows: \code{Fixed_0perc}, \code{e} concentrations,
#' \code{Mean_100perc}. For the canonical 24-column plate with a single 100\%
#' column this is \code{floor(23/2) + 2 = 13} rows, matching the other viability
#' processors. With \code{split_replicates = FALSE} the table has
#' \code{(k - m) + 2} rows.
#'
#' @param data Raw plate data.frame / matrix, including any header rows/columns.
#'   Auto-detection locates the data block (rows A-H + cols 1-12 for a 96-well
#'   plate, or rows A-P + cols 1-24 for a 384-well plate).
#'
#' @param control_0perc \strong{Required.} A single finite numeric value used as
#'   the fixed 0\% baseline. It is placed as the constant first row
#'   (\code{Fixed_0perc}) of \code{modified_table}. Unlike v1/v2 this is
#'   \emph{not} a plate column: pass a number such as \code{0} (a 0\% floor) or a
#'   measured background signal. Non-numeric, length \eqn{\neq} 1, or \code{NA}
#'   values raise an error.
#'
#' @param control_100perc \strong{Required.} The 100\% (positive / untreated)
#'   control column(s). Accepts \strong{either}:
#'   \itemize{
#'     \item numeric well index/indices in \code{[1, n_cols]}, e.g. \code{24} or
#'           \code{c(23, 24)}; \strong{or}
#'     \item character column name(s), e.g. \code{"24"} or \code{c("23", "24")}.
#'   }
#'   The mean of these columns becomes the last row (\code{Mean_100perc}); see
#'   \code{control_mean_scope} for how the mean is computed. Every 100\% column
#'   must also be present in \code{selected_columns}.
#'
#' @param split_replicates Logical. If \code{TRUE} (default), the middle
#'   (experimental) rows of the transposed table are split into two technical
#'   replicates that share the \code{Fixed_0perc} / \code{Mean_100perc}
#'   bracketing rows (v1-style). The second replicate's columns receive a
#'   \code{.2} suffix. If the number of experimental wells is odd, the last well
#'   is dropped (a one-line warning is emitted).
#'
#' @param info_table \strong{Required.} A data.frame with at least four columns,
#'   in this order: \code{log(inhibitor)}, \code{Plate_Row} (A-H or A-P),
#'   \code{Construct}, \code{Compound}. Repeated \code{Construct}+\code{Compound}
#'   pairs are treated as biological replicates and disambiguated with
#'   \code{_2}, \code{_3}, ... suffixes on the construct name (as in v1).
#'
#' @param selected_columns Optional integer vector of well/column indices in
#'   \code{[1, n_cols]} to keep. Default \code{NULL} means all columns. Must
#'   include every \code{control_100perc} column and must not contain
#'   duplicates.
#'
#' @param low_value_threshold Numeric. Values strictly below this threshold are
#'   replaced with \code{NA} before any averaging. Default \code{0}.
#'
#' @param verbose Logical. If \code{TRUE} (default), print progress messages.
#'
#' @param control_mean_scope Character; one of \code{"row"} (default for v2),
#'   \code{"construct"}, or \code{"global"}. Controls how the 100\% columns are
#'   averaged into the \code{Mean_100perc} last row. See \strong{Details} for
#'   worked definitions of each mode. Does not affect the fixed 0\% row.
#'
#' @param auto_detect Logical. If \code{TRUE} (default), auto-detect the plate
#'   block. If \code{FALSE}, fall back to fixed 96-well positions analogous to
#'   v2's fallback.
#'
#' @return A list with the same shape as \code{\link{process_viability_data_v3}}:
#' \describe{
#'   \item{original_table}{Numeric data.frame of extracted viability values
#'         (rownames = A-H or A-P; colnames = \code{"1"}..\code{"n"}).}
#'   \item{modified_table}{Transposed, control-bracketed, replicate-split table
#'         with a leading \code{log(inhibitor)} column; first row
#'         \code{Fixed_0perc}, last row \code{Mean_100perc}.}
#'   \item{processing_info}{Intermediate objects used by
#'         \code{\link{batch_viability_analysis}} for QC computation, including
#'         \code{control_0_info} (with the fixed value) and \code{control_100_info}
#'         (with the column names and scope).}
#'   \item{selected_columns_info}{Description of the column selection.}
#'   \item{version}{\code{"v2"}.}
#'   \item{data_type}{\code{"viability"}.}
#'   \item{auto_detect}{Echo of the input.}
#'   \item{behavior_mode}{\code{"v2-fixed0-colmean100"}.}
#' }
#'
#' @seealso \code{\link{process_viability_data}} (v1),
#'   \code{\link{process_viability_data_v3}} (v3),
#'   \code{\link{ratio_dose_response_v2}} (NanoBRET analogue of the control
#'   style), and \code{\link{batch_viability_analysis}} (use
#'   \code{version = "v2"} to dispatch here).
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' # v2 uses a fixed 0% scalar and averages user-chosen columns for the 100%
#' # control. The bundled fixtures are not required; build a small synthetic
#' # 96-well plate (rows A-H, cols 1-12) here.
#' set.seed(1L)
#' n_hdr <- 12L
#' header_rows <- as.data.frame(
#'   matrix("", nrow = n_hdr, ncol = 13L), stringsAsFactors = FALSE
#' )
#' col_label_row <- as.data.frame(
#'   matrix(c("", as.character(1:12)), nrow = 1L), stringsAsFactors = FALSE
#' )
#' data_rows <- as.data.frame(
#'   matrix("", nrow = 8L, ncol = 13L), stringsAsFactors = FALSE
#' )
#' logc <- seq(-9, -5, length.out = 10L)   # 10 experimental wells (cols 2-11)
#' for (i in seq_len(8L)) {
#'   frac <- 1 / (1 + 10^((logc - (-7))))
#'   vals <- 40 + (800 - 40) * frac + rnorm(10L, 0, 5)
#'   data_rows[i, 1L]     <- LETTERS[i]
#'   data_rows[i, 2L:11L] <- as.character(round(vals))
#'   data_rows[i, 12L]    <- as.character(round(800 + rnorm(1L, 0, 3)))  # 100% col
#'   data_rows[i, 13L]    <- as.character(round(40  + rnorm(1L, 0, 3)))
#' }
#' raw_plate <- rbind(header_rows, col_label_row, data_rows)
#' colnames(raw_plate) <- paste0("V", seq_len(ncol(raw_plate)))
#'
#' info_table <- data.frame(
#'   log_conc  = c(NA, logc, NA),                 # concs for the 10 exp. wells
#'   Plate_Row = LETTERS[1:8],
#'   Target    = rep("KinaseZ", 8L),
#'   Compound  = paste0("Cpd", 1:8),
#'   stringsAsFactors = FALSE
#' )
#'
#' out <- process_viability_data_v2(
#'   data               = raw_plate,
#'   control_0perc      = 0,          # fixed 0% baseline (first row)
#'   control_100perc    = 12,         # column 12 averaged into last row
#'   info_table         = info_table,
#'   selected_columns   = 1:12,
#'   control_mean_scope = "row",
#'   verbose            = FALSE
#' )
#' dim(out$modified_table)
#' out$modified_table[1, 1:3]                       # Fixed_0perc row
#' out$modified_table[nrow(out$modified_table), 1:3] # Mean_100perc row
#' }
#' @importFrom stats sd
#' @export

process_viability_data_v2 <- function(data,
                                      control_0perc       = NULL,
                                      control_100perc     = NULL,
                                      split_replicates    = TRUE,
                                      info_table          = NULL,
                                      selected_columns    = NULL,
                                      low_value_threshold = 0,
                                      verbose             = TRUE,
                                      control_mean_scope  = c("row", "construct", "global"),
                                      auto_detect         = TRUE) {

  # -- 1. Argument sanity -----------------------------------------------------
  control_mean_scope <- match.arg(control_mean_scope)

  if (is.null(info_table))
    stop("v2 requires an info_table.")
  if (ncol(info_table) < 4L)
    stop("info_table must have at least 4 columns: log(inhibitor), Plate_Row, Construct, Compound.")

  # 0% control: a single finite numeric scalar (never a column).
  if (is.null(control_0perc) ||
      !is.numeric(control_0perc) || length(control_0perc) != 1L ||
      !is.finite(control_0perc))
    stop("control_0perc must be a single numeric value (the fixed 0% baseline).")
  fixed_0_value <- as.numeric(control_0perc)

  if (is.null(control_100perc))
    stop("v2 requires control_100perc (one or more 100% control columns).")

  # -- 2. Plate-block auto-detection (96- or 384-well) ------------------------
  # Parameterised block detector: scans for either an 8-letter (A-H) or
  # 16-letter (A-P) label column with a matching numeric header (1..12 / 1..24).
  # (Same strategy as process_viability_data_v3's detect_block.)
  detect_block <- function(data) {
    if (verbose) message("Looking for 96-well (A-H, 1-12) or 384-well (A-P, 1-24) block...")

    if (!is.data.frame(data))
      data <- as.data.frame(data, stringsAsFactors = FALSE)

    # Try the larger plate first so a genuine 384-well header (1..24) is not
    # mis-detected as a 96-well plate off its leading 1..12 run. The boundary
    # check below (`run_is_truncated`) is the primary guard; ordering is a
    # belt-and-braces second line of defence.
    plate_specs <- list(
      list(n_rows = 16L, n_cols = 24L, row_letters = LETTERS[1:16]),
      list(n_rows = 8L,  n_cols = 12L, row_letters = LETTERS[1:8])
    )

    # A numeric header run of length n_cols is only a real plate header if the
    # column immediately after it does NOT continue the 1..n sequence (which
    # would mean we clipped a wider plate). `vals` are the numeric tokens of the
    # candidate run; `next_tok` is the token right after the run (or NA if none).
    run_is_truncated <- function(vals, next_tok, n_cols) {
      if (length(next_tok) == 0L || is.na(next_tok)) return(FALSE)
      nx <- suppressWarnings(as.numeric(next_tok))
      !is.na(nx) && nx == (n_cols + 1L)
    }

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
          next_tok <- if (start_col + n_cols <= length(current_colnames))
            current_colnames[start_col + n_cols] else NA
          if (length(valid) >= n_cols - 4L &&
              all(valid >= 1 & valid <= n_cols) &&
              !run_is_truncated(valid, next_tok, n_cols)) {
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
            next_tok <- if (start_col + n_cols <= length(row_values))
              row_values[start_col + n_cols] else NA
            if (length(valid) >= n_cols - 4L &&
                all(valid >= 1 & valid <= n_cols) &&
                !run_is_truncated(valid, next_tok, n_cols)) {
              header_row <- i
              header_start_col <- start_col
              break
            }
          }
          if (!is.null(header_row)) break
        }
      }
      if (is.null(header_row)) next

      # 2c. Look for the expected row letters in column header_start_col.
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

      if (verbose)
        message("  Detected ", n_rows, "x", n_cols, " plate (",
                if (n_rows == 8L) "96-well" else "384-well", ").")

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
    pos <- list(
      n_rows = 8L, n_cols = 12L,
      header_row = 1L, data_start_row = 2L, data_end_row = 9L,
      label_col = 1L, first_data_col = 2L, last_data_col = 13L,
      found_by_colnames = FALSE
    )
    n_rows <- pos$n_rows; n_cols <- pos$n_cols
    if (verbose) message("Auto-detection disabled; assuming 96-well fixed positions.")
  }

  if (!(n_rows == 8L && n_cols == 12L) && !(n_rows == 16L && n_cols == 24L))
    stop(sprintf("v2 expects a 96-well (8x12) or 384-well (16x24) plate; detected %dx%d.",
                 n_rows, n_cols))

  # -- 4. Extract the numeric plate block -------------------------------------
  block <- data[pos$data_start_row:pos$data_end_row, pos$label_col:pos$last_data_col]
  row_letters <- as.character(unlist(block[, 1L]))
  via <- as.data.frame(
    lapply(block[, -1L, drop = FALSE],
           function(x) suppressWarnings(as.numeric(as.character(x)))),
    stringsAsFactors = FALSE
  )
  rownames(via) <- row_letters
  colnames(via) <- as.character(seq_len(n_cols))

  # -- 5. Map control_100perc to column names ---------------------------------
  # Accept numeric well indices OR character column names.
  map_100_cols <- function(spec) {
    if (is.numeric(spec)) {
      idx <- as.integer(spec)
      if (any(idx < 1L) || any(idx > n_cols))
        stop(sprintf("control_100perc column index out of range; must be in [1, %d].", n_cols))
      return(as.character(idx))
    }
    if (is.character(spec)) {
      # Column names in `via` are "1".."n_cols".
      bad <- spec[!spec %in% colnames(via)]
      if (length(bad) > 0L)
        stop("control_100perc column name(s) not found in plate columns: ",
             paste(bad, collapse = ", "),
             sprintf(" (valid names are \"1\"..\"%d\").", n_cols))
      return(as.character(spec))
    }
    stop("control_100perc must be numeric well indices or character column names.")
  }
  ctrl100_names <- unique(map_100_cols(control_100perc))
  m <- length(ctrl100_names)

  # -- 6. Validate selected_columns -------------------------------------------
  if (is.null(selected_columns))
    selected_columns <- seq_len(n_cols)
  if (!is.numeric(selected_columns))
    stop("selected_columns must be a numeric vector of well indices.")
  if (any(selected_columns < 1) || any(selected_columns > n_cols))
    stop(sprintf("selected_columns must be in [1, %d].", n_cols))
  selected_columns <- as.integer(selected_columns)
  if (anyDuplicated(selected_columns) > 0L)
    stop("selected_columns must not contain duplicates.")
  missing_ctrl <- ctrl100_names[!ctrl100_names %in% as.character(selected_columns)]
  if (length(missing_ctrl) > 0L)
    stop("All control_100perc columns must be included in selected_columns; missing: ",
         paste(missing_ctrl, collapse = ", "), ".")

  # Keep only the requested wells, preserving input order.
  via_selected <- via[, as.character(selected_columns), drop = FALSE]

  # -- 7. Low-value filtering -------------------------------------------------
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

  # -- 8. Build replicate-aware info_table (biological replicates) ------------
  base_id_values <- paste(info_table[[3L]], info_table[[4L]], sep = ":")
  info_table$Base_ID <- base_id_values

  id_counts     <- table(base_id_values)
  duplicate_ids <- names(id_counts)[id_counts > 1]

  if (length(duplicate_ids) > 0) {
    suffix_counter       <- stats::setNames(rep(1, length(duplicate_ids)), duplicate_ids)
    new_construct_values <- info_table[[3L]]
    for (i in seq_along(base_id_values)) {
      cur <- base_id_values[i]
      if (cur %in% duplicate_ids) {
        if (suffix_counter[cur] > 1)
          new_construct_values[i] <- paste0(info_table[[3L]][i], "_", suffix_counter[cur])
        suffix_counter[cur] <- suffix_counter[cur] + 1
      }
    }
    info_table$Construct_Modified <- new_construct_values
    if (verbose)
      message("Found and automatically distinguished ", length(duplicate_ids),
              " biological replicate(s): ", paste(duplicate_ids, collapse = ", "))
  } else {
    info_table$Construct_Modified <- info_table[[3L]]
  }
  info_table$ID <- paste(info_table$Construct_Modified, info_table[[4L]], sep = ":")

  # Validate plate-row letters present in info_table are within the plate.
  expected_letters <- LETTERS[seq_len(n_rows)]
  info_plate_rows  <- info_table[[2L]]
  info_plate_rows_nonNA <- info_plate_rows[!is.na(info_plate_rows)]
  bad_rows <- info_plate_rows_nonNA[!(info_plate_rows_nonNA %in% expected_letters)]
  if (length(bad_rows) > 0L)
    stop(sprintf("info_table contains Plate_Row entries not in A-%s: %s",
                 expected_letters[length(expected_letters)],
                 paste(unique(bad_rows), collapse = ", ")))
  if (anyDuplicated(info_plate_rows_nonNA) > 0L)
    stop("info_table contains duplicated Plate_Row letters; each plate row must appear at most once.")

  # -- 9. Compute the Mean_100perc value(s) per scope -------------------------
  # Result: a named numeric vector `mean_100_by_letter`, one value per plate-row
  # letter present in `via_selected`. After transpose this becomes the last row.
  ctrl_block <- via_selected[, ctrl100_names, drop = FALSE]
  letters_present <- rownames(via_selected)

  mean_100_by_letter <- stats::setNames(rep(NA_real_, length(letters_present)),
                                        letters_present)

  if (control_mean_scope == "global") {
    g <- mean(as.matrix(ctrl_block), na.rm = TRUE)
    if (is.nan(g))
      stop("All 100% control cells are NA; cannot compute a global mean.")
    mean_100_by_letter[] <- g
    if (verbose)
      message(sprintf("control_mean_scope = 'global': Mean_100perc = %.4g (all columns).", g))

  } else if (control_mean_scope == "row") {
    # Per plate-row mean across the 100% columns.
    rm_vals <- rowMeans(ctrl_block, na.rm = TRUE)
    rm_vals[is.nan(rm_vals)] <- NA_real_
    mean_100_by_letter[names(rm_vals)] <- rm_vals
    if (verbose)
      message("control_mean_scope = 'row': each plate row averaged across its 100% column(s).")

  } else { # "construct"
    # Per-construct mean across the 100% columns over that construct's rows.
    plate_row_values <- info_table[[2L]]
    construct_values <- info_table$Construct_Modified
    keep <- !is.na(plate_row_values) & (plate_row_values %in% letters_present)
    plate_row_values <- plate_row_values[keep]
    construct_values <- construct_values[keep]

    for (cn in unique(construct_values)) {
      cn_letters <- plate_row_values[construct_values == cn]
      cn_letters <- cn_letters[cn_letters %in% letters_present]
      if (length(cn_letters) == 0L) next
      cval <- mean(as.matrix(ctrl_block[cn_letters, , drop = FALSE]), na.rm = TRUE)
      if (is.nan(cval)) {
        cval <- NA_real_
        if (verbose)
          warning(sprintf("Construct '%s' has all-NA 100%% control cells; Mean_100perc set to NA.", cn))
      }
      mean_100_by_letter[cn_letters] <- cval
      if (verbose)
        message(sprintf("  Construct '%s' (rows %s): Mean_100perc = %s",
                        cn, paste(cn_letters, collapse = ","),
                        if (is.na(cval)) "NA" else sprintf("%.4g", cval)))
    }
  }

  # -- 10. Build the experimental block (drop the 100% columns) ---------------
  exp_names <- setdiff(colnames(via_selected), ctrl100_names)
  if (length(exp_names) == 0L)
    stop("After removing the 100% control column(s), no experimental wells remain.")
  exp_block <- via_selected[, exp_names, drop = FALSE]

  # -- 11. Assemble the pre-transpose table: Fixed_0perc | exp | Mean_100perc -
  # Columns are wells; rows are plate-row letters. We add two synthetic well
  # columns (the controls) and place them first/last, matching NanoBRET's
  # Fixed_0perc / Mean_100perc column construction before transpose.
  assembled <- exp_block
  assembled$Fixed_0perc  <- fixed_0_value
  assembled$Mean_100perc <- mean_100_by_letter[rownames(assembled)]
  assembled <- assembled[, c("Fixed_0perc", exp_names, "Mean_100perc"), drop = FALSE]

  # -- 12. Transpose: wells become rows, plate-row letters become columns -----
  via_t <- as.data.frame(t(assembled), stringsAsFactors = FALSE)
  colnames(via_t) <- rownames(assembled)   # plate-row letters

  # -- 13. Rename plate-row columns using info_table IDs ----------------------
  letter_to_id <- stats::setNames(info_table$ID, info_table[[2L]])
  new_colnames <- letter_to_id[colnames(via_t)]
  colnames(via_t) <- ifelse(is.na(new_colnames), colnames(via_t), new_colnames)

  # -- 14. Split technical replicates (v1-style) ------------------------------
  final_table <- if (split_replicates) {
    split_replicates_func <- function(df) {
      n <- nrow(df)
      if (n < 3) { warning("Not enough rows to split replicates."); return(df) }
      ctrl <- c(1, n)                 # row 1 = Fixed_0perc, row n = Mean_100perc
      exp  <- 2:(n - 1)
      sp   <- floor(length(exp) / 2)
      if (sp < 1L) { warning("Not enough experimental rows to split replicates."); return(df) }
      r1   <- exp[1:sp]
      r2   <- exp[(sp + 1):length(exp)]
      if (length(r1) != length(r2)) {
        ml <- min(length(r1), length(r2))
        warning(sprintf(
          "Odd number of experimental wells (%d); dropping the last well to equalise replicates.",
          length(exp)))
        r1 <- r1[1:ml]; r2 <- r2[1:ml]
      }
      out <- data.frame()
      for (col in colnames(df)) {
        v1 <- df[c(ctrl[1], r1, ctrl[2]), col]
        v2 <- df[c(ctrl[1], r2, ctrl[2]), col]
        if (ncol(out) == 0) {
          out <- data.frame(v1, v2)
          colnames(out) <- c(col, paste0(col, ".2"))
        } else {
          out[[col]]               <- v1
          out[[paste0(col, ".2")]] <- v2
        }
      }
      rownames(out) <- c(rownames(df)[ctrl[1]], rownames(df)[r1], rownames(df)[ctrl[2]])
      out
    }
    split_replicates_func(via_t)
  } else {
    via_t
  }

  # -- 15. Build the leading log(inhibitor) column ----------------------------
  # Concentrations are for the experimental wells only, per replicate.
  raw_log  <- info_table[[1L]]
  conc_vals <- suppressWarnings(as.numeric(raw_log))
  conc_vals <- conc_vals[!is.na(conc_vals)]

  n_needed <- nrow(final_table)
  n_exp    <- n_needed - 2L   # rows excluding Fixed_0perc (1) and Mean_100perc (last)

  if (length(conc_vals) < n_exp)
    stop(sprintf(paste0("info_table[[1]] provides %d numeric concentration(s), but the ",
                        "experimental block needs %d (one per experimental well per ",
                        "replicate). Provide at least %d concentrations."),
                 length(conc_vals), n_exp, n_exp))
  if (length(conc_vals) > n_exp && verbose)
    message(sprintf("info_table[[1]] has %d concentrations; using the first %d for the experimental wells.",
                    length(conc_vals), n_exp))

  exp_concs <- conc_vals[seq_len(n_exp)]
  log_col   <- c(NA_real_, exp_concs, NA_real_)

  # Preserve the (post-split) row names -- Fixed_0perc / conc rows / Mean_100perc
  # -- across the cbind, which can otherwise reset them.
  final_rownames_out <- rownames(final_table)
  final_table <- cbind(log_col, final_table)
  colnames(final_table)[1L] <- colnames(info_table)[1L]
  rownames(final_table) <- final_rownames_out

  # -- 16. Assemble processing_info (consumed by batch QC) --------------------
  control_0_info <- list(
    name        = "Fixed_0perc",
    user_index  = NA_integer_,
    is_relative = FALSE,
    fixed_value = fixed_0_value
  )
  control_100_info <- list(
    name         = ctrl100_names,     # character vector of the 100% column names
    user_index   = control_100perc,
    is_relative  = TRUE,
    scope        = control_mean_scope
  )

  result <- list(
    original_table = via_selected,
    modified_table = final_table,
    processing_info = list(
      viability_data      = via_selected,
      viability_assembled = assembled,
      control_0_info      = control_0_info,
      control_100_info    = control_100_info,
      info_table          = info_table,
      final_rownames      = row_letters,
      selected_columns    = selected_columns,
      control_mean_scope  = control_mean_scope,
      fixed_0_value       = fixed_0_value,
      auto_detected       = auto_detect,
      detection_method    = if (auto_detect) {
        if (isTRUE(pos$found_by_colnames)) "column_names" else "row_header"
      } else "fixed_positions",
      table_positions     = pos
    ),
    selected_columns_info = list(
      user_indices     = selected_columns,
      description      = paste0("Selected data columns: ",
                                paste(selected_columns, collapse = ", "),
                                " (well indices in [1, ", n_cols, "]); ",
                                "100% control columns consumed: ",
                                paste(ctrl100_names, collapse = ", ")),
      all_data_columns = seq_len(n_cols),
      behavior         = "v2 (fixed 0% scalar first row; column-averaged 100% last row)"
    ),
    version       = "v2",
    data_type     = "viability",
    auto_detect   = auto_detect,
    behavior_mode = "v2-fixed0-colmean100"
  )

  return(result)
}
