#' Dose-Response Ratio Processing and Normalization Function
#'
#' @description
#' This function processes dual-readout plate data (e.g., luciferase and a normalizer),
#' calculates ratios, applies optional normalization to controls (0% and 100%),
#' performs quality control metrics, handles replicate splitting, and optionally
#' exports results to Excel.
#'
#' It is designed for high-throughput dose-response experiments where two signals
#' (e.g., reporter and normalization control) are measured across multiple wells.
#'
#' @details
#' The function expects a raw data table structured similarly to plate-reader exports:
#'
#' - Row 9: column names
#' - Rows 10-25: first measurement (e.g., luciferase)
#' - Rows 28-43: second measurement (e.g., normalizer)
#'
#' The ratio is calculated as:
#' \deqn{ratio = (normalizer / luciferase) * 1000}
#'
#' Additional processing includes:
#' \itemize{
#'   \item Filtering low luciferase values (below `low_value_threshold`)
#'   \item Handling division by zero
#'   \item Flexible control definition (fixed value or column-based)
#'   \item Construction of normalized tables with 0% and 100% controls
#'   \item Optional splitting of technical replicates
#'   \item Integration with metadata (`info_table`) for labeling and QC metrics
#'   \item Calculation of assay quality metrics:
#'     \itemize{
#'       \item Z-score
#'       \item Assay window
#'       \item Signal quality classification
#'     }
#'   \item Optional Excel export with multiple sheets
#' }
#'
#' @param data A data.frame containing the raw plate data. Must have at least 43 rows.
#'
#' @param control_0perc Defines the 0% control. Can be:
#' \itemize{
#'   \item A single numeric value (fixed baseline)
#'   \item A column name (character)
#' }
#'
#' @param control_100perc Defines the 100% control. Can be:
#' \itemize{
#'   \item Numeric indices of columns (relative to data columns, excluding row names)
#'   \item Column name(s) (character vector)
#' }
#'
#' @param label_sep Single character used to separate the construct name from
#'   the compound name in column labels (e.g. \code{"EPHA1:KK135"}).
#'   Defaults to \code{":"}.  Must match the separator used when the data
#'   were originally processed.
#'
#' @param split_replicates Logical. If TRUE, splits experimental rows into two
#' technical replicates assuming symmetrical layout.
#'
#' @param info_table Optional data.frame containing metadata. Must contain at least:
#' \itemize{
#'   \item Column 1: log(inhibitor)
#'   \item Column 2: Plate row identifiers
#'   \item Column 3: Construct
#'   \item Column 4: Compound
#' }
#'
#' @param save_to_excel Optional file path (character). If provided, results are saved
#' as an Excel file using the `openxlsx` package.
#'
#' @param verbose Logical. If TRUE, prints progress messages and warnings.
#'
#' @param low_value_threshold Numeric. Values in the luciferase table below this
#' threshold are replaced with NA to avoid noise artifacts.
#'
#' @param selected_columns Optional numeric vector specifying which data columns
#' to use (excluding the first column with row identifiers).
#'
#' @param plate_format Character override for plate format: \code{"96"} (rows
#'   A-H) or \code{"384"} (rows A-P).  \code{NULL} (default) auto-detects from
#'   the number of column-header integers found.  Use this only when
#'   auto-detection fails, e.g. a 384-well plate read with only 12 columns
#'   selected (half-plate experiment).
#'
#' @return A list containing:
#' \describe{
#'   \item{modified_ratio_table}{Final processed table (transposed, normalized,
#'   optionally split into replicates, and optionally annotated with log concentrations)}
#'
#'   \item{original_ratio_table}{Raw ratio matrix before normalization and transformation}
#'
#'   \item{general_means}{(Optional) Mean values for control columns}
#'
#'   \item{interval_means}{(Optional) Construct-level quality metrics including
#'   Z-score, assay window, and signal classification}
#'
#'   \item{construct_intervals}{Mapping of constructs to plate row intervals}
#'
#'   \item{control_info}{Details about how controls were interpreted and applied}
#'
#'   \item{selected_columns_info}{Information about column selection used}
#' }
#'
#' @section Control Handling:
#' Two normalization strategies are supported:
#'
#' \strong{1. Fixed 0% control + column-based 100% control}
#' \itemize{
#'   \item Creates a `Fixed_0perc` column
#'   \item Computes `Mean_100perc` from selected columns
#'   \item Removes original 100% control columns
#' }
#'
#' \strong{2. Column-based 0% and 100% controls}
#' \itemize{
#'   \item Keeps both control columns in the dataset
#'   \item Reorders columns for consistency
#' }
#'
#' @section Quality Metrics:
#' When both controls and `info_table` are provided, the function computes:
#'
#' \itemize{
#'   \item Z-score:
#'   \deqn{Z = 1 - (3 * (SD_{pos} + SD_{bg}) / (Mean_{pos} - Mean_{bg}))}
#'
#'   \item Assay window:
#'   \deqn{Assay\ Window = Mean_{pos} / Mean_{bg}}
#'
#'   \item Signal classification:
#'   \itemize{
#'     \item High, medium, low, or insufficient
#'   }
#' }
#'
#' @section Important Notes:
#' \itemize{
#'   \item The function assumes a fixed plate layout structure.
#'   \item Replicate splitting assumes symmetric experimental design.
#'   \item If `info_table` is NULL, no concentration column is added, which may
#'   affect compatibility with downstream analysis functions.
#'   \item Low luciferase values are automatically filtered to prevent ratio artifacts.
#' }
#'
#' @examples
#' \dontrun{
#' result <- ratio_dose_response_v2(
#'   data = my_data,
#'   control_0perc = 0,
#'   control_100perc = c(1, 2),
#'   info_table = metadata,
#'   split_replicates = TRUE,
#'   save_to_excel = "output.xlsx"
#' )
#'
#' # Access processed table
#' result$modified_ratio_table
#' }
#'
#' @importFrom stats sd
#' @importFrom utils write.csv
#'
#' @seealso
#' Useful for downstream dose-response modeling and outlier detection pipelines.
#'
#' @export

ratio_dose_response_v2 <- function(data,
                                   control_0perc = NULL,
                                   control_100perc = NULL,
                                   split_replicates = TRUE,
                                   info_table = NULL,
                                   save_to_excel = NULL,
                                   verbose = TRUE,
                                   low_value_threshold = 3000,
                                   selected_columns = NULL,
                                   plate_format = NULL,
                                   label_sep = ":") {

  # ── Internal helper: detect plate layout by content ───────────────────────
  #
  # Scans `df` for:
  #   1. A row whose integers form a consecutive sequence starting at 1
  #      (e.g. 1,2,...,8 or 1,2,...,24)  → col_header_row.
  #      This handles partial plates of any width and rejects data rows
  #      with random small integers.  The count of those integers gives
  #      n_data_cols (e.g. 12 for 96-well, 24 for 384-well).
  #   2. The expected row-label set: A-H (96-well) or A-P (384-well),
  #      overridden by `plate_format` when supplied.
  #   3. The first contiguous block of rows whose col-1 value is in the
  #      label set  → table1_rows  (donor channel, 450-80 B)
  #   4. The next such block after table1 ends  → table2_rows  (acceptor
  #      channel, 610-LP A)
  #   5. If a third block exists, a warning is issued (unexpected structure).
  #
  # Returns a list(col_header_row, table1_rows, table2_rows, n_data_cols).
  .locate_plate_tables <- function(df, plate_format = NULL) {
    n_rows <- nrow(df)

    # Helper: check whether a row is a plate column-number header.
    # A genuine header contains integers that form a consecutive sequence
    # starting at 1 (e.g. 1,2,...,8 or 1,2,...,24).  This correctly handles
    # partial plates with any number of columns >= 2, while rejecting data
    # rows whose values happen to include small integers but do not start at 1
    # and increment by 1.
    .is_col_header <- function(row_vec) {
      vals <- suppressWarnings(as.integer(as.character(row_vec)))
      vals <- sort(vals[!is.na(vals)])
      if (length(vals) < 2L) return(FALSE)
      vals[1L] == 1L && identical(vals, seq_len(length(vals)))
    }

    col_hdr_row <- NA_integer_
    n_data_cols <- NA_integer_

    for (i in seq_len(min(n_rows, 30L))) {
      if (.is_col_header(df[i, ])) {
        col_hdr_row <- i
        vals        <- suppressWarnings(as.integer(as.character(df[i, ])))
        n_data_cols <- sum(!is.na(vals) & vals >= 1L & vals <= 24L)
        break
      }
    }

    if (is.na(col_hdr_row))
      stop(
        "Could not detect the plate column-number header (a row containing ",
        "consecutive integers starting at 1, e.g. 1-12 or 1-24) in the ",
        "first 30 rows of the data frame.\n",
        "Please check that the data is a valid plate-reader export."
      )

    # Honour plate_format override; otherwise infer from n_data_cols
    expected_labels <- if (!is.null(plate_format)) {
      if (plate_format == "96")       LETTERS[1:8]
      else if (plate_format == "384") LETTERS[1:16]
      else stop("plate_format must be '96', '384', or NULL (auto-detect).")
    } else {
      if (n_data_cols <= 12L) LETTERS[1:8] else LETTERS[1:16]
    }

    find_next_block <- function(start_row) {
      if (start_row > n_rows) return(NULL)   # guard: nothing left to scan
      block_start <- NA_integer_
      block_end   <- NA_integer_
      for (i in seq(start_row, n_rows)) {
        label <- trimws(as.character(df[i, 1L]))
        if (label %in% expected_labels) {
          if (is.na(block_start)) block_start <- i
          block_end <- i
        } else if (!is.na(block_start)) {
          break
        }
      }
      if (is.na(block_start)) return(NULL)
      seq(block_start, block_end)
    }

    table1_rows <- find_next_block(col_hdr_row + 1L)
    if (is.null(table1_rows))
      stop("Could not find the first emission table (rows labelled ",
           paste(expected_labels, collapse = "/"), ") after the column-number header.")

    table2_rows <- find_next_block(max(table1_rows) + 1L)
    if (is.null(table2_rows))
      stop("Could not find the second emission table after the first one. ",
           "Expected two emission tables (donor + acceptor channels).")

    # Check for unexpected additional blocks
    table3_rows <- find_next_block(max(table2_rows) + 1L)
    if (!is.null(table3_rows))
      warning(
        "More than two emission tables detected in the data frame. ",
        "Only the first two (donor + acceptor channels) will be used. ",
        "If this is unexpected, check that the file contains exactly one plate export."
      )

    list(col_header_row = col_hdr_row,
         table1_rows    = table1_rows,
         table2_rows    = table2_rows,
         n_data_cols    = n_data_cols)
  }

  # ── Validate plate_format ─────────────────────────────────────────────────
  if (!is.null(plate_format) && !plate_format %in% c("96", "384"))
    stop("plate_format must be '96', '384', or NULL (auto-detect).")

  # ── Detect layout ─────────────────────────────────────────────────────────
  layout         <- .locate_plate_tables(data, plate_format = plate_format)
  col_header_row <- layout$col_header_row
  table1_rows    <- layout$table1_rows
  table2_rows    <- layout$table2_rows
  n_data_cols    <- layout$n_data_cols

  effective_format <- if (!is.null(plate_format)) plate_format else
    ifelse(n_data_cols <= 12L, "96", "384")

  if (verbose)
    message(sprintf(
      "Detected layout: col_header=row %d | table1=rows %d-%d | table2=rows %d-%d | %d data columns (%s-well%s)",
      col_header_row,
      min(table1_rows), max(table1_rows),
      min(table2_rows), max(table2_rows),
      n_data_cols,
      effective_format,
      if (!is.null(plate_format)) " [format overridden]" else ""
    ))

  # ── Apply column names from header row ────────────────────────────────────
  colnames(data) <- as.character(data[col_header_row, ])

  total_cols <- n_data_cols + 1L

  # ── Extract subtables ─────────────────────────────────────────────────────
  subtable1_full <- data[table1_rows, seq_len(total_cols)]   # donor channel
  subtable2_full <- data[table2_rows, seq_len(total_cols)]   # acceptor channel

  final_rownames <- subtable1_full[, 1]

  # ── Column selection ──────────────────────────────────────────────────────
  if (!is.null(selected_columns)) {
    if (!is.numeric(selected_columns))
      stop("selected_columns must be numeric indices (e.g., c(2:23))")

    data_columns_indices <- selected_columns + 1

    if (max(data_columns_indices) > ncol(subtable1_full))
      stop("Selected column index ", max(selected_columns),
           " is out of bounds. Maximum allowed: ", ncol(subtable1_full) - 1)
    if (min(data_columns_indices) < 2)
      stop("Selected column indices must be >= 1")

    if (length(selected_columns) %% 2 != 0)
      warning("Number of selected data columns is not even (", length(selected_columns),
              " columns selected). split_replicates will drop the last row to equalise lengths.")

    columns_to_keep <- c(1, data_columns_indices)
    subtable1 <- subtable1_full[, columns_to_keep, drop = FALSE]
    subtable2 <- subtable2_full[, columns_to_keep, drop = FALSE]
  } else {
    subtable1 <- subtable1_full
    subtable2 <- subtable2_full
  }

  # ── Convert to numeric ────────────────────────────────────────────────────
  convert_to_numeric_df <- function(df, rownames_vec) {
    num_df <- as.data.frame(apply(df[, -1, drop = FALSE], 2, as.numeric))
    rownames(num_df) <- rownames_vec
    num_df
  }

  subtable1_num <- convert_to_numeric_df(subtable1, final_rownames)
  subtable2_num <- convert_to_numeric_df(subtable2, final_rownames)

  # ── Control argument processing ───────────────────────────────────────────
  map_control_input <- function(control_spec, is_for_0perc = FALSE) {
    if (is.null(control_spec)) return(NULL)

    if (is.numeric(control_spec) && length(control_spec) == 1) {
      if (is_for_0perc)
        return(list(type = "value", value = control_spec, column_index = NA, column_name = NA))

      actual_col_index <- control_spec + 1
      if (actual_col_index > ncol(data))
        stop("Control column index ", control_spec, " is out of bounds. Maximum allowed: ", ncol(data) - 1)
      col_name <- as.character(data[col_header_row, actual_col_index])
      return(list(type = "column", value = NULL, column_index = actual_col_index, column_name = col_name))

    } else if (is.character(control_spec) && length(control_spec) == 1) {
      return(list(type = "column", value = NULL,
                  column_index = which(colnames(data) == control_spec)[1],
                  column_name  = control_spec))

    } else if (is.numeric(control_spec) && length(control_spec) > 1) {
      actual_indices <- control_spec + 1
      if (max(actual_indices) > ncol(data))
        stop("Control column indices out of bounds. Maximum allowed: ", ncol(data) - 1)
      col_names <- as.character(data[col_header_row, actual_indices])
      return(list(type = "columns", value = NULL,
                  column_indices = actual_indices, column_names = col_names))

    } else {
      stop(if (is_for_0perc) "control_0perc must be a single numeric value or column name"
           else               "control_100perc must be numeric positions or a column name")
    }
  }

  control_0_info     <- map_control_input(control_0perc,   is_for_0perc = TRUE)
  control_0_is_value <- !is.null(control_0_info) && control_0_info$type == "value"
  control_0_value    <- if (control_0_is_value) control_0_info$value else NULL

  if (!is.null(control_0_info)) {
    if (control_0_is_value) {
      if (verbose) message("Using fixed value ", control_0_value, " for 0% control")
    } else {
      if (verbose) message("Using column '", control_0_info$column_name, "' for 0% control")
      if (!is.null(selected_columns) && !control_0_info$column_name %in% colnames(subtable1_num))
        stop("Control column '", control_0_info$column_name, "' not found in selected columns.")
    }
  }

  control_100_info  <- map_control_input(control_100perc, is_for_0perc = FALSE)
  control_100_names <- NULL

  if (!is.null(control_100_info)) {
    control_100_names <- if (control_100_info$type == "column")  control_100_info$column_name
    else                                     control_100_info$column_names

    missing_cols <- control_100_names[!control_100_names %in% colnames(subtable1_num)]
    if (length(missing_cols) > 0 && !is.null(selected_columns))
      stop("Control columns not found in selected columns: ", paste(missing_cols, collapse = ", "))

    if (verbose && length(control_100_names) > 1)
      message("Using column names for 100% control: ", paste(control_100_names, collapse = ", "))
  }

  # ── Low-value filtering ───────────────────────────────────────────────────
  replace_low_values <- function(df, threshold = low_value_threshold) {
    for (col in colnames(df)) {
      low_vals <- df[[col]] < threshold & !is.na(df[[col]])
      if (any(low_vals)) {
        warning("Replaced ", sum(low_vals), " value(s) < ", threshold,
                " with NA in column '", col, "'")
        df[[col]][low_vals] <- NA
      }
    }
    df
  }

  subtable1_num <- replace_low_values(subtable1_num)

  if (any(subtable1_num == 0, na.rm = TRUE)) {
    warning("Division by zero detected in ratio calculation - replacing with NA")
    subtable1_num[subtable1_num == 0] <- NA
  }

  # ── Core ratio calculation ────────────────────────────────────────────────
  ratio          <- (subtable2_num / subtable1_num) * 1000
  ratio_modified <- ratio
  result         <- list()

  # ── Info table processing ─────────────────────────────────────────────────
  if (!is.null(info_table)) {
    if (ncol(info_table) < 4)
      stop("Info table must have at least 4 columns: log(inhibitor), Plate_Row, Construct, Compound")

    base_id_values <- paste(info_table[[3]], info_table[[4]], sep = label_sep)
    info_table$Base_ID <- base_id_values

    id_counts     <- table(base_id_values)
    duplicate_ids <- names(id_counts)[id_counts > 1]

    if (length(duplicate_ids) > 0) {
      suffix_counter       <- setNames(rep(1, length(duplicate_ids)), duplicate_ids)
      new_construct_values <- info_table[[3]]

      for (i in seq_along(base_id_values)) {
        cur <- base_id_values[i]
        if (cur %in% duplicate_ids) {
          if (suffix_counter[cur] > 1)
            new_construct_values[i] <- paste0(info_table[[3]][i], "_", suffix_counter[cur])
          suffix_counter[cur] <- suffix_counter[cur] + 1
        }
      }

      info_table$Construct_Modified <- new_construct_values
      if (verbose)
        message("Found and automatically distinguished ", length(duplicate_ids),
                " biological replicate(s): ", paste(duplicate_ids, collapse = ", "))
    } else {
      info_table$Construct_Modified <- info_table[[3]]
    }

    info_table$ID <- paste(info_table$Construct_Modified, info_table[[4]], sep = label_sep)
  }

  # ── Quality control calculations ──────────────────────────────────────────
  if (!is.null(info_table) && !is.null(control_0perc) && !is.null(control_100perc)) {

    plate_row_values   <- info_table[[2]]
    construct_values   <- info_table$Construct_Modified
    plate_row_to_index <- setNames(seq_along(plate_row_values), plate_row_values)
    unique_constructs  <- unique(construct_values)

    construct_groups <- lapply(setNames(unique_constructs, unique_constructs), function(cn) {
      plate_row_values[which(construct_values == cn)]
    })

    row_intervals <- lapply(construct_groups, function(plate_rows) {
      as.numeric(plate_row_to_index[plate_rows][!is.na(plate_row_to_index[plate_rows])])
    })
    row_intervals <- row_intervals[sapply(row_intervals, length) > 0]

    existing_100_cols <- control_100_names[control_100_names %in% colnames(ratio)]

    if (length(existing_100_cols) > 0 && length(row_intervals) > 0) {

      general_mean_100 <- mean(as.matrix(ratio[, existing_100_cols, drop = FALSE]), na.rm = TRUE)
      result$general_means <- data.frame(
        Type         = "General",
        Control_Type = "100%_Control",
        Mean         = general_mean_100,
        Columns_Used = paste(existing_100_cols, collapse = ", "),
        stringsAsFactors = FALSE
      )

      get_lowest_comment <- function(luc, aw, zs) {
        lvls <- c("insufficient", "low", "medium", "high")
        first_word <- function(x) {
          if (is.null(x) || is.na(x) || identical(x, "")) return(NA_character_)
          strsplit(as.character(x), " ")[[1]][1]
        }
        scores <- c(match(first_word(luc), lvls),
                    match(first_word(aw),  lvls),
                    match(first_word(zs),  lvls))
        if (all(is.na(scores))) return("insufficient")
        lvls[min(scores, na.rm = TRUE)]
      }

      interval_means_list <- list()

      for (cn in unique_constructs) {
        if (!cn %in% names(row_intervals)) next
        valid_rows <- row_intervals[[cn]]
        valid_rows <- valid_rows[valid_rows >= 1 & valid_rows <= nrow(ratio)]
        if (length(valid_rows) == 0) next

        mean_luc <- mean(as.matrix(subtable1_num[valid_rows, ]), na.rm = TRUE)

        luc_comment <- if (is.na(mean_luc))          "insufficient luciferase signal"
        else if (mean_luc > 100000)   "high (>100000)"
        else if (mean_luc > 10000)    "medium (10000<x<100000)"
        else if (mean_luc > 1000)     "low (1000<x<10000)"
        else                          "insufficient luciferase signal"

        z_score <- NA; aw <- NA; aw_comment <- NA; zs_comment <- NA
        mean_bg <- NA; mean_pos <- NA; sd_bg <- NA; sd_pos <- NA

        if (control_0_is_value && length(existing_100_cols) > 0) {
          mean_bg  <- control_0_value
          d100     <- ratio[valid_rows, existing_100_cols, drop = FALSE]
          mean_pos <- mean(as.matrix(d100), na.rm = TRUE)
          sd_pos   <- sd(as.matrix(d100),   na.rm = TRUE)
          sd_bg    <- 0

          if (!is.na(mean_pos) && !is.na(mean_bg)) {
            z_score <- if ((mean_pos - mean_bg) != 0)
              1 - (3 * (sd_pos + sd_bg) / (mean_pos - mean_bg)) else NA

            aw <- if (mean_bg != 0) mean_pos / mean_bg else NA

            aw_comment <- if (is.na(aw))    "insufficient"
            else if (aw > 3)  "high (>3)"
            else if (aw > 2)  "medium (2<x<3)"
            else if (aw > 1.5) "low (<2)"
            else              "insufficient"

            zs_comment <- if (is.na(z_score))      "insufficient"
            else if (z_score > 0.7)  "high (>0.7)"
            else if (z_score > 0.5)  "medium (0.5<x<0.7)"
            else if (z_score > 0.25) "low (<0.5)"
            else                     "insufficient"
          } else {
            aw_comment <- "insufficient"; zs_comment <- "insufficient"
          }
        }

        interval_means_list[[cn]] <- data.frame(
          Type                      = "Construct_Interval",
          Construct                 = cn,
          Average_Background        = mean_bg,
          SD_Background             = sd_bg,
          Average_Positive_Ctrl     = mean_pos,
          SD_Positive_Ctrl          = sd_pos,
          Average_luciferase_signal = mean_luc,
          Luciferase_signal_comment = luc_comment,
          Z_Score                   = z_score,
          Assay_z_Comment           = zs_comment,
          Assay_Window              = aw,
          Assay_window_Comment      = aw_comment,
          Overall_Quality           = get_lowest_comment(luc_comment, aw_comment, zs_comment),
          Rows                      = paste0(cn, " (rows ", paste(range(valid_rows), collapse = "-"), ")"),
          Rows_Count                = length(valid_rows),
          Background_Type           = ifelse(control_0_is_value, "Fixed_Value", "Column"),
          Background_Value          = ifelse(control_0_is_value, control_0_value, NA),
          stringsAsFactors          = FALSE
        )
      }

      if (length(interval_means_list) > 0) {
        result$interval_means <- do.call(rbind, interval_means_list)
        rownames(result$interval_means) <- NULL
        result$interval_means <- result$interval_means[, c(
          "Type", "Construct",
          "Average_Positive_Ctrl", "SD_Positive_Ctrl",
          "Average_Background", "SD_Background",
          "Average_luciferase_signal", "Z_Score", "Assay_Window",
          "Luciferase_signal_comment", "Assay_window_Comment",
          "Assay_z_Comment", "Overall_Quality",
          "Rows", "Rows_Count", "Background_Type", "Background_Value"
        )]
        im_clean <- result$interval_means[, -c(1, 2)]
        rownames(im_clean) <- result$interval_means$Construct
        result$interval_means <- as.data.frame(t(im_clean))
        result$construct_intervals <- row_intervals
      }
    }
  }

  # ── Control column restructuring ──────────────────────────────────────────
  existing_100_cols <- if (!is.null(control_100_names))
    control_100_names[control_100_names %in% colnames(ratio_modified)] else character(0)

  if (control_0_is_value && length(existing_100_cols) > 0) {
    mean_100perc <- if (length(existing_100_cols) == 1)
      ratio_modified[, existing_100_cols]
    else
      rowMeans(ratio_modified[, existing_100_cols, drop = FALSE], na.rm = TRUE)

    ratio_modified$Fixed_0perc  <- control_0_value
    ratio_modified$Mean_100perc <- mean_100perc
    ratio_modified <- ratio_modified[, !colnames(ratio_modified) %in% existing_100_cols]

    other_cols     <- setdiff(colnames(ratio_modified), c("Fixed_0perc", "Mean_100perc"))
    ratio_modified <- ratio_modified[, c("Fixed_0perc", other_cols, "Mean_100perc")]

    if (verbose)
      message("Created new control structure with fixed 0% = ", control_0_value,
              " and mean 100% from ", length(existing_100_cols), " column(s)")

  } else if (!control_0_is_value && !is.null(control_0_info) && length(existing_100_cols) == 1) {
    ctrl_cols <- c(control_0_info$column_name, existing_100_cols)
    missing   <- ctrl_cols[!ctrl_cols %in% colnames(ratio_modified)]
    if (length(missing) > 0)
      stop("Control columns not found: ", paste(missing, collapse = ", "))
    other_cols     <- setdiff(colnames(ratio_modified), ctrl_cols)
    ratio_modified <- ratio_modified[, c(control_0_info$column_name, other_cols, existing_100_cols)]
  }

  # ── Transpose ─────────────────────────────────────────────────────────────
  ratio_t           <- as.data.frame(t(ratio_modified))
  colnames(ratio_t) <- rownames(ratio_modified)

  # ── Rename columns with info_table IDs ───────────────────────────────────
  if (!is.null(info_table)) {
    mapping      <- setNames(info_table$ID, info_table[[2]])
    new_colnames <- mapping[colnames(ratio_t)]
    colnames(ratio_t) <- ifelse(is.na(new_colnames), colnames(ratio_t), new_colnames)
  }

  # ── Split technical replicates ────────────────────────────────────────────
  final_table <- if (split_replicates) {
    split_replicates_func <- function(df) {
      n <- nrow(df)
      if (n < 3) { warning("Not enough rows to split replicates."); return(df) }
      ctrl <- c(1, n)
      exp  <- 2:(n - 1)
      sp   <- floor(length(exp) / 2)
      r1   <- exp[1:sp]
      r2   <- exp[(sp + 1):length(exp)]
      if (length(r1) != length(r2)) {
        ml <- min(length(r1), length(r2)); r1 <- r1[1:ml]; r2 <- r2[1:ml]
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
    split_replicates_func(ratio_t)
  } else {
    ratio_t
  }

  if (is.null(info_table) && verbose)
    warning("info_table is NULL: no log(inhibitor) column added. ",
            "Result is incompatible with rout_outliers_batch().")

  # ── Add log(inhibitor) column ─────────────────────────────────────────────────────────────────────────────
  if (!is.null(info_table)) {
    n_needed <- nrow(final_table)
    raw_log  <- info_table[[1]]

    # Extract only the non-NA concentration values from the info_table.
    # This makes alignment robust regardless of whether concentrations start
    # at row 1, 2, 3, or any other row (e.g. when the first N rows are
    # background/positive-control rows with NA in the concentration column).
    conc_vals <- raw_log[!is.na(raw_log)]

    if (split_replicates) {
      # final_table structure: row 1 = ctrl_0, rows 2..(n-1) = concentrations, row n = ctrl_100.
      # Build log_col as: NA, conc_vals (padded/truncated to n-2), NA.
      n_exp    <- n_needed - 2L
      exp_vals <- if (length(conc_vals) >= n_exp) {
        conc_vals[seq_len(n_exp)]
      } else {
        c(conc_vals, rep(NA_real_, n_exp - length(conc_vals)))
      }
      log_col <- c(NA_real_, exp_vals, NA_real_)
    } else {
      # No guaranteed ctrl structure (split_replicates=FALSE).
      # Fall back to simple positional assignment: use conc_vals from row 1,
      # padding or truncating to n_needed.
      log_col <- if (length(conc_vals) >= n_needed) {
        conc_vals[seq_len(n_needed)]
      } else {
        c(conc_vals, rep(NA_real_, n_needed - length(conc_vals)))
      }
    }

    final_table <- cbind(log_col, final_table)
    colnames(final_table)[1] <- colnames(info_table)[1]
  }

  # Assign before Excel export so the sheet is written correctly
  result$original_ratio_table <- ratio

  # ── Excel export ──────────────────────────────────────────────────────────
  if (!is.null(save_to_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE))
      stop("Package 'openxlsx' is required to save Excel files.")
    tryCatch({
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Modified_Ratio_Table")
      openxlsx::writeData(wb, "Modified_Ratio_Table",
                          cbind(RowNames = rownames(final_table), final_table), rowNames = FALSE)
      if (!is.null(result$original_ratio_table)) {
        openxlsx::addWorksheet(wb, "Original_Ratio_Table")
        openxlsx::writeData(wb, "Original_Ratio_Table",
                            cbind(RowNames = rownames(result$original_ratio_table),
                                  result$original_ratio_table), rowNames = FALSE)
      }
      if (!is.null(result$general_means)) {
        openxlsx::addWorksheet(wb, "General_Means")
        openxlsx::writeData(wb, "General_Means", result$general_means)
      }
      if (!is.null(result$interval_means)) {
        openxlsx::addWorksheet(wb, "Interval_Means")
        openxlsx::writeData(wb, "Interval_Means", result$interval_means)
      }
      openxlsx::saveWorkbook(wb, save_to_excel, overwrite = TRUE)
      if (verbose) message("Excel file saved: ", save_to_excel)
    }, error = function(e) warning("Failed to save Excel file: ", e$message))
  }

  # ── Assemble result ───────────────────────────────────────────────────────
  result$modified_ratio_table <- final_table

  result$control_info <- list(
    control_0_type             = ifelse(control_0_is_value, "Fixed_Value", "Column"),
    control_0_value            = if (control_0_is_value) control_0_value else control_0perc,
    control_100_input          = control_100perc,
    control_100_actual_columns = existing_100_cols,
    new_columns_created        = if (control_0_is_value && length(existing_100_cols) > 0)
      c("Fixed_0perc", "Mean_100perc") else NULL
  )

  result$selected_columns_info <- if (!is.null(selected_columns)) {
    list(user_indices   = selected_columns,
         actual_indices = selected_columns + 1,
         description    = paste("Data columns", paste(selected_columns, collapse = ", ")))
  } else {
    list(user_indices   = seq_len(n_data_cols),
         actual_indices = seq_len(n_data_cols) + 1L,
         description    = paste0("All data columns (1:", n_data_cols, ")"))
  }

  result
}
