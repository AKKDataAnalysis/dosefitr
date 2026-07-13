#' Process and Analyze Raw Dose-Response Ratio Data
#'
#' Processes raw experimental data from dose-response assays to calculate BRET ratios,
#' perform quality assessment, and prepare data for downstream curve fitting analysis.
#' This function handles the initial data processing pipeline from raw plate reader
#' data to analysis-ready formatted data.
#'
#' @param data Data frame containing **raw dose-response experimental data** with specific
#'   structure typically exported from plate readers. Must have at least 43 rows with
#'   column names in row 9.
#' @param control_0perc Character specifying the column name for 0\% control (background control,
#'   typically vehicle-treated samples like DMSO).
#' @param control_100perc Character specifying the column name for 100\% control (positive control,
#'   typically maximum inhibition samples).
#' @param split_replicates Logical indicating whether to split experimental replicates
#'   into separate columns (default: TRUE).
#' @param info_table Data frame with experimental metadata containing at least 4 columns:
#'   log(inhibitor), Plate_Row, Construct, and Compound information.
#' @param save_to_excel Character string specifying Excel file path for saving processed results
#'   (default: NULL, no saving).
#' @param verbose Logical indicating whether to display progress messages
#' @param low_value_threshold Numeric; values below this threshold in the donor
#'   channel (subtable1) are replaced with NA. Default is 1000. Helps filter out
#'   background noise or invalid measurements.
#' @param selected_columns Numeric vector of column indices (0-based) to include
#'   in the analysis. Useful for selecting specific columns from a multi-well
#'   plate. If NULL, all columns are used.
#' @param plate_format Character override for plate format: \code{"96"} (rows
#'   A-H) or \code{"384"} (rows A-P).  \code{NULL} (default) auto-detects from
#'   the number of column-header integers found.  Use this only when
#'   auto-detection fails, e.g. a 384-well plate read with only 12 columns
#'   selected (half-plate experiment).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{original_ratio_table}: Original calculated ratio table from raw data
#'   \item \code{modified_ratio_table}: Processed and formatted table ready for analysis
#'   \item \code{general_means}: Data frame with general control means across all rows
#'   \item \code{interval_means}: Data frame with construct-specific means and quality metrics
#'   \item \code{construct_intervals}: List mapping constructs to row intervals
#' }
#'
#' @details
#' This function performs the **initial processing of raw dose-response data** from BRET assays
#' or similar experimental formats. It transforms raw experimental measurements into analysis-ready
#' data through the following pipeline:
#'
#' \strong{Raw Data Processing Pipeline:}
#' \itemize{
#'   \item \strong{Raw Data Input}: Accepts direct output from plate readers or experimental data systems
#'   \item \strong{Data Extraction}: Separates measurement subtables (typically donor and acceptor channels)
#'   \item \strong{Quality Filtering}: Removes low-intensity signals (<1000) that may represent failed wells
#'   \item \strong{Ratio Calculation}: Computes BRET ratios as (subtable2 / subtable1) * 1000
#'   \item \strong{Control Processing}: Handles control wells separately for normalization
#'   \item **Quality Assessment**: Calculates assay performance metrics (Z-score, Assay Window)
#'   \item **Data Formatting**: Transposes and structures data for downstream analysis
#' }
#'
#' \strong{Typical Raw Data Structure:}
#' The function expects **raw experimental data** in the following format:
#' \preformatted{
#'   Rows 1-8:     Instrument headers, plate layout, experimental metadata (ignored)
#'   Row 9:        Column names (well identifiers or sample names)
#'   Rows 10-25:   First measurement subtable (e.g., donor channel or timepoint 1)
#'   Rows 26-27:   Separator or additional headers (ignored)
#'   Rows 28-43:   Second measurement subtable (e.g., acceptor channel or timepoint 2)
#' }
#'
#' \strong{Typical info_table structure:}
#' \preformatted{
#' Column 1: log(inhibitor) - Numeric values of inhibitor concentrations in log scale
#' Column 2: Plate_Row - Row identifiers matching the ratio table (e.g., "A", "B", "C", etc.)
#' Column 3: Construct - Construct protein identifiers (e.g., "BRD4", "EGFR", "KRAS")
#' Column 4: Compound - Compound identifiers (e.g., "JQ1", "Gefitinib", "ARS-1620")
#' }
#'
#' \strong{Quality Metrics for Raw Data Assessment:}
#' \itemize{
#'   \item \strong{Luciferase Signal}: Assesses raw signal intensity from experimental measurements
#'   \item \strong{Z-Score}: Evaluates assay robustness from raw control data variability
#'   \item \strong{Assay Window}: Calculates dynamic range from raw control measurements
#'   \item \strong{Overall Quality}: Determines if raw data quality supports further analysis
#' }
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' plate_file  <- file.path(extdata_dir, "nanobret_plate_01.xlsx")
#' info_table  <- openxlsx::read.xlsx(
#'   file.path(extdata_dir, "nanobret_info.xlsx"), sheet = 1
#' )
#' raw <- openxlsx::read.xlsx(plate_file, sheet = 1, colNames = FALSE)
#'
#' out <- ratio_dose_response(
#'   data             = raw,
#'   control_0perc    = "1",
#'   control_100perc  = "24",
#'   info_table       = info_table,
#'   selected_columns = 1:24,
#'   verbose          = FALSE
#' )
#' dim(out$modified_ratio_table)
#' }
#' @section Raw Data Requirements:
#' The input data should be **raw experimental measurements** with:
#' \itemize{
#'   \item \strong{Direct instrument output}: Minimal preprocessing required
#'   \item \strong{Consistent structure}: Fixed row positions for data extraction
#'   \item \strong{Proper controls}: Clearly identified control wells in column names
#'   \item \strong{Metadata alignment}: Info table matching experimental design
#' }
#'
#' @section Expected Raw Data Sources:
#' \itemize{
#'   \item Plate reader exports (Tecan, BMG Labtech, PerkinElmer)
#'   \item BRET assay raw measurements
#'   \item Luminescence or fluorescence intensity data
#'   \item High-throughput screening raw data
#' }
#'
#' @seealso
#' \code{\link{fit_drc_3pl}} for the next step in the analysis pipeline
#' \code{\link{plot_dose_response}} for visualization of processed data
#'
#' @export
#'
#' @references
#' For raw data processing in dose-response assays:
#' \itemize{
#'   \item "BRET Assay Development Guide" (PerkinElmer)
#'   \item "High-Throughput Screening Data Analysis" (Inglese et al.)
#'   \item Journal of Biomolecular Screening raw data standards
#' }



ratio_dose_response <- function(data,
                                control_0perc = NULL, control_100perc = NULL,
                                split_replicates = TRUE, info_table = NULL,
                                save_to_excel = NULL, verbose = TRUE,
                                low_value_threshold = 1000,
                                selected_columns = NULL,
                                plate_format = NULL) {
  
  # -- Internal helper: detect plate layout by content ------------------------
  #
  # Scans `df` for:
  #   1. A row whose integers form a consecutive sequence starting at 1
  #      (e.g. 1,2,...,8 or 1,2,...,24)  -> col_header_row.
  #      This handles partial plates of any width and rejects data rows
  #      with random small integers.  The count of those integers gives
  #      n_data_cols (e.g. 12 for 96-well, 24 for 384-well).
  #   2. The expected row-label set: A-H (96-well) or A-P (384-well),
  #      overridden by `plate_format` when supplied.
  #   3. The first contiguous block of rows whose col-1 value is in the
  #      label set  -> table1_rows  (donor channel, 450-80 B)
  #   4. The next such block after table1 ends  -> table2_rows  (acceptor
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
      # Must start at 1 and be strictly consecutive
      vals[1L] == 1L && identical(vals, seq_len(length(vals)))
    }
    
    # Step 1: find column-number header row
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
    
    # Step 2: expected row labels - honour plate_format override if supplied
    expected_labels <- if (!is.null(plate_format)) {
      if (plate_format == "96")  LETTERS[1:8]
      else if (plate_format == "384") LETTERS[1:16]
      else stop("plate_format must be '96', '384', or NULL (auto-detect).")
    } else {
      if (n_data_cols <= 12L) LETTERS[1:8] else LETTERS[1:16]
    }
    
    # Step 3-5: find contiguous A-H/A-P blocks; warn if more than two exist
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
          break   # end of contiguous block
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
    
    list(
      col_header_row = col_hdr_row,
      table1_rows    = table1_rows,
      table2_rows    = table2_rows,
      n_data_cols    = n_data_cols
    )
  }
  
  # -- Validate plate_format --------------------------------------------------------
  if (!is.null(plate_format) && !plate_format %in% c("96", "384"))
    stop("plate_format must be '96', '384', or NULL (auto-detect).")
  
  # -- Detect layout --------------------------------------------------------
  layout         <- .locate_plate_tables(data, plate_format = plate_format)
  col_header_row <- layout$col_header_row
  table1_rows    <- layout$table1_rows
  table2_rows    <- layout$table2_rows
  n_data_cols    <- layout$n_data_cols
  
  # Determine the effective plate format for the verbose message
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
  
  # -- Apply column names from header row -------------------------------------
  colnames(data) <- as.character(data[col_header_row, ])
  
  # Total columns: 1 row-label col + n_data_cols
  total_cols <- n_data_cols + 1L
  
  # -- Extract subtables --------------------------------------------------------
  subtable1_full <- data[table1_rows, seq_len(total_cols)]   # donor channel
  subtable2_full <- data[table2_rows, seq_len(total_cols)]   # acceptor channel
  
  final_rownames <- subtable1_full[, 1]
  
  # -- Column selection --------------------------------------------------------
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
              " columns selected). This may cause issues with split_replicates.")
    
    columns_to_keep <- c(1, data_columns_indices)
    subtable1 <- subtable1_full[, columns_to_keep, drop = FALSE]
    subtable2 <- subtable2_full[, columns_to_keep, drop = FALSE]
  } else {
    subtable1 <- subtable1_full
    subtable2 <- subtable2_full
  }
  
  # -- Convert to numeric --------------------------------------------------------
  convert_to_numeric_df <- function(df, rownames_vec) {
    num_df <- as.data.frame(apply(df[, -1, drop = FALSE], 2, as.numeric))
    rownames(num_df) <- rownames_vec
    num_df
  }
  
  subtable1_num <- convert_to_numeric_df(subtable1, final_rownames)
  subtable2_num <- convert_to_numeric_df(subtable2, final_rownames)
  
  # -- Control column mapping --------------------------------------------------------
  map_control_column <- function(control_spec) {
    if (is.null(control_spec)) return(NULL)
    
    if (is.numeric(control_spec)) {
      actual_col_index <- control_spec + 1
      if (actual_col_index > ncol(data))
        stop("Control column index ", control_spec, " is out of bounds. ",
             "Maximum allowed: ", ncol(data) - 1)
      col_name <- as.character(data[col_header_row, actual_col_index])
      return(list(index = actual_col_index, name = col_name, user_index = control_spec))
    } else {
      return(list(index = which(colnames(data) == control_spec)[1],
                  name  = control_spec, user_index = NA))
    }
  }
  
  control_0_info   <- map_control_column(control_0perc)
  control_100_info <- map_control_column(control_100perc)
  
  if (!is.null(control_0_info) && !control_0_info$name %in% colnames(subtable1_num))
    stop("Control column '", control_0_info$name, "' not found in selected columns.")
  
  if (!is.null(control_100_info) && !control_100_info$name %in% colnames(subtable1_num))
    stop("Control column '", control_100_info$name, "' not found in selected columns.")
  
  # -- Low-value filtering --------------------------------------------------------
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
  
  # -- Core ratio calculation --------------------------------------------------------
  ratio          <- (subtable2_num / subtable1_num) * 1000
  ratio_modified <- ratio
  result         <- list()
  
  # -- Info table processing --------------------------------------------------------
  if (!is.null(info_table)) {
    if (ncol(info_table) < 4)
      stop("Info table must have at least 4 columns: log(inhibitor), Plate_Row, Construct, Compound")
    
    base_id_values <- paste(info_table[[3]], info_table[[4]], sep = ":")
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
    
    info_table$ID <- paste(info_table$Construct_Modified, info_table[[4]], sep = ":")
  }
  
  # -- Quality control calculations --------------------------------------------------------
  if (!is.null(info_table) && !is.null(control_0_info) && !is.null(control_100_info)) {
    
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
    
    mean_columns     <- c(control_0_info$name, control_100_info$name)
    existing_columns <- mean_columns[mean_columns %in% colnames(ratio)]
    
    if (length(existing_columns) > 0 && length(row_intervals) > 0) {
      
      general_means <- colMeans(ratio[, existing_columns, drop = FALSE], na.rm = TRUE)
      result$general_means <- data.frame(
        Type   = "General",
        Column = names(general_means),
        Mean   = as.numeric(general_means),
        Rows   = paste0("All (1-", length(table1_rows), ")"),
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
        
        if (length(existing_columns) == 2) {
          d0   <- ratio[valid_rows, control_0_info$name]
          d100 <- ratio[valid_rows, control_100_info$name]
          m0   <- mean(d0,   na.rm = TRUE); s0  <- sd(d0,   na.rm = TRUE)
          m100 <- mean(d100, na.rm = TRUE); s100 <- sd(d100, na.rm = TRUE)
          
          z_score <- if (!is.na(m100) && !is.na(m0) && (m100 - m0) != 0)
            1 - (3 * (s100 + s0) / (m100 - m0)) else NA
          
          aw <- if (!is.na(m100) && !is.na(m0) && m0 != 0) m100 / m0 else NA
          
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
        }
        
        mean_bg  <- if (length(existing_columns) >= 1) mean(ratio[valid_rows, control_0_info$name],   na.rm = TRUE) else NA
        mean_pos <- if (length(existing_columns) >= 2) mean(ratio[valid_rows, control_100_info$name], na.rm = TRUE) else NA
        
        if (length(existing_columns) >= 1) ratio_modified[valid_rows, control_0_info$name]   <- mean_bg
        if (length(existing_columns) >= 2) ratio_modified[valid_rows, control_100_info$name] <- mean_pos
        
        interval_means_list[[cn]] <- data.frame(
          Type                      = "Construct_Interval",
          Construct                 = cn,
          Average_Background        = mean_bg,
          SD_Background             = if (length(existing_columns) >= 1) sd(ratio[valid_rows, control_0_info$name],   na.rm = TRUE) else NA,
          Average_Positive_Ctrl     = mean_pos,
          SD_Positive_Ctrl          = if (length(existing_columns) >= 2) sd(ratio[valid_rows, control_100_info$name], na.rm = TRUE) else NA,
          Average_luciferase_signal = mean_luc,
          Luciferase_signal_comment = luc_comment,
          Z_Score                   = z_score,
          Assay_z_Comment           = zs_comment,
          Assay_Window              = aw,
          Assay_window_Comment      = aw_comment,
          Overall_Quality           = get_lowest_comment(luc_comment, aw_comment, zs_comment),
          Rows                      = paste0(cn, " (rows ", paste(range(valid_rows), collapse = "-"), ")"),
          Rows_Count                = length(valid_rows),
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
          "Assay_z_Comment", "Overall_Quality", "Rows", "Rows_Count"
        )]
        im_clean <- result$interval_means[, -c(1, 2)]
        rownames(im_clean) <- result$interval_means$Construct
        result$interval_means <- as.data.frame(t(im_clean))
        result$construct_intervals <- row_intervals
      }
    }
  }
  
  # -- Column reorganisation --------------------------------------------------------
  if (!is.null(control_0_info) && !is.null(control_100_info)) {
    ctrl_cols <- c(control_0_info$name, control_100_info$name)
    missing   <- ctrl_cols[!ctrl_cols %in% colnames(ratio_modified)]
    if (length(missing) > 0)
      stop("Control columns not found: ", paste(missing, collapse = ", "))
    other_cols     <- setdiff(colnames(ratio_modified), ctrl_cols)
    ratio_modified <- ratio_modified[, c(control_0_info$name, other_cols, control_100_info$name)]
  }
  
  # -- Transpose -------------------------------------------------------------
  ratio_t           <- as.data.frame(t(ratio_modified))
  colnames(ratio_t) <- rownames(ratio_modified)
  
  # -- Rename columns with info_table IDs -------------------------------------
  if (!is.null(info_table)) {
    mapping      <- setNames(info_table$ID, info_table[[2]])
    new_colnames <- mapping[colnames(ratio_t)]
    colnames(ratio_t) <- ifelse(is.na(new_colnames), colnames(ratio_t), new_colnames)
  }
  
  # -- Split technical replicates ---------------------------------------------
  final_table <- if (split_replicates) {
    split_replicates_func <- function(df) {
      n <- nrow(df)
      if (n < 3) { warning("Not enough rows to split replicates."); return(df) }
      ctrl  <- c(1, n)
      exp   <- 2:(n - 1)
      sp    <- floor(length(exp) / 2)
      r1    <- exp[1:sp]
      r2    <- exp[(sp + 1):length(exp)]
      out   <- data.frame()
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
  
  # -- Add log(inhibitor) column ----------------------------------------------
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
      # No guaranteed ctrl structure (split_replicates=FALSE or controls not specified).
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
  
  # -- Excel export --------------------------------------------------------
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
  
  # -- Assemble result --------------------------------------------------------
  result$original_ratio_table <- ratio
  result$modified_ratio_table <- final_table
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
