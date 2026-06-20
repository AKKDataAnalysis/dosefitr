#' Batch Viability Analysis from Multiple Plates
#'
#' Performs automated batch processing of viability assay data across multiple
#' Excel files. Each dataset is matched with a corresponding metadata sheet
#' and processed using \code{\link{process_viability_data}}.
#'
#' This function is designed for high-throughput screening experiments where
#' multiple plates must be processed in a standardized and reproducible way.
#'
#' @param directory Character. Path to the directory containing the input Excel
#'   files (both raw data and metadata). Defaults to the current working directory.
#'
#' @param control_0perc Integer (1-24) or NULL. Column index corresponding to the
#'   0\% viability control (e.g., background signal). This is typically a condition
#'   where no viable cells are expected. If NULL, no background normalization is applied.
#'
#' @param control_100perc Integer (1-24) or NULL. Column index corresponding to the
#'   100\% viability control (e.g., untreated cells). This represents the reference
#'   condition for maximum viability. If NULL, normalization is not performed.
#'
#' @param split_replicates Logical. If TRUE, splits experimental concentrations into
#'   two technical replicates after processing. This assumes the plate layout contains
#'   duplicated concentration series. Default is TRUE.
#'
#' @param low_value_threshold Numeric. Minimum allowed value for viability signals.
#'   Values below this threshold are replaced with NA before analysis. Useful for
#'   removing background noise or invalid measurements. Default is 0.
#'
#' @param apply_control_means Logical. If TRUE, replaces individual control values
#'   with construct-specific mean values based on the metadata table. This reduces
#'   variability in control measurements. Default is TRUE.
#'
#' @param auto_detect Logical. If TRUE, automatically detects plate structure
#'   (rows A-P and columns 1-24) within the raw data file. If FALSE, fixed positions
#'   are assumed. Default is TRUE.
#'
#' @param info_file Character. Name of the Excel file containing metadata for all
#'   plates (e.g., constructs, compounds, plate rows). Must be located inside
#'   \code{directory}. Default is "info_tables.xlsx".
#'
#' @param data_pattern Character. Regular expression used to identify data files.
#'   Defaults to \code{"_\\d+\\.xlsx$"}, which matches files ending in an underscore
#'   followed by a number (e.g., "plate_1.xlsx").
#'
#' @param output_dir Character or NULL. Directory where output files will be saved.
#'   If NULL, results are written to the input \code{directory}.
#'
#' @param generate_reports Logical. If TRUE, generates Excel reports for each plate
#'   and a consolidated batch summary file. Default is TRUE.
#'
#' @param selected_columns Integer vector or NULL. Subset of data columns (1-24)
#'   to include in the analysis. These indices refer to the original plate layout
#'   and are passed directly to \code{process_viability_data}. If NULL, all columns
#'   are used.
#'
#' @param verbose Logical. If TRUE, prints progress messages and processing details.
#'   Useful for debugging and tracking batch execution. Default is TRUE.
#'
#' @param file_map Named list or \code{NULL}.  Optional explicit mapping from
#'   info-sheet name to data filename, bypassing the automatic
#'   \code{_<N>.xlsx} pattern matching.  Filenames are resolved relative to
#'   \code{directory}.  Sheets not listed here still use auto-discovery.
#'   Default \code{NULL} uses auto-discovery for all sheets.
#'
#' @param version Character. \code{"v1"} (default) dispatches to
#'   \code{\link{process_viability_data}} (two replicates per row, columns
#'   define the replicate). \code{"v2"} dispatches to
#'   \code{process_viability_data_v2()} (one replicate per row; replicate
#'   rows share \code{Target}+\code{Compound} in the info table; output
#'   column names get \code{.2}, \code{.3}, ... suffixes).
#'   When \code{version = "v2"} and \code{split_replicates = TRUE}, the
#'   wrapper sets \code{split_replicates = FALSE} and emits a one-line
#'   message: in v2 replicates are defined by repeated rows in the info
#'   table, not by row halving.
#'
#' @return A named list where each element corresponds to a processed plate.
#' Each entry contains:
#' \itemize{
#'   \item data_file: Name of the processed data file
#'   \item info_sheet: Metadata sheet used
#'   \item sheet_number: Plate identifier extracted from file/sheet name
#'   \item result: Output list returned by process_viability_data
#' }
#'
#' @details
#' The function iterates over all metadata sheets with numeric suffixes and
#' attempts to match each one to a corresponding data file using the same
#' numeric identifier.
#'
#' For each plate:
#' \itemize{
#'   \item Raw viability data and metadata are loaded from Excel files
#'   \item Data is processed using process_viability_data
#'   \item Quality metrics are computed
#'   \item Optional Excel reports are generated
#' }
#'
#' Quality metrics are calculated per construct using plate row mappings
#' defined in the metadata table.
#'
#' @section Input Requirements:
#' \itemize{
#'   \item Metadata Excel file must contain sheets ending in numeric identifiers
#'   \item Data files must follow the naming pattern defined by \code{data_pattern}
#'   \item Metadata must contain at least four columns:
#'   log(inhibitor), plate row (A-P), construct, and compound
#' }
#'
#' @section Output:
#' \itemize{
#'   \item Per-plate Excel reports (viability_results_<n>.xlsx)
#'   \item Quality metrics files (drc_quality/viability_results_<n>.xlsx)
#'   \item Batch summary report (batch_viability_report.xlsx)
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- batch_viability_analysis(
#'   control_0perc   = 13,
#'   control_100perc = 12,
#'   selected_columns = c(2:23)
#' )
#'
#' # Custom directory and no report generation
#' results <- batch_viability_analysis(
#'   directory = "data/",
#'   generate_reports = FALSE
#' )
#'
#' # --- Explicit file mapping (when files cannot be renamed) ---
#' results <- batch_viability_analysis(
#'   directory       = "data/",
#'   control_0perc   = 13,
#'   control_100perc = 12,
#'   file_map        = list(
#'     "Sheet1" = "viability_plate_A.xlsx",
#'     "Sheet2" = "viability_plate_B.xlsx"
#'   )
#' )
#' }
#'
#' @seealso \code{\link{process_viability_data}}
#'
#' @export


batch_viability_analysis <- function(directory           = getwd(),
                                     control_0perc       = NULL,
                                     control_100perc     = NULL,
                                     split_replicates    = TRUE,
                                     low_value_threshold = 0,
                                     apply_control_means = TRUE,
                                     auto_detect         = TRUE,
                                     info_file           = "info_tables.xlsx",
                                     data_pattern        = "_\\d+\\.xlsx$",
                                     output_dir          = NULL,
                                     generate_reports    = TRUE,
                                     selected_columns    = NULL,
                                     verbose             = TRUE,
                                     file_map            = NULL,
                                     version             = c("v1", "v2")) {
  
  # -- Dependency check -------------------------------------------------------
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required. Please install it.")
  
  # -- Resolve processing version --------------------------------------------
  version <- match.arg(version)

  # -- Verify the requested processor is available ---------------------------
  if (version == "v1") {
    if (!exists("process_viability_data", mode = "function"))
      stop(paste0(
        "process_viability_data() not found. ",
        "Please source the script containing it before calling this function."))
  } else {
    if (!exists("process_viability_data_v2", mode = "function"))
      stop(paste0(
        "process_viability_data_v2() not found. ",
        "Source process_viability_data_v2.R before calling with version = 'v2'."))
  }

  # -- v2: row-based replicates render split_replicates meaningless ----------
  if (version == "v2" && isTRUE(split_replicates)) {
    if (verbose)
      message("split_replicates is ignored in v2; replicates are defined by repeated rows in info_table.")
    split_replicates <- FALSE
  }

  # -- Helper: number of unique compounds in a result table -------------------
  # In v1 each compound appears as a .2 pair, so dividing by 2 is correct.
  # In v2 each replicate is its own column with .2/.3/.../ suffix, so we count
  # unique base names instead.
  n_unique_compounds <- function(modified_tbl) {
    if (is.null(modified_tbl)) return(0L)
    if (version == "v2") {
      data_names <- colnames(modified_tbl)[-1L]
      length(unique(sub("\\.\\d+$", "", data_names)))
    } else {
      (ncol(modified_tbl) - 1L) %/% 2L
    }
  }
  
  # -- Input validation -------------------------------------------------------
  if (!is.null(control_0perc) &&
      !(is.numeric(control_0perc) && length(control_0perc) == 1L &&
        control_0perc >= 1 && control_0perc <= 24))
    stop("control_0perc must be a single integer between 1 and 24.")
  
  if (!is.null(control_100perc) &&
      !(is.numeric(control_100perc) && length(control_100perc) == 1L &&
        control_100perc >= 1 && control_100perc <= 24))
    stop("control_100perc must be a single integer between 1 and 24.")

  # -- Validate file_map ------------------------------------------------------
  if (!is.null(file_map)) {
    if (!is.list(file_map) || is.null(names(file_map)))
      stop("file_map must be a named list, e.g. list('Sheet1' = 'viability_plate_A.xlsx').")
  }

  # -- Directory / output setup -----------------------------------------------
  if (!dir.exists(directory))
    stop("Directory not found: ", directory)
  
  if (is.null(output_dir)) {
    output_dir <- directory
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }
  
  # -- Locate info file -------------------------------------------------------
  info_path <- file.path(directory, info_file)
  if (!file.exists(info_path))
    stop("Info file not found: ", info_path)
  
  # -- Discover data files ----------------------------------------------------
  data_files <- list.files(directory, pattern = data_pattern, full.names = FALSE)
  data_files <- data_files[!grepl("^~\\$", data_files)]
  
  # -- Read info sheet names --------------------------------------------------
  info_sheets   <- openxlsx::getSheetNames(info_path)
  number_sheets <- info_sheets[grepl("\\d+$", info_sheets)]
  
  if (length(number_sheets) == 0L)
    stop("No valid numbered sheets found in info file: ", info_path)
  
  if (verbose)
    message("Found ", length(number_sheets), " plate sheet(s): ",
            paste(number_sheets, collapse = ", "))
  

  # -- Warn about file_map entries that don't match any info sheet ------------
  if (!is.null(file_map)) {
    unknown_sheets <- setdiff(names(file_map), number_sheets)
    if (length(unknown_sheets) > 0L)
      warning("file_map contains sheet name(s) not found in info file: ",
              paste(unknown_sheets, collapse = ", "))
  }
  # -- Internal helper: compute per-construct quality metrics -----------------
  # Viability-appropriate metrics for raw (non-normalized) data.
  #
  # Overall_Quality is driven exclusively by CV% of each control, which are
  # scale-independent and directly reflect normalization anchor reliability.
  # Signal_to_Background is retained as a descriptive number only -- it carries
  # no quality comment and does not contribute to Overall_Quality, because its
  # absolute value depends on instrument/assay signal levels and is not
  # interpretable without context when data are not yet normalized.
  #
  # Metrics computed per construct:
  #   Mean_Background        : mean of 0% control replicates
  #   SD_Background          : SD of 0% control replicates
  #   CV_Background_pct      : CV% of 0% control  (SD/Mean * 100)
  #   Mean_Positive_Ctrl     : mean of 100% control replicates
  #   SD_Positive_Ctrl       : SD of 100% control replicates
  #   CV_Positive_Ctrl_pct   : CV% of 100% control (SD/Mean * 100)
  #   Signal_to_Background   : Mean_Positive / Mean_Background (descriptive only)
  #   CV_Background_Comment  : high (<=10%), medium (10-20%), low (>20%)
  #   CV_PosCtrl_Comment     : high (<=10%), medium (10-20%), low (>20%)
  #   Overall_Quality        : lowest of CV_Background and CV_PosCtrl assessments
  #   Rows                   : construct name + row range (informational)
  #   Rows_Count             : number of replicate rows used
  #
  # Returns a data.frame (constructs as columns, metrics as rows), or NULL when
  # the required inputs (controls + info_table) are not available.
  compute_quality_metrics <- function(result) {
    
    proc_info_list <- result$processing_info
    if (is.null(proc_info_list)) return(NULL)
    
    viability_data   <- proc_info_list$viability_data
    control_0_info   <- proc_info_list$control_0_info
    control_100_info <- proc_info_list$control_100_info
    info_tbl         <- proc_info_list$info_table
    
    if (is.null(viability_data) || is.null(control_0_info) ||
        is.null(control_100_info) || is.null(info_tbl)) return(NULL)
    
    if (!"Construct_Modified" %in% colnames(info_tbl))
      info_tbl$Construct_Modified <- info_tbl[[3]]
    
    plate_row_values  <- info_tbl[[2]]
    construct_values  <- info_tbl$Construct_Modified
    unique_constructs <- unique(construct_values)
    
    # Map plate-row letters (A-P) to numeric row indices 1-16
    row_index_map <- stats::setNames(seq_len(16L), LETTERS[seq_len(16L)])
    
    construct_groups <- lapply(
      stats::setNames(unique_constructs, unique_constructs),
      function(cn) {
        rows <- plate_row_values[construct_values == cn]
        as.integer(row_index_map[rows])
      }
    )
    construct_groups <- Filter(function(x) length(x) > 0L, construct_groups)
    
    ctrl0_col   <- control_0_info$name
    ctrl100_col <- control_100_info$name
    
    if (!ctrl0_col   %in% colnames(viability_data) ||
        !ctrl100_col %in% colnames(viability_data)) return(NULL)
    
    # CV% uses inverted thresholds: lower CV = better quality.
    # Classify by checking upper bounds from best to worst.
    cv_quality_level <- function(cv_val) {
      if (is.na(cv_val)) return("insufficient (NA)")
      if (cv_val <= 10)  return("high (<=10%)")
      if (cv_val <= 20)  return("medium (10-20%)")
      return("low (>20%)")
    }
    
    # Return the lowest quality level across all supplied label strings.
    # Quality order: insufficient < low < medium < high.
    lowest_quality <- function(...) {
      quality_order <- c("insufficient", "low", "medium", "high")
      lvls <- c(...)
      first_words <- vapply(lvls, function(x) {
        if (is.null(x) || is.na(x) || identical(x, "")) return("insufficient")
        strsplit(as.character(x), " ", fixed = TRUE)[[1L]][1L]
      }, character(1L), USE.NAMES = FALSE)
      scores <- match(first_words, quality_order)
      scores[is.na(scores)] <- 1L
      quality_order[min(scores)]
    }
    
    rows_list <- lapply(names(construct_groups), function(cn) {
      valid_rows <- construct_groups[[cn]]
      valid_rows <- valid_rows[
        !is.na(valid_rows) & valid_rows >= 1L & valid_rows <= nrow(viability_data)
      ]
      if (length(valid_rows) == 0L) return(NULL)
      
      bg_vals  <- viability_data[valid_rows, ctrl0_col]
      pos_vals <- viability_data[valid_rows, ctrl100_col]
      
      mean_bg  <- mean(bg_vals,  na.rm = TRUE)
      sd_bg    <- stats::sd(bg_vals,  na.rm = TRUE)
      mean_pos <- mean(pos_vals, na.rm = TRUE)
      sd_pos   <- stats::sd(pos_vals, na.rm = TRUE)
      
      # CV% - guard against division by zero / near-zero mean
      cv_bg  <- if (!is.na(mean_bg)  && abs(mean_bg)  > 1e-9)
        (sd_bg  / mean_bg)  * 100 else NA_real_
      cv_pos <- if (!is.na(mean_pos) && abs(mean_pos) > 1e-9)
        (sd_pos / mean_pos) * 100 else NA_real_
      
      # Signal-to-background ratio
      sb_ratio <- if (!is.na(mean_bg) && abs(mean_bg) > 1e-9)
        mean_pos / mean_bg else NA_real_
      
      # Quality comments - CV% only; S/B has no comment (descriptive only)
      cv_bg_comment  <- cv_quality_level(cv_bg)
      cv_pos_comment <- cv_quality_level(cv_pos)
      
      # Overall driven by CV% of both controls only
      overall <- lowest_quality(cv_bg_comment, cv_pos_comment)
      
      row_range <- paste(LETTERS[range(valid_rows)], collapse = "-")
      
      data.frame(
        Construct              = cn,
        Mean_Background        = round(mean_bg,   3L),
        SD_Background          = round(sd_bg,     3L),
        CV_Background_pct      = round(cv_bg,     2L),
        Mean_Positive_Ctrl     = round(mean_pos,  3L),
        SD_Positive_Ctrl       = round(sd_pos,    3L),
        CV_Positive_Ctrl_pct   = round(cv_pos,    2L),
        CV_Background_Comment  = cv_bg_comment,
        CV_PosCtrl_Comment     = cv_pos_comment,
        Overall_Quality        = overall,
        Signal_to_Background   = round(sb_ratio,  3L),
        Rows                   = paste0(cn, " (", row_range, ")"),
        Rows_Count             = length(valid_rows),
        stringsAsFactors       = FALSE
      )
    })
    
    rows_df <- do.call(rbind, Filter(Negate(is.null), rows_list))
    if (is.null(rows_df) || nrow(rows_df) == 0L) return(NULL)
    
    # Transpose: constructs as columns, metrics as rows
    metrics_cols <- setdiff(colnames(rows_df), "Construct")
    mat          <- as.data.frame(t(rows_df[, metrics_cols, drop = FALSE]))
    colnames(mat) <- rows_df$Construct
    mat
  }
  
  # -- Internal helper: consolidated batch report -----------------------------
  generate_batch_report <- function(results, out_dir) {
    
    quality_dir <- file.path(out_dir, "drc_quality")
    if (!dir.exists(quality_dir))
      dir.create(quality_dir, recursive = TRUE)
    
    report_path  <- file.path(quality_dir, "batch_viability_report.xlsx")
    wb           <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Summary")
    
    header_style <- openxlsx::createStyle(
      fontColour     = "#FFFFFF",
      fgFill         = "#4F81BD",
      halign         = "center",
      textDecoration = "Bold",
      border         = "TopBottom",
      borderColour   = "#4F81BD"
    )
    
    summary_rows <- lapply(names(results), function(sheet_name) {
      entry  <- results[[sheet_name]]
      res    <- entry$result
      
      n_compounds <- n_unique_compounds(res$modified_ratio_table)
      
      data.frame(
        Sheet_Name       = sheet_name,
        Sheet_Number     = entry$sheet_number,
        Data_File        = entry$data_file,
        Control_0perc    = ifelse(is.null(entry$control_0perc),
                                  "not set", as.character(entry$control_0perc)),
        Control_100perc  = ifelse(is.null(entry$control_100perc),
                                  "not set", as.character(entry$control_100perc)),
        Selected_Columns = ifelse(is.null(entry$selected_columns),
                                  "all (1-24)",
                                  paste(entry$selected_columns, collapse = ", ")),
        N_Compounds      = n_compounds,
        Status           = "Completed",
        stringsAsFactors = FALSE
      )
    })
    
    summary_data <- do.call(rbind, summary_rows)
    openxlsx::writeData(wb, "Summary", summary_data)
    openxlsx::addStyle(wb, "Summary", header_style,
                       rows = 1L, cols = seq_len(ncol(summary_data)))
    
    openxlsx::saveWorkbook(wb, report_path, overwrite = TRUE)
    
    
  }
  
  # -- Main processing loop ---------------------------------------------------
  results <- list()
  
  for (info_sheet in number_sheets) {
    
    sheet_number   <- gsub("^.*?(\\d+)$", "\\1", info_sheet)
    pattern        <- paste0("_", sheet_number, "\\.xlsx$")

    # -- File resolution: explicit file_map takes priority over auto-discovery --
    if (!is.null(file_map) && info_sheet %in% names(file_map)) {
      data_filename <- file_map[[info_sheet]]
      if (!file.exists(file.path(directory, data_filename)))
        stop(sprintf("file_map entry for '%s': file not found: %s",
                     info_sheet, data_filename))
      if (verbose) message(sprintf("\nProcessing %s - %s  [file_map]",
                                   info_sheet, data_filename))
    } else {
      matching_files <- data_files[grepl(pattern, data_files)]
      if (length(matching_files) == 0L) {
        if (verbose) {
          files_in_dir <- if (length(data_files) > 0L)
            paste(data_files, collapse = ", ") else "(no .xlsx files found)"
          message(sprintf(
            "[skip] Sheet '%s' -- no file matching '_%s.xlsx' found.\n       Files in directory: %s\n       Tip: rename your data file to end with '_%s.xlsx' (e.g. 'viability_plate_%s.xlsx')",
            info_sheet, sheet_number, files_in_dir, sheet_number, sheet_number))
        }
        next
      }
      data_filename <- matching_files[1L]
      if (length(matching_files) > 1L)
        warning(sprintf("Multiple files match sheet '%s'. Using: %s",
                        info_sheet, data_filename))
      if (verbose) message("\nProcessing ", info_sheet, " - ", data_filename)
    }

    data_path <- file.path(directory, data_filename)
    
    tryCatch({
      
      # -- Read info table for this plate -------------------------------------
      info_table <- openxlsx::read.xlsx(info_path, sheet = info_sheet)
      
      # -- Read raw viability data --------------------------------------------
      raw_data <- openxlsx::read.xlsx(
        data_path,
        sheet         = 1L,
        colNames      = FALSE,
        skipEmptyRows = FALSE,
        skipEmptyCols = FALSE
      )
      
      if (verbose)
        message("  [read] ", nrow(raw_data), " rows x ", ncol(raw_data),
                " cols in '", data_filename, "'")
      
      # -- Call the version-appropriate processor -----------------------------
      result <- if (version == "v2") {
        process_viability_data_v2(
          data                = raw_data,
          control_0perc       = control_0perc,
          control_100perc     = control_100perc,
          split_replicates    = split_replicates,
          info_table          = info_table,
          selected_columns    = selected_columns,
          low_value_threshold = low_value_threshold,
          verbose             = verbose,
          apply_control_means = apply_control_means,
          auto_detect         = auto_detect
        )
      } else {
        process_viability_data(
          data                = raw_data,
          control_0perc       = control_0perc,
          control_100perc     = control_100perc,
          split_replicates    = split_replicates,
          info_table          = info_table,
          selected_columns    = selected_columns,
          low_value_threshold = low_value_threshold,
          verbose             = verbose,
          apply_control_means = apply_control_means,
          auto_detect         = auto_detect
        )
      }
      
      # -- Rename $modified_table -> $modified_ratio_table --------------------
      if (!is.null(result$modified_table)) {
        result$modified_ratio_table <- result$modified_table
        result$modified_table       <- NULL
      }
      
      # -- Compute quality metrics --------------------------------------------
      quality_metrics <- compute_quality_metrics(result)
      
      
      
      # -- Optional per-plate Excel report ------------------------------------
      if (generate_reports) {
        
        quality_dir <- file.path(output_dir, "drc_quality")
        if (!dir.exists(quality_dir))
          dir.create(quality_dir, recursive = TRUE)
        
        excel_path <- file.path(quality_dir,
                                paste0("viability_results_", sheet_number, ".xlsx"))
        wb <- openxlsx::createWorkbook()
        
        header_style <- openxlsx::createStyle(
          fontColour     = "#FFFFFF",
          fgFill         = "#4F81BD",
          halign         = "center",
          textDecoration = "Bold",
          border         = "TopBottom",
          borderColour   = "#4F81BD"
        )
        
        # Sheet 1: Quality metrics (first for quick inspection)
        if (!is.null(quality_metrics)) {
          openxlsx::addWorksheet(wb, "Quality_Metrics")
          openxlsx::writeData(
            wb, "Quality_Metrics",
            cbind(Metric = rownames(quality_metrics), quality_metrics),
            rowNames = FALSE
          )
          openxlsx::addStyle(
            wb, "Quality_Metrics", header_style,
            rows = 1L,
            cols = seq_len(ncol(quality_metrics) + 1L)
          )
        }
        
        # Sheet 2: Modified table
        openxlsx::addWorksheet(wb, "Modified_Table")
        openxlsx::writeData(
          wb, "Modified_Table",
          cbind(RowNames = rownames(result$modified_ratio_table),
                result$modified_ratio_table),
          rowNames = FALSE
        )
        openxlsx::addStyle(
          wb, "Modified_Table", header_style,
          rows = 1L,
          cols = seq_len(ncol(result$modified_ratio_table) + 1L)
        )
        
        # Sheet 3: Original table
        if (!is.null(result$original_table)) {
          openxlsx::addWorksheet(wb, "Original_Table")
          openxlsx::writeData(
            wb, "Original_Table",
            cbind(RowNames = rownames(result$original_table),
                  result$original_table),
            rowNames = FALSE
          )
          openxlsx::addStyle(
            wb, "Original_Table", header_style,
            rows = 1L,
            cols = seq_len(ncol(result$original_table) + 1L)
          )
        }
        
        # Processing info
        openxlsx::addWorksheet(wb, "Processing_Info")
        proc_info_df <- data.frame(
          Parameter = c(
            "Data File", "Info Sheet", "Sheet Number",
            "Control 0%", "Control 100%",
            "Split Replicates", "Low Value Threshold",
            "Apply Control Means", "Auto Detect",
            "Selected Columns", "Processing Date"
          ),
          Value = c(
            data_filename, info_sheet, sheet_number,
            ifelse(is.null(control_0perc),  "not set",
                   as.character(control_0perc)),
            ifelse(is.null(control_100perc), "not set",
                   as.character(control_100perc)),
            as.character(split_replicates),
            as.character(low_value_threshold),
            as.character(apply_control_means),
            as.character(auto_detect),
            ifelse(is.null(selected_columns), "all (1-24)",
                   paste(selected_columns, collapse = ", ")),
            as.character(Sys.time())
          ),
          stringsAsFactors = FALSE
        )
        openxlsx::writeData(wb, "Processing_Info", proc_info_df)
        openxlsx::addStyle(wb, "Processing_Info", header_style,
                           rows = 1L, cols = seq_len(2L))
        
        openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
        
        if (verbose)
          message("  Saved: ", basename(excel_path))
      }
      
      # -- Store result -------------------------------------------------------
      n_compounds <- n_unique_compounds(result$modified_ratio_table)
      
      results[[info_sheet]] <- list(
        data_file        = data_filename,
        info_sheet       = info_sheet,
        sheet_number     = sheet_number,
        control_0perc    = control_0perc,
        control_100perc  = control_100perc,
        selected_columns = selected_columns,
        result           = result
      )
      
      if (verbose)
        message("OK - ", info_sheet, " done")
      
    }, error = function(e) {
      warning(sprintf("Failed to process sheet '%s' (%s): %s",
                      info_sheet, data_filename, e$message))
      if (verbose)
        message("X Failed to process sheet ", info_sheet)
    })
  }
  
  # -- Consolidated report ----------------------------------------------------
  if (length(results) > 0L) {
    if (generate_reports)
      generate_batch_report(results, output_dir)
  } else {
    warning("No plates were successfully processed.")
  }
  
  # -- Final summary ----------------------------------------------------------
  if (verbose) {
    message("\n=== BATCH COMPLETE ===")
    message("Plates processed: ", length(results), " / ", length(number_sheets))
    if (length(results) > 0L && generate_reports)
      message("Report saved:     ",
              file.path(output_dir, "batch_viability_report.xlsx"))
    message("======================\n")
  }
  
  attr(results, "assay_source") <- "viability"
  return(invisible(results))
}
