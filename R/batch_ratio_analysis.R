#' Batch NanoBRET Ratio Dose-Response Analysis
#'
#' Processes an entire directory of raw BMG PHERAstar Excel files in one
#' call: reads each plate, computes BRET ratios via
#' \code{\link{ratio_dose_response_v2}} (or the legacy
#' \code{ratio_dose_response}), and optionally writes per-plate and
#' summary Excel reports.
#'
#' @description
#' This is the \strong{top-level entry point} for the NanoBRET processing
#' pipeline.  It pairs each numbered sheet in the info file with the
#' matching raw data file (matched by the trailing plate number in the
#' filename), calls the chosen ratio function, and collects results into a
#' named list suitable for passing directly to
#' \code{\link{rout_outliers_batch}}.
#'
#' @section File matching:
#' The info file (default \code{info_tables.xlsx}) must contain one sheet
#' per plate, named with a trailing integer (e.g. \code{Sheet1},
#' \code{Plate_01}, \code{2}).  Each sheet is matched to a raw data file
#' whose name ends with \code{_<number>.xlsx} (controlled by
#' \code{data_pattern}).  If multiple files match, the first is used and a
#' warning is issued.
#'
#' @section Raw file parsing:
#' Raw files are parsed by the internal helper \code{.read_nanobret_raw()},
#' which uses a two-stage detection strategy:
#' \enumerate{
#'   \item \strong{Primary path (BMG PHERAstar):} searches for rows matching
#'     \code{"\\d+\\.\\s*Raw Data"} in column 2 (e.g.
#'     \code{"1. Raw Data (450-80 B)"}).  This is the fast, instrument-specific
#'     path and is tried first.
#'   \item \strong{Fallback path (generic):} if no PHERAstar title row is found,
#'     scans the first 30 rows for a row whose integers form a consecutive
#'     sequence starting at 1 (e.g. \code{1, 2, ..., 8} or \code{1, 2, ..., 24}).
#'     This handles partial plates of any width and works with any plate reader
#'     that exports a standard layout (e.g. Tecan Spark, BioTek Synergy).
#' }
#'
#' @param directory Character string.  Directory containing the raw data
#'   files and the info file (default: current working directory).
#'
#' @param control_0perc 0\% control specification passed to the ratio
#'   function.  For \code{function_version = "v1"}: a column name string.
#'   For \code{function_version = "v2"}: either a single numeric value
#'   (fixed-value mode) or a column name string.
#'
#' @param control_100perc 100\% control specification passed to the ratio
#'   function.  For \code{function_version = "v1"}: a column name string.
#'   For \code{function_version = "v2"}: a numeric vector of column
#'   positions or a column name string.
#'
#' @param split_replicates Logical.  Split technical replicates into
#'   separate \code{.2} columns (default \code{TRUE}).  Passed to the
#'   ratio function.
#'
#' @param low_value_threshold Numeric.  Donor-channel values below this
#'   threshold are replaced with \code{NA} before ratio calculation
#'   (default \code{1000}).  Passed to the ratio function.
#'
#' @param info_file Character string.  Name of the info Excel file
#'   (default \code{"info_tables.xlsx"}).  Must reside in \code{directory}.
#'
#' @param data_pattern Regular expression used to identify raw data files
#'   in \code{directory} (default \code{"_\\\\d+\\.xlsx$"}).  Files
#'   whose names begin with \code{~$} (temporary Excel lock files) are
#'   always excluded.
#'
#' @param output_dir Character string.  Directory for output files.
#'   Defaults to \code{directory}.  Created automatically if it does not
#'   exist.
#'
#' @param generate_reports Logical.  Save per-plate Excel files
#'   (\code{results_<N>.xlsx}) and a consolidated
#'   \code{batch_analysis_report.xlsx} (default \code{TRUE}).
#'
#' @param function_version Character.  \code{"v1"} to call the legacy
#'   \code{ratio_dose_response()} or \code{"v2"} (default) to call
#'   \code{\link{ratio_dose_response_v2}}.
#'
#' @param verbose Logical.  Print progress messages (default \code{TRUE}).
#'
#' @param selected_columns Integer vector of 1-based data column indices to
#'   pass to the ratio function (\code{NULL} uses all 24 columns).
#'
#' @param plate_format Character override for plate format: \code{"96"} (rows
#'   A-H) or \code{"384"} (rows A-P).  \code{NULL} (default) auto-detects from
#'   the number of column-header integers.  Use this only when auto-detection
#'   fails, e.g. a 384-well plate read with only 12 columns selected.
#'
#' @param file_map Named list or \code{NULL}.  Optional explicit mapping from
#'   info-sheet name to data filename, bypassing the automatic
#'   \code{_<N>.xlsx} pattern matching.  Filenames are resolved relative to
#'   \code{directory}.  Sheets not listed here still use auto-discovery.
#'   Default \code{NULL} uses auto-discovery for all sheets.
#'
#' @return A named list with one element per successfully processed plate.
#'   The element name is the info-sheet name (e.g. \code{"Sheet1"}).
#'   Each element is a list with:
#'   \describe{
#'     \item{\code{data_file}}{Filename of the raw data file used.}
#'     \item{\code{info_sheet}}{Name of the info sheet used.}
#'     \item{\code{sheet_number}}{Plate number extracted from the sheet
#'       name.}
#'     \item{\code{function_version}}{\code{"v1"} or \code{"v2"}.}
#'     \item{\code{control_0perc}}{Value passed as \code{control_0perc}.}
#'     \item{\code{control_100perc}}{Value passed as
#'       \code{control_100perc}.}
#'     \item{\code{selected_columns}}{Value passed as
#'       \code{selected_columns}.}
#'     \item{\code{result}}{The full return value of the ratio function
#'       (see \code{\link{ratio_dose_response_v2}} for the structure).}
#'   }
#'   Plates that fail processing are omitted from the list and a warning is
#'   issued.  If no plates succeed, a warning is issued and an empty list
#'   is returned.
#'
#' @section Output files (when generate_reports = TRUE):
#' \describe{
#'   \item{\code{results_<N>.xlsx}}{Per-plate workbook with sheets:
#'     \code{Modified_Ratio_Table}, \code{Original_Ratio_Table},
#'     \code{General_Means}, \code{Interval_Means},
#'     \code{Control_Info}, \code{Selected_Columns_Info},
#'     \code{Processing_Info}.}
#'   \item{\code{batch_analysis_report.xlsx}}{Consolidated workbook with a
#'     \code{Summary} sheet (one row per plate) and one quality-metrics
#'     sheet per plate.}
#' }
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' # Bundled example plates + info table live in inst/extdata/.
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_ratio")
#' dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
#' invisible(file.copy(
#'   list.files(extdata_dir, pattern = "^nanobret_", full.names = TRUE),
#'   work_dir, overwrite = TRUE
#' ))
#'
#' ratio_res <- batch_ratio_analysis(
#'   directory        = work_dir,
#'   info_file        = "nanobret_info.xlsx",
#'   data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
#'   control_0perc    = "1",
#'   control_100perc  = "24",
#'   selected_columns = 1:24,
#'   generate_reports = FALSE,
#'   output_dir       = tempdir(),
#'   verbose          = FALSE
#' )
#'
#' names(ratio_res)
#' head(ratio_res$plate_01$result$modified_ratio_table[, 1:4])
#' }
#' @seealso
#' \code{\link{ratio_dose_response_v2}} for the single-plate ratio
#' function called internally.
#'
#' \code{\link{rout_outliers_batch}} for the next step in the
#' pipeline.
#'
#' @export

batch_ratio_analysis <- function(directory = getwd(),
                                 control_0perc = NULL,
                                 control_100perc = NULL,
                                 split_replicates = TRUE,
                                 low_value_threshold = 1000,
                                 info_file = "info_tables.xlsx",
                                 data_pattern = "_\\d+\\.xlsx$",
                                 output_dir = NULL,
                                 generate_reports = TRUE,
                                 function_version = "v1",
                                 verbose = TRUE,
                                 selected_columns = NULL,
                                 plate_format = NULL,
                                 file_map = NULL) {
  
  # -- Validate function_version ----------------------------------------------
  valid_versions <- c("v1", "v2")
  if (!function_version %in% valid_versions)
    stop("function_version must be either 'v1' or 'v2'")
  
  # -- Validate plate_format --------------------------------------------------
  if (!is.null(plate_format) && !plate_format %in% c("96", "384"))
    stop("plate_format must be '96', '384', or NULL (auto-detect).")

  # -- Validate file_map ------------------------------------------------------
  if (!is.null(file_map)) {
    if (!is.list(file_map) || is.null(names(file_map)))
      stop("file_map must be a named list, e.g. list('Sheet1' = 'plate_A.xlsx').")
  }
  
  # -- Internal helper: parse a raw plate-reader Excel file -------------------
  #
  # Reads the file and returns a clean data frame starting from the
  # column-number header row.  No metadata padding or fixed row indices are
  # used - the ratio functions detect the plate layout entirely by content.
  #
  # Detection strategy (two-stage):
  #   1. Primary path (BMG PHERAstar): search for "Raw Data" title rows in
  #      column 2 (e.g. "1. Raw Data (450-80 B)").  Fast and instrument-specific.
  #   2. Fallback path (generic): if no title row is found, scan the first
  #      30 rows for a row whose integers form a consecutive sequence starting
  #      at 1 (e.g. 1-8, 1-12, 1-24).  Works with
  #      any plate reader that exports a standard plate layout.
  #
  # When validate = TRUE, the detected A-H/A-P blocks are checked for correct
  # row labels using the same content-based logic as .locate_plate_tables().
  .read_nanobret_raw <- function(path, plate_format = NULL, validate = TRUE, verbose = FALSE) {
    
    if (!file.exists(path))
      stop("File not found: ", path)
    
    # 1. Read entire sheet with no layout assumptions
    raw <- openxlsx::read.xlsx(path, sheet = 1L,
                               colNames      = FALSE,
                               skipEmptyRows = FALSE,
                               skipEmptyCols = FALSE)
    
    n_raw <- nrow(raw)
    
    # 2. Try primary path: BMG PHERAstar title-row detection
    col2       <- as.character(raw[, 2L])
    title_rows <- grep("^\\d+\\.\\s*Raw Data", col2)
    phera_ok   <- length(title_rows) > 0L
    col_hdr_row_raw <- NA_integer_
    
    if (phera_ok) {
      col_hdr_row_raw <- title_rows[1L] + 1L   # col-number header immediately follows title
      
      if (verbose)
        message("  [read] found ", length(title_rows), " emission table(s) in '",
                basename(path), "'")
      
    } else {
      # 3. Fallback path: scan for a row whose integers form a consecutive
      #    sequence starting at 1 (e.g. 1,2,...,8 or 1,2,...,24).
      #    This handles partial plates of any width while rejecting data rows
      #    with random small integers.
      if (verbose)
        message("  [read fallback] no 'Raw Data' title found in '",
                basename(path), "' - trying generic plate-layout detection.")
      
      .is_col_header <- function(row_vec) {
        vals <- suppressWarnings(as.integer(as.character(row_vec)))
        vals <- sort(vals[!is.na(vals)])
        if (length(vals) < 2L) return(FALSE)
        vals[1L] == 1L && identical(vals, seq_len(length(vals)))
      }
      
      for (i in seq_len(min(n_raw, 30L))) {
        if (.is_col_header(raw[i, ])) {
          col_hdr_row_raw <- i
          break
        }
      }
      
      if (is.na(col_hdr_row_raw))
        stop(
          "Could not locate the plate column-number header in '",
          basename(path), "'.\n",
          "Primary path: no 'Raw Data' title rows found in column 2.\n",
          "Fallback path: no row with consecutive integers starting at 1 ",
          "(e.g. 1-12 or 1-24) found in the first 30 rows.\n",
          "Please check that the file is a valid plate-reader export."
        )
      
      if (verbose)
        message(sprintf(
          "  [read fallback] column-number header found at raw row %d in '%s'",
          col_hdr_row_raw, basename(path)))
    }
    
    # 4. Slice from the column-number header row onward; strip all-NA rows
    result <- raw[col_hdr_row_raw:n_raw, , drop = FALSE]
    all_na <- apply(result, 1L, function(r) all(is.na(r)))
    result <- result[!all_na, , drop = FALSE]
    rownames(result) <- NULL
    
    # 5. Content-based validation of the two A-H/A-P data blocks
    if (validate) {
      
      # Determine plate format: honour plate_format override if supplied,
      # otherwise infer from the number of column-header integers.
      hdr_vals    <- suppressWarnings(as.integer(as.character(result[1L, ])))
      n_data_cols <- sum(!is.na(hdr_vals) & hdr_vals >= 1L & hdr_vals <= 24L)
      expected_labels <- if (!is.null(plate_format)) {
        if (plate_format == "96")       LETTERS[1:8]
        else if (plate_format == "384") LETTERS[1:16]
        else stop("plate_format must be '96', '384', or NULL (auto-detect).")
      } else {
        if (n_data_cols <= 12L) LETTERS[1:8] else LETTERS[1:16]
      }
      
      # Locate the two contiguous label blocks (mirrors .locate_plate_tables)
      find_block <- function(start_row) {
        if (start_row > nrow(result)) return(NULL)   # guard: nothing left to scan
        bs <- NA_integer_; be <- NA_integer_
        for (i in seq(start_row, nrow(result))) {
          lbl <- trimws(as.character(result[i, 1L]))
          if (lbl %in% expected_labels) {
            if (is.na(bs)) bs <- i
            be <- i
          } else if (!is.na(bs)) break
        }
        if (is.na(bs)) return(NULL)
        seq(bs, be)
      }
      
      blk1 <- find_block(2L)
      blk2 <- if (!is.null(blk1)) find_block(max(blk1) + 1L) else NULL
      
      val_ok <- TRUE
      msgs   <- character(0L)
      
      check_block <- function(blk, blk_name) {
        if (is.null(blk)) {
          msgs <<- c(msgs, paste(blk_name, ": block not found"))
          val_ok <<- FALSE
          return()
        }
        actual   <- trimws(as.character(result[blk, 1L]))
        expected <- expected_labels[seq_along(blk)]
        bad      <- which(actual != expected)
        if (length(bad) > 0L) {
          msgs <<- c(msgs, sprintf(
            "%s: position(s) %s have label(s) '%s' (expected '%s')",
            blk_name,
            paste(bad, collapse = ", "),
            paste(actual[bad], collapse = ", "),
            paste(expected[bad], collapse = ", ")))
          val_ok <<- FALSE
        }
      }
      
      check_block(blk1, "Table 1 (donor)")
      check_block(blk2, "Table 2 (acceptor)")
      
      if (!val_ok)
        stop(
          "Row label validation failed in '", basename(path), "'.\n",
          "Expected row labels ", paste(expected_labels, collapse = "/"),
          " in both emission tables.\n",
          paste(msgs, collapse = "\n")
        )
      
      
    }
    
    result
  }   # end .read_nanobret_raw
  
  # -- Internal helper: generate consolidated Excel report --------------------
  .generate_batch_report <- function(results, directory, verbose = TRUE) {
    
    report_path <- file.path(directory, "batch_analysis_report.xlsx")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Summary")
    
    summary_data <- data.frame(
      Sheet_Name = character(), Sheet_Number = character(),
      Data_File = character(), Function_Version = character(),
      Selected_Columns = character(), Status = character(),
      N_Constructs = integer(), Overall_Quality = character(),
      Control_Structure = character(), stringsAsFactors = FALSE
    )
    
    for (plate_sheet in names(results)) {
      plate_result      <- results[[plate_sheet]]$result
      n_constructs      <- 0L
      overall_quality   <- "Not assessed"
      control_structure <- "Unknown"
      selected_cols_info <- "All"
      
      if (!is.null(plate_result$control_info)) {
        if ("new_columns_created" %in% names(plate_result$control_info) &&
            !is.null(plate_result$control_info$new_columns_created))
          control_structure <- paste(plate_result$control_info$new_columns_created, collapse = ", ")
        else
          control_structure <- "Original columns"
      }
      
      if (!is.null(plate_result$selected_columns_info))
        selected_cols_info <- plate_result$selected_columns_info$description
      
      if (!is.null(plate_result$interval_means)) {
        n_constructs <- ncol(plate_result$interval_means)
        if ("Overall_Quality" %in% rownames(plate_result$interval_means)) {
          qualities <- as.character(plate_result$interval_means["Overall_Quality", ])
          overall_quality <- paste(unique(qualities), collapse = ", ")
        }
      }
      
      summary_data <- rbind(summary_data, data.frame(
        Sheet_Name = plate_sheet,
        Sheet_Number = results[[plate_sheet]]$sheet_number,
        Data_File = results[[plate_sheet]]$data_file,
        Function_Version = results[[plate_sheet]]$function_version,
        Selected_Columns = selected_cols_info,
        Status = "Completed",
        N_Constructs = n_constructs,
        Overall_Quality = overall_quality,
        Control_Structure = control_structure,
        stringsAsFactors = FALSE
      ))
    }
    
    openxlsx::writeData(wb, "Summary", summary_data)
    header_style <- openxlsx::createStyle(
      fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center",
      textDecoration = "Bold", border = "TopBottom", borderColour = "#4F81BD"
    )
    openxlsx::addStyle(wb, "Summary", header_style, rows = 1, cols = 1:ncol(summary_data))
    
    for (plate_sheet in names(results)) {
      plate_result <- results[[plate_sheet]]$result
      if (!is.null(plate_result$interval_means)) {
        openxlsx::addWorksheet(wb, plate_sheet)
        info_text <- c(
          paste("Data file:",        results[[plate_sheet]]$data_file),
          paste("Info sheet:",       results[[plate_sheet]]$info_sheet),
          paste("Function version:", results[[plate_sheet]]$function_version),
          paste("Control 0%:",  ifelse(is.numeric(results[[plate_sheet]]$control_0perc),
                                       paste("Fixed value", results[[plate_sheet]]$control_0perc),
                                       results[[plate_sheet]]$control_0perc)),
          paste("Control 100%:", ifelse(is.numeric(results[[plate_sheet]]$control_100perc),
                                        paste("Positions", paste(results[[plate_sheet]]$control_100perc, collapse = ", ")),
                                        paste(results[[plate_sheet]]$control_100perc, collapse = ", "))),
          paste("Selected columns:", ifelse(!is.null(results[[plate_sheet]]$selected_columns),
                                            paste(results[[plate_sheet]]$selected_columns, collapse = ", "),
                                            "All")),
          ""
        )
        openxlsx::writeData(wb, plate_sheet, info_text, startRow = 1)
        quality_data <- as.data.frame(t(plate_result$interval_means))
        quality_data <- cbind(Construct = rownames(quality_data), quality_data)
        rownames(quality_data) <- NULL
        openxlsx::writeData(wb, plate_sheet, quality_data, startRow = length(info_text) + 1)
        openxlsx::addStyle(wb, plate_sheet, header_style,
                           rows = length(info_text) + 1, cols = 1:ncol(quality_data))
      }
    }
    
    openxlsx::saveWorkbook(wb, report_path, overwrite = TRUE)
    
    if (verbose) {
      message("\n=== BATCH COMPLETE ===")
      message("Plates processed: ", length(results))
      message("Report saved:     ", report_path)
      message("======================\n")
    }
  }   # end .generate_batch_report
  
  # ========== MAIN FUNCTION BODY ==========
  
  
  
  if (is.null(output_dir)) {
    output_dir <- directory
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }
  
  if (!dir.exists(directory))
    stop("Directory not found: ", directory)
  
  info_path <- file.path(directory, info_file)
  if (!file.exists(info_path))
    stop("Info file not found: ", info_path)
  
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required. Please install it.")
  
  info_sheets   <- openxlsx::getSheetNames(info_path)
  number_sheets <- info_sheets[grepl("\\d+$", info_sheets)]
  
  if (length(number_sheets) == 0L)
    stop("No valid number sheets found in info file.")
  
  if (verbose) {
    message("Found ", length(number_sheets), " plate sheet(s): ",
            paste(number_sheets, collapse = ", "))
  }
  

  # -- Warn about file_map entries that don't match any info sheet ------------
  if (!is.null(file_map)) {
    unknown_sheets <- setdiff(names(file_map), number_sheets)
    if (length(unknown_sheets) > 0L)
      warning("file_map contains sheet name(s) not found in info file: ",
              paste(unknown_sheets, collapse = ", "))
  }
  data_files <- list.files(directory, pattern = data_pattern, full.names = FALSE)
  data_files <- data_files[!grepl("^~\\$", data_files)]
  
  results <- list()
  
  for (info_sheet in number_sheets) {
    
    sheet_number    <- gsub("^.*?(\\d+)$", "\\1", info_sheet)
    pattern         <- paste0("_", sheet_number, "\\.xlsx$")

    # -- File resolution: explicit file_map takes priority over auto-discovery --
    if (!is.null(file_map) && info_sheet %in% names(file_map)) {
      data_filename <- file_map[[info_sheet]]
      if (!file.exists(file.path(directory, data_filename)))
        stop(sprintf("file_map entry for '%s': file not found: %s",
                     info_sheet, data_filename))
      if (verbose) message(sprintf("\nProcessing %s - %s  [file_map]",
                                   info_sheet, data_filename))
    } else {
      matching_files  <- data_files[grepl(pattern, data_files)]
      if (length(matching_files) == 0L) {
        if (verbose) {
          files_in_dir <- if (length(data_files) > 0L)
            paste(data_files, collapse = ", ") else "(no .xlsx files found)"
          message(sprintf(
            "[skip] Sheet '%s' -- no file matching '_%s.xlsx' found.\n       Files in directory: %s\n       Tip: rename your data file to end with '_%s.xlsx' (e.g. 'plate_experiment_%s.xlsx')",
            info_sheet, sheet_number, files_in_dir, sheet_number, sheet_number))
        }
        next
      }
      data_filename <- matching_files[1L]
      if (length(matching_files) > 1L)
        warning("Multiple files match sheet ", info_sheet, ". Using: ", data_filename)
      if (verbose) message("\nProcessing ", info_sheet, " - ", data_filename)
    }

    data_path <- file.path(directory, data_filename)
    tryCatch({
      
      info_table <- openxlsx::read.xlsx(info_path, sheet = info_sheet)
      
      # Parse raw data - primary PHERAstar path with generic fallback;
      # validates row labels using content-based detection before passing
      # to the ratio function.
      data <- .read_nanobret_raw(data_path,
                                 plate_format = plate_format,
                                 validate     = TRUE,
                                 verbose      = verbose)
      
      if (function_version == "v1") {
        result <- ratio_dose_response(
          data                = data,
          control_0perc       = control_0perc,
          control_100perc     = control_100perc,
          split_replicates    = split_replicates,
          info_table          = info_table,
          low_value_threshold = low_value_threshold,
          verbose             = verbose,
          save_to_excel       = NULL,
          selected_columns    = selected_columns,
          plate_format        = plate_format
        )
      } else {
        result <- ratio_dose_response_v2(
          data                = data,
          control_0perc       = control_0perc,
          control_100perc     = control_100perc,
          split_replicates    = split_replicates,
          info_table          = info_table,
          low_value_threshold = low_value_threshold,
          verbose             = verbose,
          save_to_excel       = NULL,
          selected_columns    = selected_columns,
          plate_format        = plate_format
        )
      }
      
      # Save individual plate results to Excel
      if (generate_reports) {
        excel_path <- file.path(output_dir, paste0("results_", sheet_number, ".xlsx"))
        wb <- openxlsx::createWorkbook()
        
        # BUG_7 fix: define header_style once here so it is always in scope,
        # regardless of whether result$interval_means is NULL.
        header_style <- openxlsx::createStyle(
          fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "center",
          textDecoration = "Bold", border = "TopBottom", borderColour = "#4F81BD"
        )
        
        # Quality_Metrics is added first so it opens as the active sheet
        if (!is.null(result$interval_means)) {
          openxlsx::addWorksheet(wb, "Quality_Metrics")
          interval_df <- as.data.frame(result$interval_means)
          interval_with_rownames <- cbind(Metric = rownames(interval_df), interval_df)
          rownames(interval_with_rownames) <- NULL
          openxlsx::writeData(wb, "Quality_Metrics", interval_with_rownames, rowNames = FALSE)
          
          openxlsx::addStyle(wb, "Quality_Metrics", header_style,
                             rows = 1, cols = 1:ncol(interval_with_rownames))
          
          alt_row_style <- openxlsx::createStyle(fgFill = "#F2F2F2")
          n_rows <- nrow(interval_with_rownames)
          if (n_rows > 1L)
            for (row in seq(2L, n_rows, by = 2L))
              openxlsx::addStyle(wb, "Quality_Metrics", alt_row_style,
                                 rows = row, cols = 1:ncol(interval_with_rownames))
        }
        
        openxlsx::addWorksheet(wb, "Modified_Ratio_Table")
        openxlsx::writeData(wb, "Modified_Ratio_Table",
                            cbind(RowNames = rownames(result$modified_ratio_table),
                                  result$modified_ratio_table),
                            rowNames = FALSE)
        
        openxlsx::addWorksheet(wb, "Original_Ratio_Table")
        openxlsx::writeData(wb, "Original_Ratio_Table",
                            cbind(RowNames = rownames(result$original_ratio_table),
                                  result$original_ratio_table),
                            rowNames = FALSE)
        
        if (!is.null(result$general_means)) {
          openxlsx::addWorksheet(wb, "General_Means")
          openxlsx::writeData(wb, "General_Means", result$general_means)
        }
        
        if (!is.null(result$control_info)) {
          openxlsx::addWorksheet(wb, "Control_Info")
          control_df <- data.frame(
            Parameter = names(result$control_info),
            Value = sapply(result$control_info, function(x) {
              if (is.null(x))        return("NULL")
              if (length(x) == 0L)   return("Empty")
              if (is.character(x))   return(paste(x, collapse = ", "))
              if (is.numeric(x))     return(paste(x, collapse = ", "))
              if (is.logical(x))     return(as.character(x))
              return(as.character(x))
            }),
            stringsAsFactors = FALSE
          )
          openxlsx::writeData(wb, "Control_Info", control_df)
          openxlsx::addStyle(wb, "Control_Info", header_style, rows = 1, cols = 1:2)
        }
        
        if (!is.null(result$selected_columns_info)) {
          openxlsx::addWorksheet(wb, "Selected_Columns_Info")
          selected_cols_df <- data.frame(
            Parameter = names(result$selected_columns_info),
            Value = sapply(result$selected_columns_info, function(x) {
              if (is.null(x))      return("NULL")
              if (length(x) == 0L) return("Empty")
              if (is.character(x)) return(x)
              if (is.numeric(x))   return(paste(x, collapse = ", "))
              return(as.character(x))
            }),
            stringsAsFactors = FALSE
          )
          openxlsx::writeData(wb, "Selected_Columns_Info", selected_cols_df)
          openxlsx::addStyle(wb, "Selected_Columns_Info", header_style, rows = 1, cols = 1:2)
        }
        
        openxlsx::addWorksheet(wb, "Processing_Info")
        processing_info <- data.frame(
          Parameter = c("Data File", "Info Sheet", "Function Version",
                        "Control 0%", "Control 100%", "Selected Columns",
                        "Plate Format", "Split Replicates",
                        "Low Value Threshold", "Processing Date"),
          Value = c(
            data_filename, info_sheet, function_version,
            ifelse(is.numeric(control_0perc),
                   paste("Fixed value:", control_0perc),
                   ifelse(is.null(control_0perc), "Not specified", control_0perc)),
            ifelse(is.numeric(control_100perc),
                   paste("Positions:", paste(control_100perc, collapse = ", ")),
                   ifelse(is.null(control_100perc), "Not specified",
                          paste(control_100perc, collapse = ", "))),
            ifelse(is.null(selected_columns), "All columns",
                   paste("Data columns:", paste(selected_columns, collapse = ", "))),
            ifelse(is.null(plate_format), "Auto-detect", paste0(plate_format, "-well")),
            split_replicates, low_value_threshold, as.character(Sys.time())
          ),
          stringsAsFactors = FALSE
        )
        openxlsx::writeData(wb, "Processing_Info", processing_info)
        openxlsx::addStyle(wb, "Processing_Info", header_style, rows = 1, cols = 1:2)
        
        openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
        if (verbose) message("  Saved: ", basename(excel_path))
      }
      
      results[[info_sheet]] <- list(
        data_file        = data_filename,
        info_sheet       = info_sheet,
        sheet_number     = sheet_number,
        function_version = function_version,
        control_0perc    = control_0perc,
        control_100perc  = control_100perc,
        selected_columns = selected_columns,
        result           = result
      )
      
      if (verbose) message("OK - ", info_sheet, " done")
      
    }, error = function(e) {
      warning("Failed to process sheet ", info_sheet, " (", data_filename,
              ") with ", function_version, ": ", e$message)
      if (verbose) message("X Failed to process sheet ", info_sheet)
    })
  }
  
  if (length(results) > 0L && generate_reports) {
    .generate_batch_report(results, output_dir, verbose)
  } else if (length(results) > 0L) {
    if (verbose) message("Skipping report generation as requested")
  } else {
    warning("No plates were successfully processed")
  }
  
  attr(results, "assay_source") <- "nanobret"
  return(results)
}
