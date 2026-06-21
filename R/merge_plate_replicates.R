#' Merge Replicate Plates into a Single Combined Result
#'
#' Takes two or more plates from a \code{\link{batch_ratio_analysis}} (or
#' \code{\link{batch_viability_analysis}}) result list and combines their
#' \code{modified_ratio_table}s into a single merged plate entry.  Replicates
#' from each plate are renumbered sequentially so that, for example, two plates
#' each contributing two replicates of \code{LRRK2:MDKM34} produce four
#' columns labelled \code{LRRK2:MDKM34}, \code{LRRK2:MDKM34.2},
#' \code{LRRK2:MDKM34.3}, and \code{LRRK2:MDKM34.4}.
#'
#' The merged entry replaces the individual plate entries in the returned list
#' and is fully compatible with all downstream functions
#' (\code{\link{rout_outliers_batch}}, \code{\link{batch_drc_analysis}},
#' \code{\link{plot_multiple_compounds}}, etc.).
#'
#' @section Concentration matching:
#' All plates being merged must share an identical \code{log(inhibitor)}
#' column (same values in the same order).  If they differ, the function
#' stops with an informative error.  Use \code{check_concentrations = FALSE}
#' only if you are certain the differences are negligible floating-point
#' artefacts and you want to suppress the check.
#'
#' @section Column renumbering:
#' Within each plate, replicate columns are already named with a \code{.2}
#' suffix (e.g. \code{LRRK2:MDKM34} and \code{LRRK2:MDKM34.2}).  When
#' merging, the function strips all existing numeric suffixes, groups columns
#' by base name, and re-assigns suffixes \code{.2}, \code{.3}, \code{.4}, ...
#' across plates in the order the plates are supplied.  The first replicate
#' of the first plate keeps no suffix (i.e. it remains the "primary" column).
#'
#' @section Asymmetric plates:
#' Compounds that appear in only a subset of the merged plates are included
#' in the output with however many replicates they have.  No NA-padding is
#' added for missing plates.
#'
#' @param results Named list returned by \code{\link{batch_ratio_analysis}} or
#'   \code{\link{batch_viability_analysis}}.
#'
#' @param plates Character vector of plate names (element names of
#'   \code{results}) to merge.  Defaults to \code{NULL}, which merges
#'   \strong{all} plates in \code{results}.
#'
#' @param merged_name Character string.  Name of the new merged entry in the
#'   returned list (default \code{"merged"}).
#'   
#' @param output_dir Character string.  Directory where the Excel report is
#'   saved when \code{generate_reports = TRUE}.  A \code{drc_quality/}
#'   sub-folder is created inside this directory (default \code{NULL}, which
#'   uses the current working directory).
#'   
#' @param check_concentrations Logical.  If \code{TRUE} (default), stop with
#'   an error when the \code{log(inhibitor)} columns of the plates being merged
#'   are not identical.
#'
#' @param generate_reports Logical.  Save a
#'   \code{merged_results_<merged_name>.xlsx} workbook inside a
#'   \code{drc_quality/} sub-folder of the current working directory
#'   (default \code{TRUE}).  The folder is created automatically if it does
#'   not exist.  Set to \code{FALSE} to skip writing entirely.
#'
#' @param verbose Logical.  Print progress messages (default \code{TRUE}).
#'
#' @return A named list in the same format as \code{\link{batch_ratio_analysis}}
#'   output.  The individual plate entries that were merged are removed and
#'   replaced by a single entry named \code{merged_name}.  Any plates that
#'   were \emph{not} selected for merging are kept unchanged.
#'
#'   The merged entry's \code{$result$modified_ratio_table} contains all
#'   replicate columns from every merged plate, renumbered sequentially.
#'   \code{$result$original_ratio_table} is set to \code{NULL} (the concept
#'   of a single "original" table does not apply to a merge).
#'   All other \code{$result} slots (\code{general_means},
#'   \code{interval_means}, etc.) are set to \code{NULL}; they will be
#'   recomputed by downstream functions as needed.
#'
#' @section Output files:
#' When \code{generate_reports = TRUE}, a single Excel workbook is written to
#' a \code{drc_quality/} sub-folder of the current working directory
#' (created automatically if it does not exist):
#' \describe{
#'   \item{\code{drc_quality/merged_results_<merged_name>.xlsx}}{Workbook
#'     with sheets: \code{Merged_Table} (combined modified table with all
#'     replicates renumbered), one \code{Original_<plate>} sheet per merged
#'     plate (the \code{modified_ratio_table} from that plate before merging),
#'     and a \code{Provenance} sheet listing which plates were merged and
#'     their source data files.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- batch_ratio_analysis(
#'   control_0perc   = 24,
#'   control_100perc = 12
#' )
#'
#' # Merge all plates (default)
#' merged <- merge_plate_replicates(results)
#'
#' # Merge only two specific plates
#' merged <- merge_plate_replicates(results,
#'   plates      = c("Sheet1", "Sheet2"),
#'   merged_name = "plates_1_2"
#' )
#'
#' # Feed directly into batch_drc_analysis
#' drc_results <- batch_drc_analysis(
#'   batch_results = merged,
#'   normalize     = TRUE,
#'   output_dir    = "./drc_results"
#' )
#' }
#'
#' @seealso
#' \code{\link{batch_ratio_analysis}} for the upstream function whose output
#' this function consumes.
#'
#' \code{\link{batch_drc_analysis}} for the downstream function that accepts
#' the merged result.
#'
#' @export
merge_plate_replicates <- function(results,
                                   plates                = NULL,
                                   merged_name           = "merged",
                                   check_concentrations  = TRUE,
                                   output_dir            = NULL,
                                   generate_reports      = TRUE,
                                   verbose               = TRUE) {
  
  # -- Input validation -------------------------------------------------------
  if (!is.list(results) || length(results) == 0L)
    stop("'results' must be a non-empty named list from batch_ratio_analysis().")
  
  if (is.null(plates)) {
    plates <- names(results)
  } else {
    missing_plates <- setdiff(plates, names(results))
    if (length(missing_plates) > 0L)
      stop("The following plate(s) were not found in 'results': ",
           paste(missing_plates, collapse = ", "),
           "\nAvailable plates: ", paste(names(results), collapse = ", "))
  }
  
  if (length(plates) < 2L)
    stop("At least 2 plates are required for merging. ",
         "Only 1 plate selected: ", paste(plates, collapse = ", "))
  
  if (!is.character(merged_name) || length(merged_name) != 1L || nchar(merged_name) == 0L)
    stop("'merged_name' must be a single non-empty character string.")
  
  # -- Extract modified_ratio_tables ------------------------------------------
  tables <- lapply(plates, function(p) {
    tbl <- results[[p]]$result$modified_ratio_table
    if (is.null(tbl))
      stop("Plate '", p, "' has no $result$modified_ratio_table. ",
           "Has batch_ratio_analysis() been run successfully for this plate?")
    if (!is.data.frame(tbl))
      stop("$result$modified_ratio_table for plate '", p,
           "' is not a data.frame (class: ", paste(class(tbl), collapse = "/"), ").")
    if (nrow(tbl) == 0L)
      stop("$result$modified_ratio_table for plate '", p,
           "' has 0 rows. The plate may have failed to process correctly. ",
           "Check the batch_ratio_analysis() output for warnings on this plate.")
    if (ncol(tbl) < 2L)
      stop("$result$modified_ratio_table for plate '", p,
           "' has fewer than 2 columns (", ncol(tbl), " column(s) found). ",
           "Expected at least a concentration column and one compound column.")
    tbl
  })
  names(tables) <- plates
  
  # -- Identify the log(inhibitor) column -------------------------------------
  # Convention: first column of modified_ratio_table is the concentration column.
  conc_col <- colnames(tables[[1L]])[1L]
  
  if (verbose) {
    cat(strrep("=", 60), "\n")
    cat("MERGE PLATE REPLICATES\n")
    cat(strrep("=", 60), "\n")
    cat(sprintf("Plates to merge  : %s\n", paste(plates, collapse = ", ")))
    cat(sprintf("Merged name      : %s\n", merged_name))
    cat(sprintf("Concentration col: %s\n", conc_col))
    cat(strrep("-", 60), "\n")
  }
  
  # -- Warn if concentration column names differ across plates ----------------
  # Column 1 is always the concentration column regardless of its name.
  conc_col_names <- vapply(tables, function(tbl) colnames(tbl)[1L], character(1L))
  if (length(unique(conc_col_names)) > 1L) {
    warning(
      "Concentration column names differ across plates: ",
      paste(unique(conc_col_names), collapse = ", "), ". ",
      "All plates are still processed correctly using column position (column 1). ",
      "The merged table will use the name from the first plate: '", conc_col, "'."
    )
  }
  
  # -- Check concentration columns are identical ------------------------------
  if (check_concentrations) {
    ref_conc <- tables[[1L]][[1L]]          # always positional
    for (p in plates[-1L]) {
      other_conc <- tables[[p]][[1L]]       # always positional
      if (length(other_conc) != length(ref_conc) ||
          !isTRUE(all.equal(ref_conc, other_conc, tolerance = 1e-9))) {
        stop(
          "Concentration mismatch between plate '", plates[1L],
          "' and plate '", p, "'.\n",
          "  Plate '", plates[1L], "' has ", length(ref_conc),
          " concentration points: ",
          paste(head(ref_conc, 5L), collapse = ", "),
          if (length(ref_conc) > 5L) " ..." else "", "\n",
          "  Plate '", p, "' has ", length(other_conc),
          " concentration points: ",
          paste(head(other_conc, 5L), collapse = ", "),
          if (length(other_conc) > 5L) " ..." else "", "\n",
          "All plates must share an identical log(inhibitor) column.\n",
          "Set check_concentrations = FALSE to suppress this check."
        )
      }
    }
    if (verbose) message("Concentration columns match across all plates.")
  }
  
  # -- Helper: strip trailing numeric suffix from a column name ---------------
  # "LRRK2:MDKM34"    -> "LRRK2:MDKM34"
  # "LRRK2:MDKM34.2"  -> "LRRK2:MDKM34"
  # "LRRK2:MDKM34.10" -> "LRRK2:MDKM34"
  strip_suffix <- function(nm) {
    sub("\\.\\d+$", "", nm)
  }
  
  # -- Collect data columns (everything except column 1 = concentration)-------
  # Use positional exclusion so plates with different concentration column names
  # are handled correctly.
  plate_data <- lapply(plates, function(p) {
    tbl       <- tables[[p]]
    data_cols <- colnames(tbl)[-1L]        # drop column 1 (concentration) by position
    bases     <- strip_suffix(data_cols)
    list(
      plate     = p,
      data_cols = data_cols,
      bases     = bases,
      values    = tbl[, data_cols, drop = FALSE]
    )
  })
  names(plate_data) <- plates
  
  # -- Determine the global order of base compound names ----------------------
  # Preserve the order in which compounds first appear (plate 1 first, then
  # any new compounds from plate 2, etc.).
  all_bases_ordered <- unique(unlist(lapply(plate_data, `[[`, "bases")))
  
  if (verbose) {
    cat(sprintf("Unique compounds : %d\n", length(all_bases_ordered)))
    for (p in plates) {
      n_cols <- length(plate_data[[p]]$data_cols)
      n_base <- length(unique(plate_data[[p]]$bases))
      cat(sprintf("  %-20s: %d column(s), %d unique compound(s)\n",
                  p, n_cols, n_base))
    }
    cat(strrep("-", 60), "\n")
  }
  
  # -- Build the merged data columns ------------------------------------------
  # For each base compound name, collect all replicate columns across plates
  # in plate order, then renumber them 1, 2, 3, 4, ...
  merged_cols <- list()  # will hold named numeric vectors (one per output column)
  
  for (base in all_bases_ordered) {
    
    all_reps <- list()  # collect replicate vectors in order
    
    for (p in plates) {
      pd      <- plate_data[[p]]
      matches <- which(pd$bases == base)
      if (length(matches) == 0L) next
      for (idx in matches) {
        all_reps <- c(all_reps, list(pd$values[[idx]]))
      }
    }
    
    # Assign output column names: first rep has no suffix, rest get .2, .3, ...
    n_reps <- length(all_reps)
    col_names <- character(n_reps)
    col_names[1L] <- base
    if (n_reps > 1L)
      col_names[-1L] <- paste0(base, ".", seq(2L, n_reps))
    
    for (i in seq_len(n_reps)) {
      merged_cols[[col_names[i]]] <- all_reps[[i]]
    }
  }
  
  # -- Assemble the merged data.frame -----------------------------------------
  conc_values <- tables[[1L]][[1L]]          # positional: always column 1
  merged_df   <- as.data.frame(
    c(list(conc_values), merged_cols),
    check.names = FALSE
  )
  colnames(merged_df)[1L] <- conc_col
  rownames(merged_df)     <- rownames(tables[[1L]])
  
  if (verbose) {
    cat(sprintf("Merged table     : %d rows x %d columns (%d data columns)\n",
                nrow(merged_df), ncol(merged_df), ncol(merged_df) - 1L))
    # Show per-compound replicate counts
    data_col_names <- colnames(merged_df)[-1L]
    base_counts    <- table(strip_suffix(data_col_names))
    cat("Replicates per compound:\n")
    for (nm in names(base_counts)) {
      cat(sprintf("  %-40s: %d replicate(s)\n", nm, base_counts[[nm]]))
    }
    cat(strrep("=", 60), "\n")
  }
  
  # -- Build the merged result entry -----------------------------------------
  # Mirror the structure of a single batch_ratio_analysis() plate entry.
  # Fields that are plate-specific or cannot be meaningfully merged are set
  # to NULL or to a descriptive character string.
  merged_entry <- list(
    data_file        = paste(vapply(plates, function(p)
      results[[p]]$data_file %||% p,
      character(1L)), collapse = " + "),
    info_sheet       = merged_name,
    sheet_number     = merged_name,
    function_version = results[[plates[1L]]]$function_version,
    control_0perc    = results[[plates[1L]]]$control_0perc,
    control_100perc  = results[[plates[1L]]]$control_100perc,
    selected_columns = results[[plates[1L]]]$selected_columns,
    merged_from      = plates,   # extra field - records provenance
    result           = list(
      modified_ratio_table  = merged_df,
      original_ratio_table  = NULL,   # no single "original" for a merge
      general_means         = NULL,   # recomputed downstream if needed
      interval_means        = NULL,   # recomputed downstream if needed
      control_info          = NULL,
      selected_columns_info = NULL
    )
  )
  
  # -- Assemble the output list -----------------------------------------------
  # Remove the individual merged plates and insert the merged entry.
  output <- results[setdiff(names(results), plates)]
  output[[merged_name]] <- merged_entry
  
  # Propagate the assay_source attribute so batch_drc_analysis() can still
  # auto-detect the assay type after merging (list subsetting drops attributes).
  src <- attr(results, "assay_source")
  if (!is.null(src)) attr(output, "assay_source") <- src
  
  # -- Excel report -----------------------------------------------------------
  if (generate_reports) {
    
    if (!requireNamespace("openxlsx", quietly = TRUE))
      warning("Package 'openxlsx' is required to save reports. Skipping.")
    else {
      
      # Write into drc_quality/ inside output_dir (or working directory if NULL)
      base_dir    <- if (!is.null(output_dir)) output_dir else getwd()
      quality_dir <- file.path(base_dir, "drc_quality")
      if (!dir.exists(quality_dir)) {
        dir.create(quality_dir, recursive = TRUE)
        if (verbose) message("Created output directory: ", quality_dir)
      }
      
      excel_path <- file.path(quality_dir,
                              paste0("merged_results_", merged_name, ".xlsx"))
      
      wb <- openxlsx::createWorkbook()
      
      header_style <- openxlsx::createStyle(
        fontColour     = "#FFFFFF",
        fgFill         = "#4F81BD",
        halign         = "center",
        textDecoration = "Bold",
        border         = "TopBottom",
        borderColour   = "#4F81BD"
      )
      
      # -- Sheet 1: Merged_Table ---------------------------------------------
      openxlsx::addWorksheet(wb, "Merged_Table")
      openxlsx::writeData(
        wb, "Merged_Table",
        cbind(RowNames = rownames(merged_df), merged_df),
        rowNames = FALSE
      )
      openxlsx::addStyle(wb, "Merged_Table", header_style,
                         rows = 1L, cols = seq_len(ncol(merged_df) + 1L))
      
      # -- Sheets 2..N: one per source plate (original modified_ratio_table) --
      for (p in plates) {
        orig_tbl <- tables[[p]]
        
        # Sheet names must be <= 31 chars and unique; truncate plate name if needed
        sheet_label <- paste0("Original_", p)
        if (nchar(sheet_label) > 31L)
          sheet_label <- paste0("Orig_", substr(p, 1L, 26L))
        
        openxlsx::addWorksheet(wb, sheet_label)
        openxlsx::writeData(
          wb, sheet_label,
          cbind(RowNames = rownames(orig_tbl), orig_tbl),
          rowNames = FALSE
        )
        openxlsx::addStyle(wb, sheet_label, header_style,
                           rows = 1L, cols = seq_len(ncol(orig_tbl) + 1L))
      }
      
      # -- Last sheet: Provenance --------------------------------------------
      openxlsx::addWorksheet(wb, "Provenance")
      
      # Per-plate rows
      prov_rows <- lapply(plates, function(p) {
        data_col_names <- colnames(tables[[p]])[-1L]
        base_counts    <- table(strip_suffix(data_col_names))
        data.frame(
          Plate        = p,
          Data_File    = results[[p]]$data_file    %||% "unknown",
          Info_Sheet   = results[[p]]$info_sheet   %||% p,
          N_Columns    = length(data_col_names),
          N_Compounds  = length(base_counts),
          stringsAsFactors = FALSE
        )
      })
      prov_df <- do.call(rbind, prov_rows)
      
      # Summary rows appended below
      summary_rows <- data.frame(
        Plate       = c("--- MERGED ---", "Merged name", "Merged date"),
        Data_File   = c("",
                        merged_name,
                        as.character(Sys.time())),
        Info_Sheet  = "",
        N_Columns   = c(NA_integer_,
                        ncol(merged_df) - 1L,
                        NA_integer_),
        N_Compounds = c(NA_integer_,
                        length(all_bases_ordered),
                        NA_integer_),
        stringsAsFactors = FALSE
      )
      
      openxlsx::writeData(wb, "Provenance",
                          rbind(prov_df, summary_rows),
                          rowNames = FALSE)
      openxlsx::addStyle(wb, "Provenance", header_style,
                         rows = 1L, cols = seq_len(ncol(prov_df)))
      
      openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
      
      if (verbose)
        message("Merged report saved: ", basename(excel_path))
    }
  }
  
  return(invisible(output))
}

# -- Null-coalescing operator (internal use) ----------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
