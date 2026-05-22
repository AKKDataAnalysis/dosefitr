#' Batch Save Dose-Response Curve Plots
#'
#' Generates and saves dose-response curve plots for all valid compounds across
#' multiple plates from batch DRC analysis results. Plots are saved to disk with
#' flexible directory organization and customizable aesthetics.
#'
#' This function scans all plates, identifies successfully fitted compounds,
#' and generates publication-quality plots using \code{plot_dose_response()}.
#'
#' @param batch_drc_results A list of DRC results. Can be either:
#'   \itemize{
#'     \item Output from \code{batch_drc_analysis()} (with \code{$drc_results})
#'     \item A direct list of plate-level DRC results
#'   }
#' @param output_dir Character. Directory where plots will be saved.
#'   Default is \code{"DRC_Plots"}.
#' @param organize_by Character. Directory organization strategy:
#'   \itemize{
#'     \item \code{"plate"}: plots grouped by plate (default)
#'     \item \code{"compound"}: plots grouped by compound
#'     \item \code{"flat"}: all plots in a single directory
#'   }
#' @param compounds_to_plot Optional character vector of compound names to include.
#'   If \code{NULL}, all compounds are plotted.
#' @param plates_to_plot Optional character vector of plate names to include.
#'   If \code{NULL}, all plates are processed.
#' @param format Character. File format for saved plots (e.g. \code{"png"}, \code{"pdf"}).
#' @param width Numeric. Plot width in inches.
#' @param height Numeric. Plot height in inches.
#' @param dpi Numeric. Plot resolution in dots per inch.
#' @param point_color Character. Color of data points in plots.
#' @param verbose Logical. If \code{TRUE}, prints progress and summary messages.
#' @param show_ic50_line Logical. Whether to display IC50 vertical line.
#' @param plot_title Logical. Whether to include plot titles.
#' @param point_size Numeric. Size of data points in plots.
#' @param y_limits Numeric vector of length 2 specifying the y-axis limits
#'   (e.g. \code{c(0, 100)}).  If \code{NULL} (default), each plot auto-scales
#'   to its own data range.  Passed directly to \code{plot_dose_response()}.
#' @param y_axis_title Character string for the y-axis label.  If \code{NULL}
#'   (default), auto-detected from the batch result: \code{"Cell Viability (\%)"}
#'   or \code{"Luminescence"} for viability assays, \code{"Normalized BRET ratio [\%]"}
#'   or \code{"BRET ratio"} for NanoBRET assays (depending on whether
#'   \code{normalize} was \code{TRUE} or \code{FALSE}).
#' @param save_panel Logical. If \code{TRUE} (default), saves one combined panel
#'   image per plate containing all individual compound plots assembled with
#'   \pkg{patchwork}.  Set to \code{FALSE} to skip panel generation.
#' @param subplot_title Character. Controls what text is used as the title of
#'   each compound sub-plot inside the panel. One of \code{"auto"} (default),
#'   \code{"full"} (e.g. \code{"KinaseA:Cpd1"}), \code{"compound"} (e.g.
#'   \code{"Cpd1"}), or \code{"construct"} (e.g. \code{"KinaseA"}).
#'   \code{"auto"} mirrors the same logic used for individual plots: shows only
#'   the compound name when all compounds in the batch share a single construct,
#'   and the full \code{Construct:Compound} string otherwise.
#' @param panel_ncol Integer. Number of columns in the panel grid (default \code{4}).
#' @param panel_width_per_col Numeric. Width in inches per panel column (default \code{4}).
#' @param panel_height_per_row Numeric. Height in inches per panel row (default \code{4}).
#' @param ... Additional arguments passed to \code{plot_dose_response()}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input structure and extracts plate-level results
#'   \item Identifies compounds with successful model fits
#'   \item Optionally filters by plate and/or compound
#'   \item Creates directory structure for output
#'   \item Generates plots using \code{plot_dose_response()}
#'   \item Saves plots to disk with safe filenames
#' }
#'
#' Compound and construct names are automatically parsed from input strings
#' (e.g. \code{"Construct | Compound"} or \code{"Construct:Compound"} formats).
#'
#' Filenames are sanitized to remove invalid filesystem characters.
#'
#' @return
#' Invisibly returns a list with summary information:
#' \itemize{
#'   \item \code{total} Number of compounds processed
#'   \item \code{successes} Number of successfully generated plots
#'   \item \code{failures} Number of failed plots
#'   \item \code{failed_compounds} Character vector of failed entries
#'   \item \code{error_messages} List of error messages
#'   \item \code{output_dir} Output directory path
#'   \item \code{organization} Directory structure used
#'   \item \code{point_color} Point color used in plots
#'   \item \code{timestamp} Time of execution
#' }
#'
#' @examples
#' \dontrun{
#' # Run batch DRC analysis first
#' results <- batch_drc_analysis(data_list)
#'
#' # Save all plots organized by plate
#' batch_save_all_drc_plots(results)
#'
#' # Save only selected compounds
#' batch_save_all_drc_plots(
#'   results,
#'   compounds_to_plot = c("DrugA", "DrugB")
#' )
#'
#' # Organize plots by compound
#' batch_save_all_drc_plots(
#'   results,
#'   organize_by = "compound",
#'   format = "pdf"
#' )
#' }
#'
#' @seealso
#' \code{\link{batch_drc_analysis}},
#' \code{\link{plot_dose_response}}
#'
#' @export

batch_save_all_drc_plots <- function(batch_drc_results,
                                     output_dir = "DRC_Plots",
                                     organize_by = "plate",
                                     compounds_to_plot = NULL,
                                     plates_to_plot = NULL,
                                     format = "png",
                                     width = 10,
                                     height = 10,
                                     dpi = 600,
                                     point_color = "black",
                                     verbose = TRUE,
                                     show_ic50_line = FALSE,
                                     plot_title = FALSE,
                                     point_size = 2,
                                     y_limits        = NULL,
                                     y_axis_title    = NULL,
                                     save_panel      = TRUE,
                                     panel_ncol      = 4L,
                                     panel_width_per_col  = 6,
                                     panel_height_per_row = 6,
                                     subplot_title = "auto",
                                     ...) {
  
  # ============================================================================
  # 1. VALIDATION AND SETUP
  # ============================================================================
  
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  if (save_panel && !requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' is required for panel assembly. ",
            "Install it with install.packages('patchwork'). ",
            "Individual plots will still be saved; panels will be skipped.")
    save_panel <- FALSE
  }
  
  # Helper function for safe filename generation
  safe_filename <- function(string) {
    if (is.null(string) || is.na(string)) return("unknown")
    # Remove any characters that could cause filesystem issues
    s <- gsub("[^[:alnum:]._-]", "_", string)
    s <- gsub("_+", "_", s)
    s <- gsub("^_|_$", "", s)
    if (nchar(s) == 0) return("unknown")
    return(s)
  }
  
  # Helper function to extract compound name properly
  extract_compound_name <- function(compound_string) {
    if (is.null(compound_string)) return("Unknown")
    
    # Remove replicate suffix if present (.1, .2, etc.)
    name <- gsub("\\.\\d+$", "", compound_string)
    
    # Handle "Construct | Compound" format
    if (grepl(" \\| ", name)) {
      parts <- strsplit(name, " \\| ")[[1]]
      return(trimws(parts[2]))
    }
    
    # Handle "Construct:Compound" format
    if (grepl(":", name)) {
      parts <- strsplit(name, ":")[[1]]
      return(trimws(parts[2]))
    }
    
    return(name)
  }
  
  # Helper function to extract construct name
  extract_construct_name <- function(compound_string) {
    if (is.null(compound_string)) return("Unknown")
    
    # Remove replicate suffix if present (.1, .2, etc.)
    name <- gsub("\\.\\d+$", "", compound_string)
    
    # Handle "Construct | Compound" format
    if (grepl(" \\| ", name)) {
      parts <- strsplit(name, " \\| ")[[1]]
      return(trimws(parts[1]))
    }
    
    # Handle "Construct:Compound" format
    if (grepl(":", name)) {
      parts <- strsplit(name, ":")[[1]]
      return(trimws(parts[1]))
    }
    
    return("Unknown")
  }
  
  # Extract drc_results if batch_drc_results is the wrapper object
  if (is.list(batch_drc_results)) {
    if ("drc_results" %in% names(batch_drc_results)) {
      if (verbose) message("Detected batch_drc_analysis wrapper. Extracting drc_results...")
      drc_results <- batch_drc_results$drc_results
    } else {
      drc_results <- batch_drc_results
    }
  } else {
    stop("batch_drc_results must be a list")
  }
  
  # Get plate names
  plate_names <- names(drc_results)
  if (is.null(plate_names) || length(plate_names) == 0) {
    stop("No plates found in drc_results")
  }
  
  # Filter plates if specified
  if (!is.null(plates_to_plot)) {
    plate_names <- intersect(plate_names, plates_to_plot)
    if (length(plate_names) == 0) {
      stop("No valid plates specified")
    }
  }
  
  if (verbose) {
    message("Found ", length(plate_names), " plates to process")
    message("Output directory: ", output_dir)
  }
  
  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  subplot_title <- match.arg(subplot_title, c("auto", "full", "compound", "construct"))
  
  # Helper: derive a display label from a raw "Construct:Compound" string
  # given an explicit mode (never "auto" — resolve that before calling).
  .label_for_mode <- function(compound_string, mode) {
    parts <- strsplit(compound_string, ":", fixed = TRUE)[[1L]]
    switch(mode,
           full      = compound_string,
           compound  = if (length(parts) >= 2L) parts[[2L]] else compound_string,
           construct = if (length(parts) >= 2L) parts[[1L]] else compound_string
    )
  }
  
  # Resolved at scan-time (after compounds_list is built):
  # single_construct_batch and effective_subplot_mode are set below.
  
  # ============================================================================
  # 2. SCAN FOR VALID COMPOUNDS ACROSS ALL PLATES
  # ============================================================================
  
  # Auto-detect y_axis_title from batch metadata when not supplied.
  if (is.null(y_axis_title)) {
    assay_src  <- batch_drc_results$metadata$assay_type
    normalized <- isTRUE(batch_drc_results$metadata$normalize)
    y_axis_title <- if (!is.null(assay_src) && assay_src == "viability") {
      if (normalized) "Cell Viability (%)" else "Luminescence"
    } else {
      if (normalized) "Normalized BRET ratio [%]" else "BRET ratio"
    }
  }
  
  if (verbose) message("\nScanning for valid compounds...")
  
  compounds_list <- list()
  
  for (plate_name in plate_names) {
    plate <- drc_results[[plate_name]]
    
    # Check if plate has drc_result
    if (is.null(plate$drc_result)) {
      if (verbose > 1) message("  Skipping ", plate_name, ": no drc_result")
      next
    }
    
    # Get detailed results
    detailed <- plate$drc_result$detailed_results
    if (is.null(detailed) || !is.list(detailed)) {
      if (verbose > 1) message("  Skipping ", plate_name, ": no detailed_results")
      next
    }
    
    # For each compound in the plate
    for (i in seq_along(detailed)) {
      result <- detailed[[i]]
      
      # Check if fit was successful
      if (!isTRUE(result$success)) next
      
      # Get compound and construct names
      compound_name <- extract_compound_name(result$compound)
      construct_name <- extract_construct_name(result$compound)
      
      # Filter by compound if specified
      if (!is.null(compounds_to_plot) && !compound_name %in% compounds_to_plot) next
      
      # Store compound info
      compounds_list <- append(compounds_list, list(list(
        plate = plate_name,
        construct = construct_name,
        compound = compound_name,
        compound_full = result$compound,
        index = i,
        results_obj = plate$drc_result  # Pass the full results object for plot_dose_response
      )))
    }
  }
  
  if (length(compounds_list) == 0) {
    stop("No valid compounds found to plot")
  }
  
  if (verbose) {
    message("Found ", length(compounds_list), " valid compounds")
    message("  - Plates: ", paste(unique(sapply(compounds_list, function(x) x$plate)), collapse = ", "))
    message("  - Compounds: ", length(unique(sapply(compounds_list, function(x) x$compound))))
  }
  
  # ── Auto-detect: resolve the best title mode from batch composition ─────────
  # Used by plot_title = TRUE (individual plots) and subplot_title = "auto" (panel).
  #   1 construct, N compounds  -> "compound"   (construct is constant, redundant)
  #   N constructs, 1 compound  -> "construct"  (compound is constant, redundant)
  #   N constructs, N compounds -> "full"        (need both to distinguish)
  #   1 construct,  1 compound  -> "compound"   (degenerate; either works)
  all_constructs      <- unique(sapply(compounds_list, function(x) x$construct))
  all_compounds       <- unique(sapply(compounds_list, function(x) x$compound))
  single_construct_batch <- length(all_constructs) == 1L
  single_compound_batch  <- length(all_compounds)  == 1L
  
  auto_mode <- if (single_construct_batch) {
    "compound"    # one construct — compound name is the distinguishing label
  } else if (single_compound_batch) {
    "construct"   # one compound  — construct name is the distinguishing label
  } else {
    "full"        # many of both  — need Construct:Compound
  }
  
  # Resolve effective panel mode now that we know the batch composition.
  effective_subplot_mode <- if (subplot_title == "auto") auto_mode else subplot_title
  
  # ============================================================================
  # 3. CREATE DIRECTORY STRUCTURE
  # ============================================================================
  
  if (organize_by == "plate") {
    # Create subfolders for each plate
    for (plate_name in unique(sapply(compounds_list, function(x) x$plate))) {
      plate_dir <- file.path(output_dir, safe_filename(plate_name))
      if (!dir.exists(plate_dir)) {
        dir.create(plate_dir, recursive = TRUE)
      }
    }
  } else if (organize_by == "compound") {
    # Create subfolders for each compound
    for (compound_name in unique(sapply(compounds_list, function(x) x$compound))) {
      compound_dir <- file.path(output_dir, safe_filename(compound_name))
      if (!dir.exists(compound_dir)) {
        dir.create(compound_dir, recursive = TRUE)
      }
    }
  }
  
  # ============================================================================
  # 4. GENERATE ALL PLOTS
  # ============================================================================
  
  if (verbose) message("\nGenerating plots...")
  
  total <- length(compounds_list)
  
  # Use a local environment for mutable state inside tryCatch handlers.
  # This avoids <<- leaking variables into the caller's global environment.
  state <- new.env(parent = emptyenv())
  state$successes     <- 0L
  state$failures      <- 0L
  state$failed_list   <- character()
  state$error_messages <- list()
  state$plate_plots   <- list()   # plate_name -> list of ggplot objects for panel
  
  # Progress bar
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }
  
  for (i in seq_along(compounds_list)) {
    info <- compounds_list[[i]]
    
    # Filenames always use the full Construct_Compound name to guarantee
    # uniqueness when the same compound is tested against multiple constructs.
    base_name <- safe_filename(info$compound_full)
    
    if (organize_by == "plate") {
      output_path <- file.path(output_dir, safe_filename(info$plate),
                               paste0(base_name, ".", format))
    } else if (organize_by == "compound") {
      output_path <- file.path(output_dir, safe_filename(info$compound),
                               paste0(safe_filename(info$plate), "_", base_name, ".", format))
    } else {
      output_path <- file.path(output_dir,
                               paste0(safe_filename(info$plate), "_", base_name, ".", format))
    }
    
    # Create directory if it doesn't exist
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
    
    tryCatch({
      # Individual plot — respects the user's plot_title setting.
      # When plot_title = TRUE, auto-detect: use compound-only if the whole
      # batch has a single construct, otherwise use the full Construct:Compound.
      indiv_title <- if (isTRUE(plot_title)) {
        .label_for_mode(info$compound_full, auto_mode)
      } else {
        plot_title   # FALSE or a custom character string — pass through as-is
      }
      p_single <- plot_dose_response(
        results = info$results_obj,
        compound_index = info$index,
        save_plot = NULL,
        plot_width = width,
        plot_height = height,
        plot_dpi = dpi,
        point_color = point_color,
        show_ic50_line = show_ic50_line,
        verbose = FALSE,
        plot_title = indiv_title,
        point_size = point_size,
        y_limits     = y_limits,
        y_axis_title = y_axis_title,
        ...
      )
      
      # Save individual file
      ggplot2::ggsave(output_path, plot = p_single,
                      width = width, height = height,
                      dpi = dpi, bg = "white")
      
      state$successes <- state$successes + 1L
      
      # Panel version always shows a title using effective_subplot_mode,
      # regardless of what the user chose for the individual files.
      if (save_panel) {
        panel_label <- .label_for_mode(info$compound_full, effective_subplot_mode)
        p_panel <- if (isTRUE(plot_title) && identical(indiv_title, panel_label)) {
          p_single   # already has the right title — reuse without a second call
        } else {
          plot_dose_response(
            results = info$results_obj,
            compound_index = info$index,
            save_plot = NULL,
            plot_width = width,
            plot_height = height,
            plot_dpi = dpi,
            point_color = point_color,
            show_ic50_line = show_ic50_line,
            verbose = FALSE,
            plot_title = panel_label,
            point_size = point_size,
            y_limits     = y_limits,
            y_axis_title = y_axis_title,
            ...
          )
        }
        state$plate_plots[[info$plate]] <- c(state$plate_plots[[info$plate]], list(p_panel))
      }
    }, error = function(e) {
      state$failures      <- state$failures + 1L
      state$failed_list   <- c(state$failed_list, paste(info$plate, info$compound, sep = "/"))
      state$error_messages <- c(state$error_messages,
                                list(paste(info$plate, info$compound, ":", e$message)))
    })
    
    if (verbose) setTxtProgressBar(pb, i)
  }
  
  # Extract state back to plain variables for the rest of the function
  successes     <- state$successes
  failures      <- state$failures
  failed_list   <- state$failed_list
  error_messages <- state$error_messages
  plate_plots   <- state$plate_plots
  
  if (verbose) close(pb)
  
  # ============================================================================
  # 5. PANEL ASSEMBLY (one combined image per plate)
  # ============================================================================
  
  panel_files <- character()
  
  if (save_panel && length(plate_plots) > 0L) {
    if (verbose) message("\nAssembling panels...")
    
    for (plate_name in names(plate_plots)) {
      plot_list <- plate_plots[[plate_name]]
      n_plots   <- length(plot_list)
      if (n_plots == 0L) next
      
      n_cols_panel <- min(panel_ncol, n_plots)
      n_rows_panel <- ceiling(n_plots / n_cols_panel)
      panel_w      <- n_cols_panel * panel_width_per_col
      panel_h      <- n_rows_panel * panel_height_per_row + 0.6  # +0.6 for title
      
      combined <- patchwork::wrap_plots(plot_list, ncol = n_cols_panel)
      
      # Save panel next to the individual plots
      panel_filename <- paste0(safe_filename(plate_name), "_panel.", format)
      panel_path <- switch(organize_by,
                           plate   = file.path(output_dir, safe_filename(plate_name), panel_filename),
                           compound = file.path(output_dir, panel_filename),
                           flat     = file.path(output_dir, panel_filename)
      )
      
      tryCatch({
        ggplot2::ggsave(panel_path, combined,
                        width = panel_w, height = panel_h,
                        dpi = dpi, bg = "white")
        panel_files <- c(panel_files, panel_path)
        if (verbose) message(sprintf("  Panel saved: %s  (%d plots, %.0fx%.0f in)",
                                     basename(panel_path), n_plots,
                                     panel_w, panel_h))
      }, error = function(e) {
        warning(sprintf("Failed to save panel for plate '%s': %s",
                        plate_name, e$message))
      })
    }
  }
  
  # ============================================================================
  # 6. SUMMARY AND RETURN
  # ============================================================================
  
  if (verbose) {
    message("\n")
    message("========================================")
    message("PLOT GENERATION COMPLETE")
    message("========================================")
    message("Total compounds: ", total)
    message("Successful: ", successes)
    message("Failed: ", failures)
    message("Point color: ", point_color)
    message("Output directory: ", normalizePath(output_dir))
    
    if (failures > 0) {
      message("\nFailed compounds:")
      for (f in failed_list) {
        message("  - ", f)
      }
      
      if (verbose > 1) {
        message("\nError details:")
        for (err in error_messages) {
          message("  - ", err)
        }
      }
    }
    
    # Show directory structure
    message("\nDirectory structure:")
    if (organize_by == "plate") {
      for (plate in unique(sapply(compounds_list, function(x) x$plate))) {
        plate_files <- list.files(file.path(output_dir, safe_filename(plate)), pattern = paste0("\\.", format, "$"))
        message("  ", plate, "/ : ", length(plate_files), " files")
      }
    } else if (organize_by == "compound") {
      for (compound in unique(sapply(compounds_list, function(x) x$compound))) {
        compound_files <- list.files(file.path(output_dir, safe_filename(compound)), pattern = paste0("\\.", format, "$"))
        message("  ", compound, "/ : ", length(compound_files), " files")
      }
    } else {
      all_files <- list.files(output_dir, pattern = paste0("\\.", format, "$"))
      message("  Flat structure: ", length(all_files), " files in root")
    }
  }
  
  # ============================================================================
  # 7. RETURN INVISIBLE SUMMARY
  # ============================================================================
  
  invisible(list(
    total = total,
    successes = successes,
    failures = failures,
    failed_compounds = failed_list,
    error_messages = error_messages,
    output_dir = output_dir,
    organization = organize_by,
    point_color = point_color,
    panel_files = panel_files,
    timestamp = Sys.time()
  ))
}




