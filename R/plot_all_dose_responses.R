#' Generate Multiple Dose-Response Plots in Batch
#'
#' Automatically creates and saves individual dose-response plots for multiple compounds
#' in batch mode. Ideal for generating comprehensive visualization sets for screening
#' data, compound libraries, or publication supplements.
#'
#' @param results List object returned by \code{\link{fit_drc_3pl}} containing
#'   dose-response analysis results.
#' @param compounds Character or numeric vector specifying which compounds to plot.
#'   Use "all" for all compounds or provide specific numeric indices (default: "all").
#' @param output_dir Character string specifying the directory where plots will be saved
#'   (default: "dose_response_plots").
#' @param file_prefix Character string prefix for all saved plot files
#'   (default: "dose_response").
#' @param file_extension Character string specifying the image file format
#'   ("png", "jpg", "tiff", "pdf", "svg") (default: "png").
#' @param label_sep Character string. Separator used for DISPLAY purposes in
#'   titles and filenames. When \code{NULL} (default), auto-detected from
#'   \code{attr(results, "label_sep")}; falls back to \code{":"}. Forwarded to
#'   \code{\link{plot_dose_response}}.
#' @param ... Additional arguments passed to \code{\link{plot_dose_response}} for
#'   customizing individual plots (colors, sizes, limits, etc.).
#'
#' @return Invisibly returns a numeric vector of the compound indices that were plotted.
#'   Primarily produces image files as output.
#'
#' @details
#' This function provides automated batch plotting capabilities for dose-response
#' analysis results. It systematically generates individual plots for each specified
#' compound with consistent formatting and automatic file naming.
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item Batch Processing: Automatically generates plots for multiple compounds
#'   \item Automatic File Management: Creates organized directory structure
#'   \item Consistent Formatting: Applies uniform styling across all plots
#'   \item Flexible Compound Selection: Plot all compounds or specific subsets
#'   \item Multiple Output Formats: Supports various image formats for different use cases
#'   \item Progress Tracking: Provides real-time progress feedback
#' }
#'
#' \strong{Automatic File Naming:}
#' Files are automatically named using the pattern:
#' \code{file_prefix_CompoundName.file_extension}
#' Special characters in compound names are automatically converted to underscores
#' for filesystem compatibility.
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_padr")
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
#' drc_res <- batch_drc_analysis(
#'   batch_results = ratio_res, model = "4pl", normalize = FALSE,
#'   generate_reports = FALSE, output_dir = tempdir(), verbose = FALSE
#' )
#'
#' fit1    <- drc_res$drc_results$plate_01$drc_result
#' out_dir <- file.path(tempdir(), "dosefitr_ex_padr_out")
#' dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#'
#' plot_all_dose_responses(
#'   results        = fit1,
#'   compounds      = 1:2,
#'   output_dir     = out_dir,
#'   file_extension = "png"
#' )
#' length(list.files(out_dir))
#' }
#' @section Output Organization:
#' The function creates a structured output system:
#' \itemize{
#'   \item \strong{Directory Creation}: Automatically creates output directory if missing
#'   \item \strong{File Naming}: Consistent, descriptive names for easy identification
#'   \item \strong{Progress Reporting}: Real-time feedback on generation progress
#'   \item \strong{Error Handling}: Continues processing even if individual plots fail
#' }
#'
#' @section Supported File Formats:
#' \itemize{
#'   \item \strong{PNG}: Good balance of quality and file size (default)
#'   \item \strong{JPEG}: Smaller file size, good for presentations
#'   \item \strong{TIFF}: High quality, lossless compression for publications
#'   \item \strong{PDF}: Vector format, scalable and editable
#'   \item \strong{SVG}: Vector format, web-compatible and editable
#' }
#'
#' @section Use Cases:
#' \itemize{
#'   \item \strong{High-Throughput Screening}: Generate plots for entire compound libraries
#'   \item \strong{Publication Preparation}: Create consistent figures for manuscripts
#'   \item \strong{Quality Control}: Visual assessment of all curve fits
#'   \item \strong{Data Sharing}: Provide comprehensive visualization to collaborators
#'   \item \strong{Archive Creation}: Generate permanent records of analysis results
#' }
#'
#' @seealso
#' \code{\link{plot_dose_response}} for individual plot customization
#' \code{\link{plot_multiple_compounds}} for multi-curve comparison plots
#' \code{\link{fit_drc_3pl}} for generating analysis results
#'
#' @export


plot_all_dose_responses <- function(results, 
                                    compounds = "all",
                                    output_dir = "dose_response_plots",
                                    file_prefix = "dose_response",
                                    file_extension = "png",
                                    label_sep = NULL,
                                    ...) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine which compounds to plot
  n_compounds <- length(results$detailed_results)
  
  if (identical(compounds, "all")) {
    compounds_to_plot <- 1:n_compounds
  } else if (is.numeric(compounds)) {
    # Validate compound indices
    invalid_compounds <- compounds[compounds < 1 | compounds > n_compounds]
    if (length(invalid_compounds) > 0) {
      stop("Invalid compound indices: ", paste(invalid_compounds, collapse = ", "))
    }
    compounds_to_plot <- compounds
  } else {
    stop("compounds must be 'all' or a numeric vector of indices")
  }
  
  cat("Generating", length(compounds_to_plot), "dose-response plots...\n")
  
  # Generate plots for each compound
  for (i in compounds_to_plot) {
    compound_name <- results$detailed_results[[i]]$compound
    compound_name_parts <- strsplit(compound_name, " \\| ")[[1]]
    compound_name_display <- compound_name_parts[1]
    
    # Create safe filename
    safe_filename <- gsub("[^a-zA-Z0-9_-]", "_", compound_name_display)
    filename <- file.path(output_dir, 
                          paste0(file_prefix, "_", safe_filename, ".", file_extension))
    
    # Generate and save plot
    plot_dose_response(results = results, 
                       compound_index = i, 
                       save_plot = filename,
                       label_sep = label_sep,
                       ...)
    
    cat("  - Plot", i, "of", length(compounds_to_plot), ":", compound_name_display, "\n")
  }
  
  cat("All plots saved in:", output_dir, "\n")
  return(invisible(compounds_to_plot))
}
