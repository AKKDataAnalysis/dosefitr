#' Plot Dose-Response Curves with ggplot2
#'
#' Generates publication-quality dose-response plots using ggplot2 from analysis results.
#' Features professional styling, error bars, fitted curves, IC50 visualization, and
#' multiple export formats. Optimized for scientific publications.
#'
#' @param results List object returned by \code{\link{fit_drc_3pl}} containing
#'   dose-response analysis results.
#' @param compound_index Numeric index specifying which compound to plot (default: 1).
#' @param y_limits Numeric vector of length 2 specifying y-axis limits (default: c(0, 150)).
#' @param point_color Color for data points (default: "black").
#' @param line_color Color for fitted curve (default: "black").
#' @param ic50_line_color Color for IC50 vertical line (default: "gray").
#' @param point_size Size multiplier for data points (default: 2).
#' @param line_width Line width for fitted curve (default: 2).
#' @param error_bar_width Width of error bar ends (default: 0.01).
#' @param show_ic50_line Logical indicating whether to show vertical IC50 line (default: FALSE).
#' @param show_legend Logical indicating whether to show parameter legend (default: FALSE).
#' @param show_grid Logical indicating whether to show background grid (default: FALSE).
#' @param save_plot Defines whether to save the plot: \code{NULL} (do not save, default),
#'   \code{TRUE} (automatically saves as PNG with default name), or a file path with extension
#'   (\code{.png}, \code{.pdf}, \code{.jpeg}, \code{.tiff}, \code{.svg}, \code{.eps}) to save in a specific format.
#' @param plot_width Plot width in inches for saved plots (default: 10).
#' @param plot_height Plot height in inches for saved plots (default: 10).
#' @param plot_dpi Resolution for saved raster images (default: 600).
#' @param axis_label_size Font size for axis labels (default: 20).
#' @param axis_text_size Font size for axis numbers (default: 18).
#' @param x_axis_title Custom x-axis title. If NULL, uses default expression.
#' @param y_axis_title Custom y-axis title. If NULL, uses default based on normalization.
#' @param enforce_bottom_threshold Logical indicating whether bottom threshold enforcement
#'   was used in analysis (default: NULL, auto-detected from results).
#' @param bottom_threshold Numeric value for bottom threshold (default: 60).
#' @param verbose Logical indicating whether to show verbose messages (default: FALSE).
#' @param plot_title Controls the plot title. \code{FALSE} (default) = no title;
#'   \code{TRUE} = automatic title (construct + compound name);
#'   character = custom title text.
#'
#' @importFrom ggplot2 aes
#'
#' @return Returns a ggplot object with comprehensive metadata stored as attributes:
#' \itemize{
#'   \item \code{compound_name}: Name of the plotted compound
#'   \item \code{compound_index}: Index of the plotted compound
#'   \item \code{model_success}: Whether model fitting was successful
#'   \item \code{summary_data}: Data frame with summarized plotting data
#'   \item \code{plot_config}: Configuration settings used for plotting
#'   \item \code{y_limits_used}: Y-axis limits actually used
#'   \item \code{data_points}: Number of data points plotted
#'   \item \code{concentration_levels}: Number of concentration levels
#'   \item \code{file_saved}: Path to saved file if plot was saved
#'   \item \code{plot_dimensions}: Dimensions of the plot (width, height, dpi)
#'   \item \code{timestamp}: Time when plot was generated
#'   \item \code{ic50_excluded}: Whether IC50 was excluded due to threshold
#'   \item \code{log_ic50}: Log IC50 value if available
#' }
#'
#' @details
#' This function creates professional-quality dose-response plots using ggplot2,
#' providing superior visual quality and customization compared to base R graphics.
#' The function automatically handles data validation, curve fitting visualization,
#' and professional styling suitable for scientific publications.
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Professional Styling}: Clean, publication-ready appearance with bold axis labels
#'   \item \strong{Smart IC50 Handling}: Only shows IC50 line when valid value exists
#'   \item \strong{Error Bar Management}: Automatic handling of replicate data
#'   \item \strong{Threshold Awareness}: Detects and handles IC50 exclusion due to bottom thresholds
#'   \item \strong{Self-Contained}: No external package loading required
#'   \item \strong{Comprehensive Metadata}: Detailed information about the plot and data
#' }
#'
#' \strong{Plot Elements:}
#' \itemize{
#'   \item Data points with optional error bars (standard deviation)
#'   \item Fitted dose-response curve (when model converged)
#'   \item Vertical IC50 line (only when valid IC50 exists)
#'   \item Left-aligned parameter legend with IC50 and R2 values
#'   \item Professional axis formatting with customizable titles
#'   \item Optional background grid for better readability
#' }
#'
#' \strong{Supported Export Formats:}
#' \itemize{
#'   \item PNG (high-resolution, recommended for publications)
#'   \item JPEG (good for presentations)
#'   \item TIFF (lossless compression)
#'   \item PDF (vector format, scalable)
#'   \item SVG (vector format, editable)
#' }
#'
#' @examples
#' \dontrun{
#' # Perform dose-response analysis first
#' analysis_results <- fit_drc_3pl(my_data, normalize = TRUE)
#'
#' # Basic plot for first compound
#' p <- plot_dose_response(analysis_results)
#' print(p)
#'
#' # Customized plot with specific styling
#' p <- plot_dose_response(
#'   results = analysis_results,
#'   compound_index = 2,
#'   point_color = "blue",
#'   line_color = "darkred",
#'   show_grid = TRUE,
#'   y_limits = c(0, 200)
#' )
#' print(p)
#'
#' # Save plot automatically with compound name
#' p <- plot_dose_response(
#'   results = analysis_results,
#'   compound_index = 1,
#'   save_plot = TRUE  # Saves as "dose_response_CompoundName.png"
#' )
#'
#' # Save plot to specific file with custom dimensions
#' p <- plot_dose_response(
#'   results = analysis_results,
#'   save_plot = "my_plot.pdf",
#'   plot_width = 8,
#'   plot_height = 6
#' )
#'
#' # Access plot metadata
#' p <- plot_dose_response(analysis_results, compound_index = 3)
#' meta <- attr(p, "metadata")
#' print(meta$compound_name)
#' print(meta$data_points)
#' print(meta$log_ic50)
#' }
#'
#' @section Plot Customization:
#' The function provides extensive customization options while maintaining professional quality:
#' \itemize{
#'   \item \strong{Visual Elements}: Customize colors, sizes, and styling of all plot components
#'   \item \strong{Layout Control}: Adjust dimensions, resolution, and aspect ratio
#'   \item \strong{Content Toggle}: Show/hide legend, grid, IC50 line, and error bars
#'   \item \strong{Text Formatting}: Control font sizes and axis labels
#'   \item \strong{Export Quality}: Multiple formats with quality and resolution control
#' }
#'
#' @section Automatic Features:
#' \itemize{
#'   \item Automatic input validation and error handling
#'   \item Smart detection of threshold-based IC50 exclusion
#'   \item Automatic directory creation for saved plots
#'   \item Intelligent error bar display (only when meaningful data exists)
#'   \item Graceful handling of failed model fits with clear indication
#'   \item Professional axis formatting with appropriate scientific notation
#'   \item Compound name extraction and clean display
#' }
#'
#' @seealso
#' \code{\link{fit_drc_3pl}} for generating analysis results
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#' \code{\link[ggplot2]{ggsave}} for plot export functionality
#'
#' @export
#'
#' @references
#' For visualization best practices in scientific publications:
#' \itemize{
#'   \item Nature Scientific Figures Guidelines
#'   \item ggplot2: Elegant Graphics for Data Analysis (Springer)
#'   \item R Graphics Cookbook (O'Reilly)
#' }



plot_dose_response <- function(results, compound_index = 1, y_limits = c(0, 150),
                               point_color = "black", line_color = "black",
                               ic50_line_color = "gray", point_size = 2,
                               line_width = 2, error_bar_width = 0.01,
                               show_ic50_line = TRUE, show_legend = FALSE,
                               show_grid = FALSE, save_plot = NULL,
                               plot_width = 10, plot_height = 10, plot_dpi = 600,
                               axis_label_size = 20, axis_text_size = 18,
                               x_axis_title = NULL, y_axis_title = NULL,
                               plot_title = TRUE,
                               enforce_bottom_threshold = NULL, bottom_threshold = 60,
                               verbose = FALSE) {

  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Please install it with: install.packages('ggplot2')")
  }

  # Input validation
  validate_inputs <- function(results, compound_index) {
    if (missing(results)) {
      stop("Argument 'results' is required")
    }

    if (!is.list(results) || !"detailed_results" %in% names(results)) {
      stop("Invalid 'results' object. Must be output from fit function")
    }

    if (length(results$detailed_results) == 0) {
      stop("No compounds found in results object")
    }

    if (compound_index < 1 || compound_index > length(results$detailed_results)) {
      stop("Compound index ", compound_index, " out of range. Must be between 1 and ",
           length(results$detailed_results))
    }
  }

  validate_inputs(results, compound_index)

  # Extract compound data
  result <- results$detailed_results[[compound_index]]

  if (!is.list(result) || !"data" %in% names(result)) {
    stop("Invalid result structure for compound ", compound_index)
  }

  # Use threshold settings from results if not explicitly provided
  if (is.null(enforce_bottom_threshold)) {
    enforce_bottom_threshold <- if (!is.null(results$threshold_settings) &&
                                    !is.null(results$threshold_settings$enforce_bottom_threshold)) {
      results$threshold_settings$enforce_bottom_threshold
    } else {
      FALSE
    }
  }

  if (is.null(bottom_threshold)) {
    bottom_threshold <- if (!is.null(results$threshold_settings) &&
                            !is.null(results$threshold_settings$bottom_threshold)) {
      results$threshold_settings$bottom_threshold
    } else {
      60
    }
  }

  # Clean and validate data
  clean_data <- stats::na.omit(result$data)
  required_cols <- c("log_inhibitor", "response")
  if (nrow(clean_data) == 0 || !all(required_cols %in% names(clean_data))) {
    stop("No valid data points available or missing required columns")
  }

  # Calculate summary statistics (mean +- SD per concentration)
  calculate_summary_stats <- function(data) {
    summary_data <- do.call(rbind, lapply(split(data, data$log_inhibitor), function(sub_df) {
      n_rep <- nrow(sub_df)
      data.frame(
        log_inhibitor = unique(sub_df$log_inhibitor),
        mean_response = mean(sub_df$response, na.rm = TRUE),
        sd_response = if(n_rep > 1) sd(sub_df$response, na.rm = TRUE) else 0,
        n_replicates = nrow(sub_df)
      )
    }))
    stats::na.omit(summary_data)
  }

  summary_data <- calculate_summary_stats(clean_data)

  if (nrow(summary_data) == 0) {
    stop("No valid summary data available for plotting")
  }

  # Extract compound name (remove plate info if present)
  compound_name_display <- strsplit(result$compound, " \\| ")[[1]][1]

  # Check if IC50 was excluded due to threshold
  ic50_excluded <- FALSE
  if (!is.null(enforce_bottom_threshold) && isTRUE(enforce_bottom_threshold)) {
    comp_name <- strsplit(result$compound, " \\| ")[[1]][1]

    if (!is.null(results$summary_table) && "Compound" %in% names(results$summary_table)) {
      summary_row <- results$summary_table[results$summary_table$Compound == comp_name, ]

      if (nrow(summary_row) > 0 && !is.null(summary_row$IC50)) {
        ic50_value <- summary_row$IC50
        ic50_excluded <- is.na(ic50_value) ||
          as.character(ic50_value) == "NA" ||
          as.character(ic50_value) == "<NA>" ||
          as.character(ic50_value) == ""
      }
    }
  }

  # Setup plot configuration
  setup_plot_config <- function() {
    x_lab <- if (!is.null(x_axis_title)) {
      x_axis_title
    } else {
      expression(paste("Log"[10], " Concentration [M]"))
    }

    y_lab <- if (!is.null(y_axis_title)) {
      y_axis_title
    } else {
      ifelse(results$normalized, "Normalized BRET ratio [%]", "BRET ratio")
    }

    list(
      x_lab = x_lab,
      y_lab = y_lab,
      point_color = point_color,
      line_color = line_color,
      point_size = point_size,
      line_width = line_width,
      error_bar_width = error_bar_width,
      axis_label_size = axis_label_size,
      axis_text_size = axis_text_size
    )
  }

  plot_config <- setup_plot_config()

  # Generate fitted curve data
  generate_fitted_curve <- function(model) {
    if (is.null(model)) return(NULL)

    x_range <- range(summary_data$log_inhibitor, na.rm = TRUE)
    if (!all(is.finite(x_range))) return(NULL)

    x_seq <- seq(x_range[1], x_range[2], length.out = 300)

    predictions <- tryCatch({
      predict(model, newdata = data.frame(log_inhibitor = x_seq))
    }, error = function(e) NULL)

    if (!is.null(predictions)) {
      data.frame(log_inhibitor = x_seq, response = predictions)
    }
  }

  # Get IC50 value for vertical line
  get_ic50_value <- function(model) {
    if (is.null(model)) return(NA)

    tryCatch({
      coefs <- stats::coef(model)
      if ("LogIC50" %in% names(coefs)) coefs["LogIC50"] else NA
    }, error = function(e) NA)
  }

  # Create legend text
  create_legend_content <- function(model = NULL) {
    if (!show_legend) return(NULL)

    if (!is.null(model) && isTRUE(result$success)) {
      log_ic50 <- get_ic50_value(model)
      ic50_value <- if (is.finite(log_ic50)) 10^log_ic50 else NA
      r_squared <- round(result$goodness_of_fit$R_squared, 3)

      legend_text <- c()

      if (!is.null(ic50_excluded) && ic50_excluded) {
        legend_text <- c(legend_text, "LogIC50 = NA (threshold)")
        legend_text <- c(legend_text, "IC50 = NA (threshold)")
      } else if (is.finite(log_ic50)) {
        legend_text <- c(legend_text, paste("LogIC50 =", round(log_ic50, 3)))
        legend_text <- c(legend_text, paste("IC50 =", sprintf("%.2e", ic50_value)))
      } else {
        legend_text <- c(legend_text, "LogIC50 = NA")
        legend_text <- c(legend_text, "IC50 = NA")
      }

      legend_text <- c(legend_text, paste("R2 =", r_squared))

      return(legend_text)

    } else {
      return("Model did not converge")
    }
  }

  # Model status and data preparation
  model_success <- !is.null(result$model) && isTRUE(result$success)
  curve_data <- if (model_success) generate_fitted_curve(result$model) else NULL
  log_ic50 <- if (model_success) get_ic50_value(result$model) else NA
  legend_content <- create_legend_content(if (model_success) result$model else NULL)

  # ============================================================================
  # Determine the title based on argument plot_title
  # ============================================================================

  final_title <- NULL

  if (is.character(plot_title)) {
    final_title <- plot_title
  } else if (isTRUE(plot_title)) {
    if (model_success) {
      final_title <- compound_name_display
    } else {
      final_title <- paste(compound_name_display, "(Model failed)")
    }
  }

  # Create base plot with professional styling
  p <- ggplot2::ggplot() +
    ggplot2::labs(
      x = plot_config$x_lab,
      y = plot_config$y_lab,
      title = final_title
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = axis_label_size, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = axis_text_size, color = "black"),
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      plot.title = ggplot2::element_text(size = axis_label_size + 2, face = "bold", hjust = 0.5, color = "black"),
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(color = ifelse(show_grid, "grey90", "white")),
      panel.grid.minor = ggplot2::element_line(color = ifelse(show_grid, "grey95", "white")),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  # Add experimental data points
  p <- p +
    ggplot2::geom_point(
      data = summary_data,
      ggplot2::aes(x = log_inhibitor, y = mean_response),
      color = point_color,
      size = point_size * 2
    )

  # Add error bars for replicates
  valid_mask <- !is.na(summary_data$sd_response) &
    is.finite(summary_data$sd_response) &
    summary_data$n_replicates > 1 &
    summary_data$sd_response > 1e-10

  if (any(valid_mask)) {
    valid_data <- summary_data[valid_mask, ]
    p <- p +
      ggplot2::geom_errorbar(
        data = valid_data,
        ggplot2::aes(
          x = log_inhibitor,
          ymin = mean_response - sd_response,
          ymax = mean_response + sd_response
        ),
        width = error_bar_width * 10,
        color = point_color,
        linewidth = 0.8
      )
  }

  # Add fitted curve if model was successful
  if (!is.null(curve_data)) {
    p <- p +
      ggplot2::geom_line(
        data = curve_data,
        ggplot2::aes(x = log_inhibitor, y = response),
        color = line_color,
        linewidth = line_width / 2
      )
  }

  # Add IC50 line only if valid IC50 exists
  if (show_ic50_line && is.finite(log_ic50) &&
      !(!is.null(ic50_excluded) && ic50_excluded)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = log_ic50,
        linetype = "dashed",
        color = ic50_line_color,
        linewidth = 0.8
      )
  }

  # Add legend as text annotation (left-aligned)
  if (show_legend && !is.null(legend_content)) {
    x_range <- range(summary_data$log_inhibitor, na.rm = TRUE)
    y_range <- y_limits

    x_pos <- x_range[1] + diff(x_range) * 0.02
    y_pos <- y_range[1] + diff(y_range) * 0.1

    for (i in seq_along(legend_content)) {
      p <- p +
        ggplot2::annotate(
          "text",
          x = x_pos,
          y = y_pos - diff(y_range) * 0.05 * (i - 1),
          label = legend_content[i],
          hjust = 0,  # Left alignment
          vjust = 0,
          size = axis_text_size * 0.3,
          color = "black"
        )
    }
  }

  # Save plot if requested
  if (!is.null(save_plot)) {
    if (is.character(save_plot)) {
      filename <- save_plot
    } else if (is.logical(save_plot) && save_plot) {
      safe_name <- gsub("[^a-zA-Z0-9._-]", "_", compound_name_display)
      filename <- paste0("dose_response_", safe_name, ".png")
    } else {
      stop("save_plot must be either a file path or TRUE for auto-naming")
    }

    plot_dir <- dirname(filename)
    if (plot_dir != "." && !dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }

    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = plot_width,
      height = plot_height,
      dpi = plot_dpi,
      bg = "white"
    )

    if (verbose) {
      message("Plot saved as: ", normalizePath(filename))
    }
  }

  # Return plot with comprehensive metadata
  metadata <- list(
    compound_name = compound_name_display,
    compound_index = compound_index,
    model_success = model_success,
    summary_data = summary_data,
    plot_config = plot_config,
    y_limits_used = y_limits,
    data_points = nrow(clean_data),
    concentration_levels = nrow(summary_data),
    file_saved = if (!is.null(save_plot)) filename else NULL,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    timestamp = Sys.time(),
    ic50_excluded = ic50_excluded,
    log_ic50 = if (model_success) log_ic50 else NA,
    title_mode = if (is.character(plot_title)) "custom" else if (isTRUE(plot_title)) "auto" else "none"
  )

  attr(p, "metadata") <- metadata
  return(p)
}
