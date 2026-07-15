#' Plot Dose-Response Curves with ggplot2
#'
#' Generates publication-quality dose-response plots using ggplot2 from analysis results.
#' Features professional styling, error bars, fitted curves, IC50 visualization, and
#' multiple export formats. Optimized for scientific publications.
#'
#' @param results List object returned by \code{\link{fit_drc_3pl}} containing
#'   dose-response analysis results.
#' @param compound_index Numeric index specifying which compound to plot (default: 1).
#' @param x_limits Numeric vector of length 2 specifying the x-axis limits
#'   in log10 molar units.  \code{NULL} (default) uses the data range.
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
#' @param label_sep Character string. Separator used for DISPLAY purposes in
#'   titles, filenames, and metadata. When \code{NULL} (default), auto-detected
#'   from \code{attr(results, "label_sep")}; falls back to \code{":"} if the
#'   attribute is absent. This only affects what the user sees — the internal
#'   data separator used for parsing compound names is always read from the
#'   attribute and is never changed. For example, \code{label_sep = "/"} renders
#'   \code{"Kinase/Compound"} in the title while the data still stores
#'   \code{"Kinase:Compound"} internally.
#' @param axis_line_width Numeric. Line width of the manually drawn x/y axis
#'   lines (default: 0.8).
#' @param axis_vjust Numeric or NULL. Vertical justification (\code{vjust}) of
#'   the axis titles. \code{NULL} (default) leaves the ggplot2 default unchanged.
#' @param tick_length Numeric or NULL. Axis tick length in centimetres.
#'   \code{NULL} (default) preserves the \code{theme_minimal()} default.
#' @param error_linewidth Numeric. Line width of the error bars (default: 0.8).
#' @param point_alpha Numeric between 0 and 1. Opacity of the data points
#'   (default: 1, fully opaque). Does not affect error-bar opacity.
#' @param legend_spacing Numeric or NULL. Spacing (in points) between legend
#'   items, applied via \code{legend.spacing}. \code{NULL} (default) leaves the
#'   theme default unchanged.
#' @param aspect_ratio Numeric or NULL. Panel aspect ratio
#'   (\code{aspect.ratio} in \code{theme()}). \code{NULL} (default) leaves it unset.
#' @param byrow Logical. Whether legend keys are filled by row (\code{TRUE}) or
#'   by column (\code{FALSE}, default). Applied via \code{theme(legend.byrow=)}
#'   (the modern ggplot2 mechanism; the old \code{guide_legend(byrow=)} argument
#'   was removed in ggplot2 4.0).
#' @param axis_line_color Character string. Colour of the axis lines and tick
#'   marks (default: \code{"black"}).
#' @param transparent_background Logical. If \code{TRUE}, the plot and panel
#'   backgrounds are set to fully transparent (\code{element_rect(fill = NA)})
#'   instead of white. Default: \code{FALSE}.
#' @param panel_border Logical. If \code{TRUE}, draws a full rectangular border
#'   around the plot panel (using \code{axis_line_color}). Default: \code{FALSE}.
#' @param plot_title_size Numeric or NULL. Font size for the plot title.
#'   \code{NULL} (default) uses \code{axis_label_size + 2}.
#' @param legend_text_size Numeric or NULL. Font size for legend text.
#'   \code{NULL} (default) leaves the ggplot2 default unchanged.
#' @param legend_title_size Numeric or NULL. Font size for the legend title.
#'   \code{NULL} (default) leaves the ggplot2 default unchanged.
#' @param curve_alpha Numeric between 0 and 1. Opacity of the fitted
#'   dose-response curve (default: 1, fully opaque).
#' @param plot_margin Margin or NULL. Plot margin applied via
#'   \code{theme(plot.margin = )}. Accepts a \code{ggplot2::margin()} object.
#'   \code{NULL} (default) uses the built-in margin
#'   (t = 12, r = 8, b = 8, l = 8 pt).
#' @param axis_title_color Character string. Colour of the axis titles
#'   (default: \code{"black"}).
#' @param axis_text_color Character string. Colour of the axis text (tick labels)
#'   (default: \code{"black"}).
#' @param ic50_linetype Character or integer. Line type for the IC50 vertical
#'   reference line (default: \code{"dashed"}). See \code{?par} for valid values.
#' @param ic50_linewidth Numeric. Line width of the IC50 vertical reference line
#'   (default: 0.8).
#' @param ic50_line_alpha Numeric between 0 and 1. Opacity of the IC50 vertical
#'   reference line (default: 1, fully opaque).
#' @param error_alpha Numeric between 0 and 1. Opacity of the error bars
#'   (default: 1, fully opaque).
#' @param point_shape Integer or character. Shape of the data points
#'   (default: 16, filled circle). See \code{?points} for valid values.
#' @param grid_color Character string. Colour of the major grid lines when
#'   \code{show_grid = TRUE} (default: \code{"grey90"}).
#' @param grid_minor_color Character string. Colour of the minor grid lines when
#'   \code{show_grid = TRUE} (default: \code{"grey95"}).
#' @param grid_linewidth Numeric. Line width of the major grid lines
#'   (default: 0.5).
#' @param axis_expand Numeric vector of length 2. Expansion constants for the
#'   continuous axis scales, passed to \code{scale_*_continuous(expand = )}.
#'   Default \code{c(0, 0)} removes all padding between data and axis edge.
#' @param title_hjust Numeric. Horizontal justification of the plot title
#'   (default: 0.5, centred). Use 0 for left-aligned, 1 for right-aligned.
#' @param point_size_scale Numeric. Multiplier applied to \code{point_size} for
#'   the actual \code{geom_point} size (default: 2). The effective point size is
#'   \code{point_size * point_size_scale}.
#' @param legend_annotation_scale Numeric. Scale factor applied to
#'   \code{axis_text_size} for the in-plot legend text annotations
#'   (default: 0.3). The effective annotation size is
#'   \code{axis_text_size * legend_annotation_scale}.
#' @param verbose Logical indicating whether to show verbose messages (default: FALSE).
#' @param excluded_point_color Character. Colour of data points that were excluded
#'   from fitting by hook-effect detection. Default: `"red"`. Only visible when
#'   the result contains `data_full` with `excluded = TRUE` rows (i.e. when
#'   `hook_effect` was enabled in `batch_drc_analysis()`).
#' @param excluded_point_shape Integer. ggplot2 shape code for excluded points.
#'   Default: `4` (\eqn{\times} symbol).
#' @param excluded_point_size Numeric or `NULL`. Size of excluded points.
#'   `NULL` (default) inherits `point_size`.
#' @param excluded_point_alpha Numeric between 0 and 1. Opacity of excluded points.
#'   Default: `1`.
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
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_pdr")
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
#' fit1 <- drc_res$drc_results$plate_01$drc_result
#' p    <- plot_dose_response(fit1, compound_index = 1, verbose = FALSE)
#' inherits(p, "ggplot")
#' }
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
                               x_limits = NULL,
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
                               label_sep = NULL,
                               axis_line_width = 0.8,
                               axis_vjust = NULL,
                               tick_length = NULL,
                               error_linewidth = 0.8,
                               point_alpha = 1,
                               legend_spacing = NULL,
                               aspect_ratio = NULL,
                               byrow = FALSE,
                               axis_line_color = "black",
                               transparent_background = FALSE,
                               panel_border = FALSE,
                               plot_title_size = NULL,
                               legend_text_size = NULL,
                               legend_title_size = NULL,
                               curve_alpha = 1,
                               plot_margin = NULL,
                               axis_title_color = "black",
                               axis_text_color = "black",
                               ic50_linetype = "dashed",
                               ic50_linewidth = 0.8,
                               ic50_line_alpha = 1,
                               error_alpha = 1,
                               point_shape = 16,
                               grid_color = "grey90",
                               grid_minor_color = "grey95",
                               grid_linewidth = 0.5,
                               axis_expand = c(0, 0),
                               title_hjust = 0.5,
                               point_size_scale = 2,
                               legend_annotation_scale = 0.3,
                               verbose = FALSE,
                               excluded_point_color = "red",
                               excluded_point_shape = 4,
                               excluded_point_size  = NULL,
                               excluded_point_alpha = 1) {
  
  # Null-coalescing operator
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

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
  
# Resolve label_sep: the separator used for DISPLAY purposes (titles, labels,
# filenames). 
if (is.null(label_sep)) {
  label_sep <- attr(results, "label_sep")
  if (is.null(label_sep) || !is.character(label_sep) ||
      length(label_sep) != 1L || is.na(label_sep) || nchar(label_sep) == 0L) {
    label_sep <- ":"
  }
}

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

  # Compute excluded-point summary (hook-effect points: plotted but not fitted)
  excluded_summary <- NULL
  if (!is.null(result$data_full) && is.data.frame(result$data_full) &&
      "excluded" %in% names(result$data_full) &&
      any(result$data_full$excluded, na.rm = TRUE)) {
    excluded_raw <- result$data_full[result$data_full$excluded == TRUE, , drop = FALSE]
    if (nrow(excluded_raw) > 0L && all(c("log_inhibitor", "response") %in% names(excluded_raw))) {
      excluded_summary <- calculate_summary_stats(excluded_raw)
    }
  }
  
  if (nrow(summary_data) == 0) {
    stop("No valid summary data available for plotting")
  }
  
  # Extract compound name (remove plate info if present)
  compound_name_display <- strsplit(result$compound, " \\| ")[[1]][1]

  # Create a display version of the compound name where the internal
  # separator is replaced with the user-facing label_sep. 
  .find_data_sep <- function(name) {
    # Try the attribute first (most reliable when available)
    attr_sep <- attr(results, "label_sep")
    if (!is.null(attr_sep) && nchar(attr_sep) == 1L &&
        grepl(attr_sep, name, fixed = TRUE)) return(attr_sep)
    # Fall back to common separators in priority order
    for (s in c(":", "/", " | ")) {
      if (grepl(s, name, fixed = TRUE)) return(s)
    }
    return(NULL)  # no separator found
  }
  data_sep_found <- .find_data_sep(compound_name_display)
  compound_name_label <- if (!is.null(data_sep_found) && data_sep_found != label_sep) {
    parts <- strsplit(compound_name_display, data_sep_found, fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      paste(parts[1], paste(parts[-1], collapse = data_sep_found), sep = label_sep)
    } else {
      compound_name_display
    }
  } else {
    compound_name_display
  }
  
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
      final_title <- compound_name_label
    } else {
      final_title <- paste(compound_name_label, "(Model failed)")
    }
  }
  
  # Create base plot with professional styling
  p <- ggplot2::ggplot() +
    ggplot2::labs(
      x = plot_config$x_lab,
      y = plot_config$y_lab,
      title = final_title
    ) +
    ggplot2::scale_y_continuous(expand = axis_expand) +
    ggplot2::scale_x_continuous(expand = axis_expand) +
    ggplot2::coord_cartesian(
      xlim = if (!is.null(x_limits) && length(x_limits) == 2L) x_limits else NULL,
      ylim = if (!is.null(y_limits) && length(y_limits) == 2L) y_limits else NULL,
      clip = "on") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = axis_label_size, face = "bold",
                                           color = axis_title_color, vjust = axis_vjust),
      axis.text = ggplot2::element_text(size = axis_text_size, color = axis_text_color),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = axis_line_color, linewidth = axis_line_width),
      plot.title = ggplot2::element_text(size = if (!is.null(plot_title_size)) plot_title_size else axis_label_size + 2,
                                         face = "bold", hjust = title_hjust, color = axis_title_color),
      legend.position = ifelse(show_legend, "right", "none"),
      legend.byrow = byrow,
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
      panel.grid.major = ggplot2::element_line(color = ifelse(show_grid, grid_color, "white"), linewidth = grid_linewidth),
      panel.grid.minor = ggplot2::element_line(color = ifelse(show_grid, grid_minor_color, "white")),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 12, r = 8, b = 8, l = 8, unit = "pt")
    )

  # Optional theme tweaks (only applied when explicitly set, so the defaults
  # leave the appearance unchanged).
  if (!is.null(tick_length)) {
    p <- p + ggplot2::theme(axis.ticks.length = ggplot2::unit(tick_length, "cm"))
  }
  if (!is.null(aspect_ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect_ratio)
  }
  if (!is.null(legend_spacing)) {
    p <- p + ggplot2::theme(legend.key.spacing.y = ggplot2::unit(legend_spacing, "pt"))
  }
  if (!is.null(plot_margin)) {
    p <- p + ggplot2::theme(plot.margin = plot_margin)
  }
  if (transparent_background) {
    p <- p + ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, color = NA),
      plot.background = ggplot2::element_rect(fill = NA, color = NA)
    )
  }
  if (panel_border) {
    p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = axis_line_color, fill = NA, linewidth = 0.5)
    )
  }
  
  # Draw axis lines manually so they stop exactly at the data limits.
  x_range_data <- range(summary_data$log_inhibitor, na.rm = TRUE)
  x_lo <- x_range_data[1] - diff(x_range_data) * 0.02   # mirrors expand mult
  x_hi <- x_range_data[2] + diff(x_range_data) * 0.02
  # Resolve y_limits: if NULL (auto), derive from the data so segment endpoints
  # are always finite scalars.
  y_seg_limits <- if (!is.null(y_limits) && length(y_limits) == 2L) {
    y_limits
  } else {
    y_range_data <- range(summary_data$mean_response, na.rm = TRUE)
    y_pad <- diff(y_range_data) * 0.05
    c(y_range_data[1] - y_pad, y_range_data[2] + y_pad)
  }
  axis_segs <- data.frame(
    x    = c(x_lo,              x_lo),
    xend = c(x_hi,              x_lo),
    y    = c(y_seg_limits[1],   y_seg_limits[1]),
    yend = c(y_seg_limits[1],   y_seg_limits[2])
  )
  p <- p +
    ggplot2::geom_segment(
      data = axis_segs,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      colour = axis_line_color, linewidth = axis_line_width,
      inherit.aes = FALSE)
  
  # Add experimental data points
  p <- p +
    ggplot2::geom_point(
      data = summary_data,
      ggplot2::aes(x = log_inhibitor, y = mean_response),
      color = point_color,
      size = point_size * point_size_scale,
      shape = point_shape,
      alpha = point_alpha
    )
  
  # Add error bars for replicates
  valid_mask <- !is.na(summary_data$sd_response) &
    is.finite(summary_data$sd_response) &
    summary_data$n_replicates > 1 &
    summary_data$sd_response > 1e-10
  
  if (any(valid_mask)) {
    valid_data <- summary_data[valid_mask, ]
    # Clip error bar whiskers at the y-axis limits so they never
    # extend below/above the visible plot area.
    y_lo_clip <- y_seg_limits[1]
    y_hi_clip <- y_seg_limits[2]
    valid_data$ymin_clipped <- pmax(valid_data$mean_response - valid_data$sd_response, y_lo_clip)
    valid_data$ymax_clipped <- pmin(valid_data$mean_response + valid_data$sd_response, y_hi_clip)
    p <- p +
      ggplot2::geom_errorbar(
        data = valid_data,
        ggplot2::aes(
          x = log_inhibitor,
          ymin = ymin_clipped,
          ymax = ymax_clipped
        ),
        width = error_bar_width * 10,
        color = point_color,
        linewidth = error_linewidth,
        alpha = error_alpha
      )
  }
  
  # Add excluded (hook-effect) points as red x symbols
  if (!is.null(excluded_summary) && nrow(excluded_summary) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data = excluded_summary,
        ggplot2::aes(x = log_inhibitor, y = mean_response),
        color = excluded_point_color,
        shape = excluded_point_shape,
        size  = (excluded_point_size %||% point_size) * point_size_scale,
        alpha = excluded_point_alpha
      )
    # Error bars for excluded points (only when SD is meaningful)
    excl_valid_mask <- !is.na(excluded_summary$sd_response) &
      is.finite(excluded_summary$sd_response) &
      excluded_summary$n_replicates > 1L &
      excluded_summary$sd_response > 1e-10
    if (any(excl_valid_mask)) {
      excl_valid_data <- excluded_summary[excl_valid_mask, ]
      excl_valid_data$ymin_clipped <- pmax(
        excl_valid_data$mean_response - excl_valid_data$sd_response, y_seg_limits[1])
      excl_valid_data$ymax_clipped <- pmin(
        excl_valid_data$mean_response + excl_valid_data$sd_response, y_seg_limits[2])
      p <- p +
        ggplot2::geom_errorbar(
          data = excl_valid_data,
          ggplot2::aes(
            x    = log_inhibitor,
            ymin = ymin_clipped,
            ymax = ymax_clipped
          ),
          width     = error_bar_width * 10,
          color     = excluded_point_color,
          linewidth = error_linewidth,
          alpha     = error_alpha
        )
    }
  }

  # Add fitted curve if model was successful
  if (!is.null(curve_data)) {
    p <- p +
      ggplot2::geom_line(
        data = curve_data,
        ggplot2::aes(x = log_inhibitor, y = response),
        color = line_color,
        linewidth = line_width / 2,
        alpha = curve_alpha
      )
  }
  
  # Add IC50 line only if valid IC50 exists
  if (show_ic50_line && is.finite(log_ic50) &&
      !(!is.null(ic50_excluded) && ic50_excluded)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = log_ic50,
        linetype = ic50_linetype,
        color = ic50_line_color,
        linewidth = ic50_linewidth,
        alpha = ic50_line_alpha
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
          size = axis_text_size * legend_annotation_scale,
          color = axis_text_color
        )
    }
  }
  
  # Save plot if requested
  if (!is.null(save_plot)) {
    if (is.character(save_plot)) {
      filename <- save_plot
    } else if (is.logical(save_plot) && save_plot) {
      safe_name <- gsub("[^a-zA-Z0-9._-]", "_", compound_name_label)
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
    compound_name = compound_name_label,
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

