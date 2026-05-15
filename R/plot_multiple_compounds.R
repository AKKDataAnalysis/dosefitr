#' Plot Multiple Dose-Response Curves
#'
#' The `plot_multiple_compounds()` function generates a consolidated plot of
#' fitted dose-response curves for multiple compounds, allowing visual
#' comparison between different responses. It provides extensive customization
#' options for colors, shapes, titles, legends, gridlines, and file export.
#'
#' @param results The list returned by [`batch_drc_analysis()`].
#'   Must contain a `drc_results` element (one entry per plate), each with
#'   a `drc_result` sub-list that holds `detailed_results`.
#' @param compound_indices Numeric vector specifying which compounds to include in the plot
#'   If `NULL` (default), all available compounds in `results$detailed_results`
#'   are included.
#' @param target_compound A character string specifying a compound name, target name, or
#'   combination in "target:compound" format. Supports partial matching (e.g., "GPCR" will
#'   match all compounds with "GPCR" in the target name).
#' @param position An integer specifying a single compound position to plot (alternative
#'   to `target_compound` or `compound_indices`).
#' @param x_limits Numeric vector of length 2 specifying the x-axis limits.
#'   If \code{NULL} (default), limits are calculated automatically from the
#'   data with a 5\% margin.  The unit of the supplied values is controlled
#'   by \code{x_limits_scale}.
#' @param x_limits_scale Character string declaring the unit of \code{x_limits}.
#'   One of \code{"log10"} (default - values are already log\eqn{_{10}} molar,
#'   e.g. \code{c(-9, -5)}), \code{"molar"} (raw molar, e.g.
#'   \code{c(1e-9, 25e-6)}), \code{"uM"} (micromolar, e.g. \code{c(0.001, 25)}),
#'   or \code{"nM"} (nanomolar, e.g. \code{c(1, 25000)}).  Ignored when
#'   \code{x_limits = NULL}.
#' @param y_limits Numeric vector of length 2 specifying the y-axis limits
#' @param point_shapes Numeric vector of point shapes for different compounds
#' @param colors Either a logical (`TRUE` for automatic colors), a character vector of
#'   color names/hex codes, or a palette name (see Details for available palettes). If
#'   `TRUE` or NULL, uses the palette specified by `color_palette` or default hue palette.
#' @param color_palette A character string specifying a named palette to use when `colors`
#'   is `TRUE` or NULL. See Details for comprehensive list of available palettes including
#'   scientific journals, colorblind-friendly options, and gradient schemes.
#' @param shape_by_compound Logical; if TRUE, uses different point shapes to distinguish
#'   compounds in addition to colors. Useful for colorblind accessibility or grayscale
#'   printing.
#' @param legend_position Position of the legend: one of `"right"`, `"left"`,
#'   `"top"`, `"bottom"`, or `"none"`. Default: `"right"`.
#' @param show_grid Logical indicating whether to show background grid lines
#' @param transparent_background Logical.  If \code{TRUE}, the plot, panel,
#'   and legend backgrounds are set to fully transparent (\code{element_rect(fill = NA)})
#'   instead of white. Default: \code{FALSE}.
#' @param show_legend Logical indicating whether to display the legend
#' @param save_plot Defines whether to save the plot: \code{NULL} (do not save, default),
#'   \code{TRUE} (automatically saves as PNG with default name), or a file path with extension
#'   (\code{.png}, \code{.pdf}, \code{.jpeg}, \code{.tiff}, \code{.svg}, \code{.eps}) to save in a specific format.
#' @param plot_width,plot_height Plot dimensions (in inches) when saving.
#' @param plot_dpi Resolution (in DPI) for saved plots. Default: `600`.
#' @param axis_text_size Numeric value for axis tick label font size
#' @param axis_title_size Numeric; size of axis titles. Default is 14.
#' @param axis_title_color Character string specifying color for axis titles. Default is "black".
#' @param axis_line_color Character string.  Colour of the axis lines and tick
#'   marks.  Default: \code{"black"}.
#' @param show_border Logical.  If \code{TRUE}, draws a full rectangular border
#'   around the plot panel (using \code{axis_line_color}).  Default: \code{FALSE}.
#' @param x_axis_title Character string for a custom x-axis label.  If \code{NULL}
#'   (default), the label is \code{Log10 Concentration [M]}.
#' @param y_axis_title Character string for custom y-axis title. If NULL, automatically
#'   determined based on whether normalized data is detected.
#' @param axis_text_color Character string specifying color for axis text. Default is "black".
#' @param show_error_bars Logical indicating whether to display error bars around data points
#' @param error_bar_width Numeric value controlling the width of error bars
#' @param curve_linewidth Numeric.  Line width of the fitted dose-response curves.
#'   Default: \code{1}.
#' @param curve_alpha Numeric between 0 and 1.  Opacity of the fitted curves.
#'   \code{1} = fully opaque, \code{0} = invisible.  Default: \code{0.7}.
#' @param show_ic50_lines Logical.  If \code{TRUE}, draws a vertical dashed line
#'   at each compound's IC50 position on the x-axis, coloured to match the
#'   corresponding curve.  Lines are only drawn for compounds with a valid
#'   (non-NA, in-range) IC50.  Default: \code{FALSE}.
#' @param plot_title Custom plot title. If `NULL`, a smart title is generated
#'   automatically based on the compound names.
#' @param plot_title_size Numeric.  Font size for the plot title.  If \code{NULL}
#'   (default), uses \code{axis_title_size + 2}.
#' @param legend_text_size Numeric value for legend text font size
#' @param legend_title_size Numeric value for legend title font size
#' @param legend_ncol Numeric value specifying number of columns in legend
#' @param legend_label_wrap Maximum character width before legend labels
#'   automatically wrap to new lines. Default: `25`.
#' @param legend_title Title for the legend (displayed above symbols).
#' @param verbose Logical; if TRUE (default), prints informative messages about processing
#'   steps, compound matches, and color assignments.
#' @param point_size Numeric value for the point size
#' @param plate Character. Plate name to plot. When \code{results} is the
#'   output of \code{batch_drc_analysis()}, the function detects this
#'   automatically and extracts the correct \code{drc_result}. If
#'   \code{plate = NULL} (default) and a batch result is supplied, the
#'   first plate is used and a message lists all available plates. Pass an
#'   explicit plate name (e.g. \code{plate = "plate_01"}) to select a
#'   specific plate. Has no effect when \code{results} is already a plain
#'   \code{drc_result} object.
#'
#'
#'@importFrom ggplot2 aes
#'
#'
#' @return A ggplot2 object with the generated plot, containing metadata attributes:
#'   \item{selected_compounds}{Character vector of selected compound names}
#'   \item{smart_legend_names}{Automatically generated legend labels}
#'   \item{n_compounds}{Number of compounds plotted}
#'   \item{match_type}{How compounds were selected (position, indices, pattern match)}
#'   \item{colors}{Vector of colors used in the plot}
#'   \item{color_palette_used}{Name of the palette used}
#'   \item{available_palettes}{List of all available palette options}
#'   \item{...}{Additional metadata about plot settings}
#'
#' @details
#' This function overlays fitted dose-response curves (based on nonlinear models)
#' together with empirical mean +- SD values for each concentration, allowing
#' direct visual comparison across multiple compounds or experimental conditions.
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item \strong{Flexible compound selection}: Select by indices, pattern matching on target/compound names, or position
#'   \item \strong{Extensive color palette options}: Over 50 predefined palettes including scientific journals, colorblind-friendly schemes, and gradients
#'   \item \strong{Smart point selection}: Automatically chooses optimal point shapes based on compound count
#'   \item \strong{Adaptive sizing}: Point sizes and legend elements adjust based on number of compounds
#'   \item \strong{Intelligent text wrapping}: Automatically wraps long compound names in legend
#'   \item \strong{Intelligent legend labels}: Automatically simplifies labels when targets or compounds are consistent
#'   \item \strong{Professional styling}: Clean, publication-ready appearance with customizable elements
#'   \item \strong{Self-contained}: No external package loading required
#' }
#'
#' \strong{Color Palette Options:}
#'
#' The function provides extensive color palette options suitable for scientific publications.
#' Palettes can be specified via the `color_palette` parameter or directly via `colors`.
#'
#' \emph{Base Palettes:}
#' \itemize{
#'   \item \code{"hue"}, \code{"ggplot2"}, \code{"default"} - Default ggplot2 hue palette
#' }
#'
#' \emph{ColorBrewer Palettes:}
#' \itemize{
#'   \item \strong{Qualitative} (good for distinct categories): \code{"set1"}, \code{"set2"}, \code{"set3"}, \code{"dark2"}, \code{"paired"}, \code{"accent"}, \code{"pastel1"}, \code{"pastel2"}
#'   \item \strong{Sequential} (for ordered data): \code{"blues"}, \code{"reds"}, \code{"greens"}, \code{"purples"}, \code{"oranges"}, \code{"greys"}
#'   \item \strong{Diverging} (for deviation from midpoint): \code{"spectral"}, \code{"rdylbu"}, \code{"rdylgn"}, \code{"piyg"}, \code{"prgn"}
#' }
#'
#' \emph{Viridis Palettes} (perceptually uniform, colorblind-friendly):
#' \code{"viridis"}, \code{"magma"}, \code{"inferno"}, \code{"plasma"}
#'
#' \emph{Scientific Journal Palettes:}
#' \itemize{
#'   \item \code{"nature"} - Classic Nature journal colors
#'   \item \code{"science"} - Vibrant Science magazine palette
#'   \item \code{"cell"} - Cell journal distinctive colors
#'   \item \code{"plos"} - PLOS ONE accessible palette
#'   \item \code{"elife"} - eLife journal color scheme
#' }
#'
#' \emph{Publisher Palettes:}
#' \itemize{
#'   \item \code{"bmc"} (BioMed Central), \code{"frontiers"}, \code{"wiley"}, \code{"elsevier"},
#'   \item \code{"oxford"} (Oxford University Press), \code{"springer"} (Springer Nature),
#'   \item \code{"acs"} (American Chemical Society), \code{"rsc"} (Royal Society of Chemistry)
#' }
#'
#' \emph{Colorblind-Friendly Palettes:}
#' \itemize{
#'   \item \code{"okabe_ito"} - Okabe & Ito palette (best for colorblindness)
#'   \item \code{"colorblind"} - Generic colorblind-friendly palette
#'   \item \code{"cud"} - Color Universal Design palette
#'   \item \code{"tol"} - Paul Tol's colorblind-friendly scheme
#' }
#'
#' \emph{Gradient Palettes} (for ordered/continuous data):
#' \itemize{
#'   \item \code{"blue_red"}, \code{"green_red"}, \code{"purple_orange"},
#'   \item \code{"cool_warm"}, \code{"blue_yellow"}
#' }
#'
#' \emph{Corporate Palettes:}
#' \code{"ibm"}, \code{"google"}, \code{"microsoft"}, \code{"twitter"}
#'
#' \emph{Classic R Palettes:}
#' \code{"rainbow"}, \code{"heat"}, \code{"terrain"}, \code{"topo"}, \code{"cm"}
#'
#' \strong{Automatic Adjustments:}
#' \itemize{
#'   \item \strong{Point shapes}: Uses most distinguishable shapes first (16, 17, 15, 18, 8, etc.), recycles intelligently
#'   \item \strong{Point size}: Larger points for few compounds (3.5), smaller for many (2.5)
#'   \item \strong{Legend text}: Smaller text for many compounds (9pt), larger for few (11pt)
#'   \item \strong{Legend columns}: Single column for <=10 compounds, two columns for >10 compounds
#'   \item \strong{Legend labels}: Automatically simplified - shows only compound names when target is consistent, only targets when compound is consistent, or full "target:compound" format otherwise
#'   \item \strong{X-axis limits}: Automatically calculated from data with 5% margin
#'   \item \strong{Color assignment}: Intelligently selects appropriate palettes based on number of compounds and usage context
#' }
#'
#' \strong{Color Selection Logic:}
#' The function determines colors in the following priority order:
#' \enumerate{
#'   \item If `colors` is a vector of color names/hex codes, those colors are used directly
#'   \item If `color_palette` is specified, colors are generated from that named palette
#'   \item If `colors = TRUE`, colors are generated from the default hue palette
#'   \item Otherwise, the default hue palette is used
#' }
#'
#' @note
#' Required packages: ggplot2, scales, grDevices.
#' Optional packages for extended palettes: RColorBrewer, viridisLite.
#'
#' When using RColorBrewer palettes, the package will automatically interpolate colors
#' if more colors are needed than the palette's maximum. When using viridis palettes,
#' the function falls back to colorRampPalette if viridisLite is not installed.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic plot for all compounds
#' p1 <- plot_multiple_compounds(results)
#' print(p1)
#'
#' # Example 2: Select specific compounds by index
#' p2 <- plot_multiple_compounds(results, compound_indices = c(1, 3, 5))
#'
#' # Example 3: Enable automatic coloring and save the plot
#' p3 <- plot_multiple_compounds(
#'   results,
#'   compound_indices = 1:4,
#'   colors = TRUE,
#'   save_plot = "plots/multi_colored_curves.png"
#' )
#'
#' # Example 4: Customize legend labels and point shapes
#' p4 <- plot_multiple_compounds(
#'   results,
#'   compound_indices = 1:3,
#'   legend_labels = c("Compound A", "Compound B", "Compound C"),
#'   point_shapes = c(15, 17, 19),
#'   legend_title = "Treatments"
#' )
#'
#' # Example 5: Disable error bars and use black-and-white mode
#' p5 <- plot_multiple_compounds(
#'   results,
#'   colors = FALSE,
#'   show_error_bars = FALSE,
#'   show_grid = TRUE,
#'   plot_title = "Curves without error bars"
#' )
#'
#' # Example 6: Place legend below the plot in multiple columns
#' p6 <- plot_multiple_compounds(
#'   results,
#'   legend_position = "bottom",
#'   legend_ncol = 3,
#'   colors = TRUE
#' )
#'
#' # Example 7: Plot without any legend
#' p7 <- plot_multiple_compounds(results, show_legend = FALSE)
#'
#' # Example 8: Fine-tune fonts, limits, and title
#' p8 <- plot_multiple_compounds(
#'   results,
#'   compound_indices = 1:2,
#'   y_limits = c(0, 120),
#'   axis_label_size = 16,
#'   axis_text_size = 13,
#'   plot_title = "Comparison of Two Compounds"
#' )
#'
#'
#'
#' # Extract metadata for reproducibility
#' meta <- attr(p, "metadata")
#' cat("Plotted", meta$n_compounds, "compounds\n")
#' cat("Point shapes:", meta$point_shapes, "\n")
#' cat("X-axis range:", round(meta$x_limits, 2), "\n")
#'
#' # Save styling information
#' styling_info <- data.frame(
#'   Compound = meta$compound_names,
#'   PointShape = meta$point_shapes,
#'   Color = meta$colors
#' )
#' write.csv(styling_info, "plot_styling.csv", row.names = FALSE)
#' }
#'
#' @seealso
#' \code{\link{fit_drc_3pl}} for generating input data
#' \code{\link{fit_drc_4pl}} for generating input data
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#'
#' @export



plot_multiple_compounds <- function(results,
                                    compound_indices = NULL,
                                    target_compound = NULL,
                                    position = NULL,
                                    x_limits = NULL,
                                    x_limits_scale = "log10",
                                    y_limits = c(0, 150),
                                    colors = NULL,
                                    color_palette = NULL,
                                    point_shapes = NULL,
                                    show_error_bars = TRUE,
                                    show_ic50_lines = FALSE,
                                    legend_position = "right",
                                    show_legend = TRUE,
                                    shape_by_compound = FALSE,
                                    show_grid = FALSE,
                                    show_border = FALSE,
                                    transparent_background = FALSE,
                                    save_plot = NULL,
                                    plot_width = 10,
                                    plot_height = 8,
                                    plot_dpi = 600,
                                    plot_title = NULL,
                                    plot_title_size = NULL,
                                    legend_title = "",
                                    legend_text_size = NULL,
                                    legend_title_size = 11,
                                    x_axis_title = NULL,
                                    y_axis_title = NULL,
                                    verbose = TRUE,
                                    axis_text_color = "black",
                                    axis_text_size = 12,
                                    axis_title_color = "black",
                                    axis_title_size = 14,
                                    axis_line_color = "black",
                                    curve_linewidth = 1,
                                    curve_alpha = 0.7,
                                    error_bar_width = 0.05,
                                    legend_label_wrap = 25,
                                    legend_ncol = NULL,
                                    point_size = NULL,
                                    plate = NULL) {

  # ============================================================================
  # 0. PLATE EXTRACTION (when batch_drc_analysis result is supplied)
  # ============================================================================
  # Detect whether `results` is a batch_drc_analysis output by checking for
  # the $drc_results slot. If so, extract the requested plate's drc_result.
  # When `plate` is NULL, the first plate is used automatically.
  if (!is.null(results$drc_results)) {
    batch_obj <- results$drc_results
    available_plates <- names(batch_obj)
    if (is.null(plate)) {
      plate <- available_plates[1]
      if (verbose) {
        message("batch_drc_analysis result detected. Using plate: '", plate,
                "'. Available plates: ", paste(available_plates, collapse = ", "),
                ". Use the `plate` argument to select a different plate.")
      }
    } else {
      if (!plate %in% available_plates) {
        stop("Plate '", plate, "' not found in batch results. ",
             "Available plates: ", paste(available_plates, collapse = ", "))
      }
    }
    results <- batch_obj[[plate]]$drc_result
    if (is.null(results)) {
      stop("No drc_result found for plate '", plate, "'.")
    }
  } else if (!is.null(plate)) {
    # results is already a $drc_results sub-list (user passed that level directly)
    if (!plate %in% names(results)) {
      stop("Plate '", plate, "' not found. ",
           "Available plates: ", paste(names(results), collapse = ", "))
    }
    results <- results[[plate]]$drc_result
    if (is.null(results)) {
      stop("No drc_result found for plate '", plate, "'.")
    }
  }

  # ============================================================================
  # 1. DEPENDENCY CHECK AND INITIAL VALIDATION
  # ============================================================================

  # Check required packages
  required_packages <- c("ggplot2", "scales", "grDevices")
  missing_packages <- sapply(required_packages, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  })

  if (any(missing_packages)) {
    stop("The following packages are required: ",
         paste(required_packages[missing_packages], collapse = ", "))
  }

  # Load packages
  library(ggplot2)
  library(scales)

  # Validate input
  if (is.null(results$detailed_results) || length(results$detailed_results) == 0) {
    stop("Results object is empty or invalid")
  }

  if (verbose) {
    message("Processing multiple compounds for plotting...")
    message("Total compounds available: ", length(results$detailed_results))
  }

  # ============================================================================
  # 2. HELPER FUNCTIONS
  # ============================================================================

  detect_input_type <- function(input) {
    if (is.null(input)) return("none")
    if (grepl(":", input)) {
      parts <- strsplit(input, ":")[[1]]
      if (length(parts) == 2 && nchar(parts[1]) > 0 && nchar(parts[2]) > 0) {
        return("both")
      }
    }
    return("unknown")
  }

  find_matches <- function(pattern, compound_list, exact_first = TRUE) {
    all_names <- sapply(compound_list, function(x) {
      name <- strsplit(x$compound, " \\| ")[[1]][1]
      return(name)
    })

    searchable <- list()
    for (i in seq_along(all_names)) {
      name <- all_names[i]

      if (grepl(":", name)) {
        parts <- strsplit(name, ":")[[1]]
        target <- parts[1]
        compound <- ifelse(length(parts) >= 2, parts[2], name)

        compound_clean <- trimws(gsub("\\.2$", "", compound))
        target_clean <- trimws(gsub("\\.2$", "", target))

        searchable[[i]] <- list(
          full = name,
          target = target_clean,
          compound = compound_clean,
          target_compound = paste(target_clean, compound_clean, sep = ":")
        )
      } else {
        searchable[[i]] <- list(
          full = name,
          target = name,
          compound = name,
          target_compound = name
        )
      }
    }

    if (exact_first) {
      for (i in seq_along(searchable)) {
        if (searchable[[i]]$full == pattern) {
          return(list(matches = list(compound_list[[i]]),
                      indices = i,
                      type = "exact_full"))
        }
      }
    }

    if (grepl(":", pattern)) {
      pattern_parts <- strsplit(pattern, ":")[[1]]
      pattern_target <- pattern_parts[1]
      pattern_compound <- ifelse(length(pattern_parts) >= 2, pattern_parts[2], "")

      for (i in seq_along(searchable)) {
        if (searchable[[i]]$target == pattern_target &&
            searchable[[i]]$compound == pattern_compound) {
          return(list(matches = list(compound_list[[i]]),
                      indices = i,
                      type = "exact_target_compound"))
        }
      }
    }

    matches <- list()
    match_indices <- c()

    for (i in seq_along(searchable)) {
      if (grepl(pattern, searchable[[i]]$target, ignore.case = TRUE)) {
        matches[[length(matches) + 1]] <- compound_list[[i]]
        match_indices <- c(match_indices, i)
      }
      else if (grepl(pattern, searchable[[i]]$compound, ignore.case = TRUE)) {
        matches[[length(matches) + 1]] <- compound_list[[i]]
        match_indices <- c(match_indices, i)
      }
      else if (grepl(pattern, searchable[[i]]$full, ignore.case = TRUE)) {
        matches[[length(matches) + 1]] <- compound_list[[i]]
        match_indices <- c(match_indices, i)
      }
    }

    if (length(matches) > 0) {
      return(list(matches = matches, indices = match_indices, type = "partial"))
    }

    return(NULL)
  }

  smart_label_wrap <- function(labels, width = legend_label_wrap) {
    sapply(labels, function(label) {
      if (is.na(label) || nchar(label) <= width) return(label)

      if (grepl("[-_]", label)) {
        parts <- strsplit(label, "[-_]")[[1]]
        if (length(parts) > 1) {
          current_length <- 0
          break_points <- numeric(0)
          for (j in seq_along(parts)) {
            current_length <- current_length + nchar(parts[j]) + 1
            if (current_length <= width + 1) {
              break_points <- c(break_points, j)
            }
          }

          if (length(break_points) > 0) {
            break_point <- max(break_points)
            if (break_point < length(parts)) {
              line1 <- paste(parts[1:break_point], collapse = "-")
              line2 <- paste(parts[(break_point + 1):length(parts)], collapse = "-")
              return(paste(line1, line2, sep = "\n"))
            }
          }
        }
      }

      words <- strsplit(label, " ")[[1]]
      if (length(words) > 1) {
        lines <- character(0)
        current_line <- words[1]

        for (j in 2:length(words)) {
          test_line <- paste(current_line, words[j])
          if (nchar(test_line) <= width) {
            current_line <- test_line
          } else {
            lines <- c(lines, current_line)
            current_line <- words[j]
          }
        }
        lines <- c(lines, current_line)

        if (length(lines) <= 3) {
          return(paste(lines, collapse = "\n"))
        } else {
          mid <- ceiling(length(words) / 2)
          line1 <- paste(words[1:mid], collapse = " ")
          line2 <- paste(words[(mid + 1):length(words)], collapse = " ")
          return(paste(line1, line2, sep = "\n"))
        }
      }

      if (nchar(label) > width) {
        break_point <- ceiling(nchar(label) / 2)
        part1 <- substr(label, 1, break_point)
        part2 <- substr(label, break_point + 1, nchar(label))
        return(paste(part1, part2, sep = "\n"))
      }

      return(label)
    }, USE.NAMES = FALSE)
  }

  # ============================================================================
  # 3. COLOR PALETTE FUNCTIONS
  # ============================================================================

  generate_colors_from_palette <- function(n, palette = "hue") {

    # Scientific journal palettes
    nature_palettes <- list(
      nature = function(n) {
        # Classic Nature colors
        base_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      science = function(n) {
        # Vibrant Science magazine colors
        base_colors <- c("#3070B0", "#B03070", "#30B070", "#B07030", "#7030B0", "#70B030")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      cell = function(n) {
        # Cell journal colors
        base_colors <- c("#DC143C", "#4682B4", "#2E8B57", "#FF8C00", "#9370DB", "#20B2AA")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      plos = function(n) {
        # PLOS ONE colors
        base_colors <- c("#3498DB", "#E74C3C", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      elife = function(n) {
        # eLife journal colors
        base_colors <- c("#F04E4E", "#4EA5F0", "#4EF0A5", "#F0A54E", "#A54EF0", "#F04EA5")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      }
    )

    # Corporate palettes
    corporate_palettes <- list(
      ibm = function(n) {
        # IBM Design Library colors
        base_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000", "#000000")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      google = function(n) {
        # Google Colors
        base_colors <- c("#4285F4", "#EA4335", "#FBBC05", "#34A853", "#FF6D00", "#46BDC6")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      microsoft = function(n) {
        # Microsoft colors
        base_colors <- c("#F65314", "#7CBB00", "#00A1F1", "#FFBB00", "#A0A0A0", "#505050")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      twitter = function(n) {
        # Twitter/X colors
        base_colors <- c("#1DA1F2", "#14171A", "#657786", "#AAB8C2", "#E1E8ED", "#F5F8FA")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      }
    )

    # Colorblind-friendly palettes
    colorblind_palettes <- list(
      okabe_ito = function(n) {
        # Okabe & Ito palette (recommended for colorblindness)
        base_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      colorblind = function(n) {
        # Generic colorblind-friendly palette
        base_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      cud = function(n) {
        # Color Universal Design
        base_colors <- c("#0072B2", "#E69F00", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      tol = function(n) {
        # Paul Tol's colorblind-friendly palette
        if (n <= 1) return(c("#4477AA"))
        if (n == 2) return(c("#4477AA", "#EE6677"))
        if (n == 3) return(c("#4477AA", "#EE6677", "#228833"))
        if (n == 4) return(c("#4477AA", "#EE6677", "#228833", "#CCBB44"))
        if (n == 5) return(c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"))
        if (n == 6) return(c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"))
        # For n > 6, interpolate
        base <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377")
        return(colorRampPalette(base)(n))
      }
    )

    # Publisher-specific palettes
    journal_palettes <- list(
      bmc = function(n) {
        # BioMed Central
        base_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      frontiers = function(n) {
        # Frontiers journals
        base_colors <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      wiley = function(n) {
        # Wiley journals
        base_colors <- c("#5A9BD5", "#ED7D31", "#A5A5A5", "#FFC000", "#4472C4", "#70AD47")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      elsevier = function(n) {
        # Elsevier journals
        base_colors <- c("#F39800", "#DC143C", "#004080", "#009944", "#8B4513", "#4B0082")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      oxford = function(n) {
        # Oxford University Press
        base_colors <- c("#002147", "#8B0000", "#006A4E", "#FF6B35", "#5E2A84", "#008080")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      springer = function(n) {
        # Springer Nature
        base_colors <- c("#B22222", "#0066CC", "#228B22", "#FF8C00", "#9400D3", "#20B2AA")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      acs = function(n) {
        # American Chemical Society
        base_colors <- c("#0066CC", "#CC0000", "#009966", "#FF9900", "#660099", "#FF6600")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      },
      rsc = function(n) {
        # Royal Society of Chemistry
        base_colors <- c("#B31B1B", "#005F8C", "#2E8B57", "#FF7F0E", "#9467BD", "#17BECF")
        if (n <= length(base_colors)) {
          return(base_colors[1:n])
        } else {
          return(colorRampPalette(base_colors)(n))
        }
      }
    )

    # Gradient palettes
    gradient_palettes <- list(
      blue_red = function(n) colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
                                                "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(n),
      green_red = function(n) colorRampPalette(c("#1A9850", "#91CF60", "#D9EF8B", "#FEE08B", "#FC8D59", "#D73027"))(n),
      purple_orange = function(n) colorRampPalette(c("#542788", "#998EC3", "#D8DAEB", "#FEE0B6", "#F1A340", "#B35806"))(n),
      cool_warm = function(n) colorRampPalette(c("#3A5F8F", "#7BA0C0", "#C7D8E8", "#F1E3CF", "#EAAA7D", "#CC673B"))(n),
      blue_yellow = function(n) colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                                                   "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(n)
    )

    palettes <- c(
      # Base palettes
      hue = function(n) scales::hue_pal()(n),
      ggplot2 = function(n) scales::hue_pal()(n),
      default = function(n) scales::hue_pal()(n),

      # ColorBrewer qualitative palettes
      set1 = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 9)
          colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      set2 = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 8)
          colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      set3 = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 12)
          colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      dark2 = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 8)
          colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      paired = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 12)
          colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      accent = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 8)
          colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      pastel1 = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 9)
          colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },
      pastel2 = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          max_n <- min(n, 8)
          colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2"))(n)
        } else {
          scales::hue_pal()(n)
        }
      },

      # ColorBrewer sequential palettes
      blues = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(n)
        } else {
          colorRampPalette(c("#F7FBFF", "#08306B"))(n)
        }
      },
      reds = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(n)
        } else {
          colorRampPalette(c("#FFF5F0", "#67000D"))(n)
        }
      },
      greens = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(n)
        } else {
          colorRampPalette(c("#F7FCF5", "#00441B"))(n)
        }
      },
      purples = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))(n)
        } else {
          colorRampPalette(c("#FCFBFD", "#3F007D"))(n)
        }
      },
      oranges = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))(n)
        } else {
          colorRampPalette(c("#FFF5EB", "#7F2704"))(n)
        }
      },
      greys = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(n)
        } else {
          colorRampPalette(c("#FFFFFF", "#000000"))(n)
        }
      },

      # ColorBrewer diverging palettes
      spectral = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n)
        } else {
          colorRampPalette(c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B",
                             "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD"))(n)
        }
      },
      rdylbu = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n)
        } else {
          colorRampPalette(c("#D73027", "#FC8D59", "#FEE090",
                             "#E0F3F8", "#91BFDB", "#4575B4"))(n)
        }
      },
      rdylgn = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(n)
        } else {
          colorRampPalette(c("#D73027", "#FC8D59", "#FEE08B",
                             "#D9EF8B", "#91CF60", "#1A9850"))(n)
        }
      },
      piyg = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG"))(n)
        } else {
          colorRampPalette(c("#C51B7D", "#E9A3C9", "#FDE0EF",
                             "#E6F5D0", "#A6DBA0", "#008837"))(n)
        }
      },
      prgn = function(n) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          colorRampPalette(RColorBrewer::brewer.pal(11, "PRGn"))(n)
        } else {
          colorRampPalette(c("#762A83", "#9970AB", "#C2A5CF",
                             "#E7D4E8", "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837"))(n)
        }
      },

      # Viridis palettes (perceptually uniform)
      viridis = function(n) {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::viridis(n)
        } else {
          colorRampPalette(c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"))(n)
        }
      },
      magma = function(n) {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::magma(n)
        } else {
          colorRampPalette(c("#000004", "#2D1263", "#721F81", "#B63679", "#F8765C", "#FCFDBF"))(n)
        }
      },
      inferno = function(n) {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::inferno(n)
        } else {
          colorRampPalette(c("#000004", "#1F0C48", "#550F6D", "#A52C60", "#E7683A", "#FCFDBF"))(n)
        }
      },
      plasma = function(n) {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::plasma(n)
        } else {
          colorRampPalette(c("#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786",
                             "#D8576B", "#ED7953", "#FA9E3B", "#FDC926", "#F0F921"))(n)
        }
      },

      # Classic R palettes
      rainbow = function(n) rainbow(n),
      heat = function(n) heat.colors(n),
      terrain = function(n) terrain.colors(n),
      topo = function(n) topo.colors(n),
      cm = function(n) cm.colors(n)
    )

    # Combine all palette groups
    palettes <- c(palettes, nature_palettes, corporate_palettes, colorblind_palettes, journal_palettes, gradient_palettes)

    # Use the requested palette if it exists
    if (!is.null(palette) && palette %in% names(palettes)) {
      return(palettes[[palette]](n))
    }

    # Fallback: try RColorBrewer if palette name is not in the list above
    if (!is.null(palette) && requireNamespace("RColorBrewer", quietly = TRUE)) {
      if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
        max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
        return(colorRampPalette(RColorBrewer::brewer.pal(max_colors, palette))(n))
      }
    }

    # Final fallback to default hue palette
    warning("Palette '", palette, "' not recognized. Using default hue palette.")
    return(scales::hue_pal()(n))
  }

  list_available_palettes <- function() {
    palettes <- c(
      # Base palettes
      "hue", "ggplot2",

      # ColorBrewer qualitative
      "set1", "set2", "set3", "dark2", "paired", "accent", "pastel1", "pastel2",

      # ColorBrewer sequential
      "blues", "reds", "greens", "purples", "oranges", "greys",

      # ColorBrewer diverging
      "spectral", "rdylbu", "rdylgn", "piyg", "prgn",

      # Viridis
      "viridis", "magma", "inferno", "plasma",

      # Scientific journals
      "nature", "science", "cell", "plos", "elife",
      "bmc", "frontiers", "wiley", "elsevier", "oxford", "springer", "acs", "rsc",

      # Corporate
      "ibm", "google", "microsoft", "twitter",

      # Colorblind-friendly
      "okabe_ito", "colorblind", "cud", "tol",

      # Gradient
      "blue_red", "green_red", "purple_orange", "cool_warm", "blue_yellow",

      # Classic R
      "rainbow", "heat", "terrain", "topo", "cm"
    )
    return(palettes)
  }

  # ============================================================================
  # 4. BUILD COMPOUND LIST
  # ============================================================================

  build_compound_list <- function(results) {
    compound_list <- list()

    for (i in seq_along(results$detailed_results)) {
      result <- results$detailed_results[[i]]

      if (is.null(result$model) || !isTRUE(result$success)) next

      compound_name <- strsplit(result$compound, " \\| ")[[1]][1]

      if (grepl(":", compound_name)) {
        name_parts <- strsplit(compound_name, ":")[[1]]
        target <- name_parts[1]
        compound <- ifelse(length(name_parts) >= 2, name_parts[2], target)
      } else {
        target <- compound_name
        compound <- compound_name
      }

      target_clean <- trimws(gsub("\\.2$", "", target))
      compound_clean <- trimws(gsub("\\.2$", "", compound))

      compound_list[[i]] <- list(
        index = i,
        full_name = compound_name,
        target = target_clean,
        compound = compound_clean,
        target_compound = paste(target_clean, compound_clean, sep = ":"),
        result = result,
        success = TRUE
      )
    }

    return(compound_list)
  }

  compound_list <- build_compound_list(results)

  if (length(compound_list) == 0) {
    stop("No compounds with successful model fits found")
  }

  if (verbose) {
    message("Valid compounds with successful fits: ", length(compound_list))
  }

  # ============================================================================
  # 5. SELECT COMPOUNDS
  # ============================================================================

  selected_compounds <- list()
  match_type <- "unknown"
  input_used <- NULL
  selected_indices <- c()

  if (!is.null(position)) {
    if (position > length(compound_list)) {
      stop("Position ", position, " is out of range. Only ", length(compound_list), " compounds available.")
    }
    selected_compounds[[1]] <- compound_list[[position]]
    selected_indices <- position
    match_type <- "position"
    input_used <- paste("Position", position)

    if (verbose) {
      message("Selected by position ", position, ": ", compound_list[[position]]$target_compound)
    }
  }
  else if (!is.null(target_compound)) {
    input_type <- detect_input_type(target_compound)
    input_used <- target_compound

    if (verbose) message("Input '", target_compound, "' detected as type: ", input_type)

    result <- find_matches(target_compound, compound_list)

    if (!is.null(result)) {
      selected_compounds <- result$matches
      selected_indices <- result$indices
      match_type <- result$type

      if (verbose) {
        message("Found ", length(selected_compounds), " matches (match type: ", match_type, ")")
        for (i in seq_along(selected_compounds)) {
          message("  Match ", i, ": ", selected_compounds[[i]]$target_compound)
        }
      }
    } else {
      stop("No match found for '", target_compound, "'.")
    }
  }
  else if (!is.null(compound_indices)) {
    if (max(compound_indices) > length(compound_list)) {
      stop("Compound indices exceed available compounds")
    }

    for (idx in compound_indices) {
      if (idx <= length(compound_list) && !is.null(compound_list[[idx]])) {
        selected_compounds[[length(selected_compounds) + 1]] <- compound_list[[idx]]
        selected_indices <- c(selected_indices, idx)
      }
    }
    match_type <- "indices"
    input_used <- paste(compound_indices, collapse = ", ")

    if (verbose) {
      message("Selected by indices: ", length(selected_compounds), " compounds")
    }
  }
  else {
    selected_compounds <- compound_list
    selected_indices <- seq_along(compound_list)
    match_type <- "all"

    if (verbose) message("No selection criteria provided. Using all valid compounds.")
  }

  if (length(selected_compounds) == 0) {
    stop("No compounds selected for plotting.")
  }

  n_compounds <- length(selected_compounds)

  if (verbose) {
    message("Selected ", n_compounds, " compounds for plotting (match type: ", match_type, ")")
  }

  # ============================================================================
  # 6. EXTRACT DATA
  # ============================================================================

  extract_compound_data <- function(selected_compounds) {
    curve_data_list <- list()
    point_data_list <- list()

    for (i in seq_along(selected_compounds)) {
      comp_info <- selected_compounds[[i]]
      result <- comp_info$result

      data <- result$data
      valid_data <- data[is.finite(data$log_inhibitor) & is.finite(data$response), ]

      if (nrow(valid_data) < 2) {
        warning("Compound ", comp_info$target_compound, ": insufficient data - skipping")
        next
      }

      x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)
      x_seq <- seq(x_range[1], x_range[2], length.out = 100)

      curve_df <- data.frame(
        log_inhibitor = x_seq,
        response = predict(result$model, newdata = data.frame(log_inhibitor = x_seq)),
        compound = comp_info$target_compound,
        compound_index = i,
        target = comp_info$target,
        compound_name = comp_info$compound,
        full_name = comp_info$full_name
      )

      curve_data_list[[i]] <- curve_df

      conc_levels <- unique(valid_data$log_inhibitor)

      point_stats <- do.call(rbind, lapply(conc_levels, function(conc) {
        conc_data <- valid_data[valid_data$log_inhibitor == conc, ]
        if (nrow(conc_data) > 0) {
          data.frame(
            log_inhibitor = conc,
            mean_response = mean(conc_data$response, na.rm = TRUE),
            sd_response = sd(conc_data$response, na.rm = TRUE),
            n_points = nrow(conc_data),
            compound = comp_info$target_compound,
            compound_index = i,
            target = comp_info$target,
            compound_name = comp_info$compound,
            full_name = comp_info$full_name
          )
        }
      }))

      point_data_list[[i]] <- point_stats
    }

    curve_data_list <- Filter(Negate(is.null), curve_data_list)
    point_data_list <- Filter(Negate(is.null), point_data_list)

    if (length(curve_data_list) == 0) {
      stop("No valid data to plot")
    }

    list(
      curves = do.call(rbind, curve_data_list),
      points = if (length(point_data_list) > 0) do.call(rbind, point_data_list) else data.frame()
    )
  }

  plot_data <- extract_compound_data(selected_compounds)
  n_valid_compounds <- length(unique(plot_data$curves$compound))

  if (verbose) {
    message("Plot data prepared: ", nrow(plot_data$curves), " curve points")
    if (nrow(plot_data$points) > 0) {
      message("  and ", nrow(plot_data$points), " data points for ", n_valid_compounds, " compounds")
    }
  }

  # ============================================================================
  # 7. ANALYZE DATA
  # ============================================================================

  unique_targets <- unique(plot_data$curves$target)
  unique_compounds <- unique(plot_data$curves$compound_name)

  if (verbose) {
    message("Analysis for title generation:")
    message("  Unique targets: ", paste(unique_targets, collapse = ", "))
    message("  Unique compounds: ", paste(unique_compounds, collapse = ", "))
    message("  Match type: ", match_type)
    message("  Input used: ", ifelse(is.null(input_used), "None", input_used))
  }

  is_normalized <- if (!is.null(results$normalized)) results$normalized else FALSE

  y_axis_title_final <- if (is.null(y_axis_title)) {
    ifelse(is_normalized, "Normalized BRET ratio [%]", "BRET ratio")
  } else {
    y_axis_title
  }

  # ============================================================================
  # 8. GENERATE TITLE
  # ============================================================================

  generate_intelligent_title <- function(match_type, input_used, unique_targets,
                                         unique_compounds, n_compounds) {

    if (!is.null(plot_title)) {
      if (verbose) message("Using custom plot title: ", plot_title)
      return(plot_title)
    }

    if (n_compounds == 1) {
      if (length(unique_compounds) == 1 && unique_compounds[1] != "") {
        title <- unique_compounds[1]
        if (verbose) message("Title strategy: Single compound -> ", title)
        return(title)
      }
    }

    if (length(unique_compounds) == 1 && unique_compounds[1] != "") {
      title <- unique_compounds[1]
      if (verbose) message("Title strategy: Single compound across targets -> ", title)
      return(title)
    }

    if (length(unique_targets) == 1 && unique_targets[1] != "") {
      title <- unique_targets[1]
      if (verbose) message("Title strategy: Single target -> ", title)
      return(title)
    }

    if (match_type %in% c("exact_full", "exact_target_compound", "partial")) {
      if (!is.null(input_used)) {
        if (grepl(":", input_used)) {
          parts <- strsplit(input_used, ":")[[1]]
          if (length(parts) >= 2) {
            title <- parts[2]
          } else {
            title <- input_used
          }
        } else {
          title <- input_used
        }
        if (verbose) message("Title strategy: Matched input -> ", title)
        return(title)
      }
    }

    if (n_compounds > 1) {
      title <- paste("Multiple Compounds (", n_compounds, ")", sep = "")
    } else {
      title <- "Dose-Response Curve"
    }

    if (verbose) message("Title strategy: Default -> ", title)
    return(title)
  }

  plot_title_final <- generate_intelligent_title(match_type, input_used, unique_targets,
                                                 unique_compounds, n_valid_compounds)

  # ============================================================================
  # 9. SET LEGEND TITLE AND LABELS
  # ============================================================================

  # Generate smart legend labels
  generate_smart_legend_labels <- function(compound_labels, unique_targets, unique_compounds) {

    # Extract target and compound from each label (format "target:compound")
    extract_parts <- function(label) {
      if (grepl(":", label)) {
        parts <- strsplit(label, ":")[[1]]
        return(list(target = parts[1], compound = parts[2]))
      } else {
        return(list(target = label, compound = label))
      }
    }

    parts_list <- lapply(compound_labels, extract_parts)
    targets <- sapply(parts_list, function(x) x$target)
    compounds <- sapply(parts_list, function(x) x$compound)

    all_targets_same  <- length(unique(targets))   == 1
    all_compounds_same <- length(unique(compounds)) == 1

    smart_labels <- character(length(compound_labels))

    if (all_targets_same && !all_compounds_same) {
      # Same target, different compounds -> show only compound names
      if (verbose) message("Legend strategy: Same target, different compounds -> showing only compound names")
      smart_labels <- compounds
    } else if (!all_targets_same && all_compounds_same) {
      # Different targets, same compound -> show only target names
      if (verbose) message("Legend strategy: Different targets, same compound -> showing only target names")
      smart_labels <- targets
    } else {
      # Both differ (or both equal) -> show full target:compound
      if (verbose) message("Legend strategy: Showing full target:compound format")
      smart_labels <- compound_labels
    }

    return(smart_labels)
  }

  compound_labels <- unique(plot_data$curves$compound)

  smart_legend_names <- generate_smart_legend_labels(compound_labels, unique_targets, unique_compounds)

  wrapped_labels <- smart_label_wrap(smart_legend_names, legend_label_wrap)

  legend_title_final <- legend_title

  # Factor compounds to preserve order
  plot_data$curves$compound <- factor(plot_data$curves$compound,
                                      levels = compound_labels)
  if (nrow(plot_data$points) > 0) {
    plot_data$points$compound <- factor(plot_data$points$compound,
                                        levels = compound_labels)
  }

  # ============================================================================
  # 10. CONFIGURE SHAPES AND COLORS
  # ============================================================================

  # Point shape selection
  if (is.null(point_shapes)) {
    optimal_shapes <- c(16, 17, 15, 18, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14)

    if (n_valid_compounds <= 6) {
      point_shapes <- optimal_shapes[1:n_valid_compounds]
    } else {
      point_shapes <- rep(optimal_shapes, length.out = n_valid_compounds)
    }
  } else {
    point_shapes <- rep(point_shapes, length.out = n_valid_compounds)
  }

  determine_colors <- function(colors_param, color_palette_param, n_colors) {

    # colors = TRUE: use automatic hue palette
    if (is.logical(colors_param) && isTRUE(colors_param)) {
      if (verbose) message("Using automatic colors (colors = TRUE)")
      return(generate_colors_from_palette(n_colors, "hue"))
    }

    # Custom color vector supplied via 'colors'
    if (!is.null(colors_param) && is.character(colors_param)) {
      if (length(colors_param) < n_colors) {
        warning("Number of colors provided (", length(colors_param),
                ") is less than number of compounds (", n_colors,
                "). Recycling colors.")
        return(rep(colors_param, length.out = n_colors))
      } else {
        return(colors_param[1:n_colors])
      }
    }

    # Named palette specified via 'color_palette'
    if (!is.null(color_palette_param)) {
      if (verbose) {
        message("Using color palette: ", color_palette_param)

        # Informational messages about the palette type
        if (color_palette_param %in% c("set1", "set2", "set3", "dark2", "paired", "accent", "pastel1", "pastel2")) {
          message("  Note: ColorBrewer qualitative palette - good for distinct categories")
        } else if (color_palette_param %in% c("blues", "reds", "greens", "purples", "oranges", "greys")) {
          message("  Note: ColorBrewer sequential palette - good for ordered data")
        } else if (color_palette_param %in% c("spectral", "rdylbu", "rdylgn", "piyg", "prgn")) {
          message("  Note: ColorBrewer diverging palette - good for showing deviation from a midpoint")
        } else if (color_palette_param %in% c("viridis", "magma", "inferno", "plasma")) {
          message("  Note: Perceptually uniform palette - good for scientific visualization")
        } else if (color_palette_param %in% c("nature", "science", "cell", "plos", "elife")) {
          message("  Note: Journal-specific palette (", color_palette_param, ")")
        } else if (color_palette_param %in% c("bmc", "frontiers", "wiley", "elsevier", "oxford", "springer", "acs", "rsc")) {
          message("  Note: Publisher-specific palette (", color_palette_param, ")")
        } else if (color_palette_param %in% c("ibm", "google", "microsoft", "twitter")) {
          message("  Note: Corporate brand palette (", color_palette_param, ")")
        } else if (color_palette_param %in% c("okabe_ito", "colorblind", "cud", "tol")) {
          message("  Note: Colorblind-friendly palette (", color_palette_param, ")")
        } else if (color_palette_param %in% c("blue_red", "green_red", "purple_orange", "cool_warm", "blue_yellow")) {
          message("  Note: Gradient palette for heatmaps/continuous data")
        }
      }
      return(generate_colors_from_palette(n_colors, color_palette_param))
    }

    # Fallback to default hue palette
    if (verbose) message("Using default hue palette")
    return(generate_colors_from_palette(n_colors, "hue"))
  }

  colors_final <- determine_colors(colors, color_palette, n_valid_compounds)

  if (verbose && n_valid_compounds <= 10) {
    message("Colors assigned:")
    for (i in seq_along(compound_labels)) {
      message("  ", compound_labels[i], ": ", colors_final[i])
    }
  } else if (verbose && n_valid_compounds > 10) {
    message("Colors assigned to ", n_valid_compounds, " compounds (first 10 shown):")
    for (i in 1:min(10, n_valid_compounds)) {
      message("  ", compound_labels[i], ": ", colors_final[i])
    }
  }

  # ============================================================================
  # 11. CREATE THE PLOT
  # ============================================================================

  calculate_x_limits <- function(curve_data, point_data) {
    all_x <- c(curve_data$log_inhibitor)
    if (nrow(point_data) > 0) {
      all_x <- c(all_x, point_data$log_inhibitor)
    }
    valid_x <- all_x[is.finite(all_x)]

    if (length(valid_x) == 0) return(c(-10, -2))

    x_range <- range(valid_x)
    x_margin <- diff(x_range) * 0.05
    return(c(x_range[1] - x_margin, x_range[2] + x_margin))
  }

  # Resolve x_limits to log10 molar scale.
  # x_limits_scale controls the unit of the user-supplied values:
  #   "log10"  - already log10 molar (default, no conversion)
  #   "molar"  - raw molar (e.g. 1e-9); converted via log10()
  #   "uM"     - micromolar (e.g. 0.001 to 25); converted via log10(x * 1e-6)
  #   "nM"     - nanomolar; converted via log10(x * 1e-9)
  resolve_x_limits <- function(xl) {
    if (is.null(xl)) return(calculate_x_limits(plot_data$curves, plot_data$points))
    if (length(xl) != 2 || !all(is.finite(xl)))
      stop("`x_limits` must be a finite numeric vector of length 2.")
    scale <- match.arg(x_limits_scale, c("log10", "molar", "uM", "nM"))
    converted <- switch(scale,
                        log10 = xl,
                        molar = log10(xl),
                        uM    = log10(xl * 1e-6),
                        nM    = log10(xl * 1e-9)
    )
    if (scale != "log10" && verbose)
      message(sprintf("`x_limits` converted from %s to log10 molar: [%.2f, %.2f]",
                      scale, converted[1], converted[2]))
    converted
  }

  calculate_y_limits <- function(curve_data, point_data) {
    all_y <- c(curve_data$response)
    if (nrow(point_data) > 0) {
      # include the full error-bar extent when available
      if (!is.null(point_data$sd_response)) {
        all_y <- c(all_y,
                   point_data$mean_response + point_data$sd_response,
                   point_data$mean_response - point_data$sd_response)
      } else {
        all_y <- c(all_y, point_data$mean_response)
      }
    }
    valid_y <- all_y[is.finite(all_y)]
    if (length(valid_y) == 0) return(c(0, 150))
    y_range  <- range(valid_y)
    y_margin <- diff(y_range) * 0.08
    return(c(y_range[1] - y_margin, y_range[2] + y_margin))
  }

  x_limits_final <- resolve_x_limits(x_limits)
  y_limits_final <- if (!is.null(y_limits)) y_limits else
    calculate_y_limits(plot_data$curves, plot_data$points)

  x_axis_title_final <- if (!is.null(x_axis_title)) x_axis_title else
    expression(paste("Log"[10], " Concentration [M]"))

  point_size <- if (!is.null(point_size)) point_size else if (n_valid_compounds > 15) 3 else if (n_valid_compounds > 8) 3.5 else 4

  if (is.null(legend_ncol)) {
    legend_ncol <- if (n_valid_compounds > 10) 2 else 1
  }

  # Smart defaults for legend text size
  if (is.null(legend_text_size)) {
    legend_text_size <- if (n_valid_compounds > 15) 9 else if (n_valid_compounds > 8) 10 else 11
  }

  p <- ggplot2::ggplot() +
    ggplot2::labs(
      x = x_axis_title_final,
      y = y_axis_title_final,
      title = plot_title_final,
      color = legend_title_final
    ) +
    ggplot2::coord_cartesian(xlim = x_limits_final, ylim = y_limits_final)

  p <- p +
    ggplot2::geom_line(data = plot_data$curves,
                       ggplot2::aes(x = log_inhibitor, y = response,
                                    group = compound, color = compound),
                       linewidth = curve_linewidth, alpha = curve_alpha)

  if (nrow(plot_data$points) > 0) {
    if (shape_by_compound) {
      p <- p +
        ggplot2::geom_point(data = plot_data$points,
                            ggplot2::aes(x = log_inhibitor, y = mean_response,
                                         shape = compound, color = compound),
                            size = point_size)

      p <- p + ggplot2::scale_shape_manual(
        values = setNames(point_shapes[1:n_valid_compounds], compound_labels),
        guide = "none"
      )
    } else {
      p <- p +
        ggplot2::geom_point(data = plot_data$points,
                            ggplot2::aes(x = log_inhibitor, y = mean_response,
                                         color = compound),
                            size = point_size, shape = 16)
    }

    if (show_error_bars && "sd_response" %in% colnames(plot_data$points)) {
      p <- p +
        ggplot2::geom_errorbar(data = plot_data$points,
                               ggplot2::aes(x = log_inhibitor,
                                            ymin = mean_response - sd_response,
                                            ymax = mean_response + sd_response,
                                            color = compound),
                               width = error_bar_width,
                               linewidth = 0.5, alpha = 0.6)
    }
  }

  # Optional IC50 vertical reference lines
  if (show_ic50_lines) {
    ic50_rows <- lapply(seq_along(selected_compounds), function(i) {
      res       <- selected_compounds[[i]]$result
      log_ic50  <- if (!is.null(res$parameters)) res$parameters$Value[3] else NA_real_
      if (is.na(log_ic50) || !is.finite(log_ic50)) return(NULL)
      data.frame(
        log_ic50 = log_ic50,
        compound = selected_compounds[[i]]$target_compound,
        stringsAsFactors = FALSE
      )
    })
    ic50_data <- do.call(rbind, Filter(Negate(is.null), ic50_rows))
    if (!is.null(ic50_data) && nrow(ic50_data) > 0) {
      p <- p + ggplot2::geom_vline(
        data = ic50_data,
        ggplot2::aes(xintercept = log_ic50, color = compound),
        linetype = "dashed", linewidth = 0.5, alpha = 0.8
      )
    }
  }

  # Apply colors and smart legend labels
  p <- p + ggplot2::scale_color_manual(
    values = setNames(colors_final[1:n_valid_compounds], compound_labels),
    labels = setNames(wrapped_labels, compound_labels)
  )

  # Configure legend guide
  guide_args <- list(
    ncol = legend_ncol,
    override.aes = list(
      size = point_size,
      linetype = 1,
      linewidth = 0.8,
      fill = NA
    ),
    title = legend_title_final
  )

  # Add shapes to override when shape_by_compound is active
  if (nrow(plot_data$points) > 0 && shape_by_compound) {
    guide_args$override.aes$shape <- point_shapes[1:n_valid_compounds]
  }

  p <- p + ggplot2::guides(color = do.call(ggplot2::guide_legend, guide_args))

  # ============================================================================
  # 12. CONFIGURE PLOT THEME
  # ============================================================================

  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = ifelse(show_legend, legend_position, "none"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                         size = if (!is.null(plot_title_size)) plot_title_size else axis_title_size + 2),
      axis.text = ggplot2::element_text(color = axis_text_color, size = axis_text_size),
      axis.title = ggplot2::element_text(color = axis_title_color, size = axis_title_size, face = "bold"),
      axis.line.x.bottom = ggplot2::element_line(color = axis_line_color, linewidth = 0.5),
      axis.line.y.left   = ggplot2::element_line(color = axis_line_color, linewidth = 0.5),
      axis.line.x.top    = ggplot2::element_blank(),
      axis.line.y.right  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = axis_line_color),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      legend.text = ggplot2::element_text(size = legend_text_size, lineheight = 0.8),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_blank()
    )

  if (transparent_background) {
    base_theme <- base_theme + ggplot2::theme(
      legend.key        = ggplot2::element_rect(fill = NA, color = NA),
      legend.background = ggplot2::element_rect(fill = NA, color = NA),
      panel.background  = ggplot2::element_rect(fill = NA, color = NA),
      plot.background   = ggplot2::element_rect(fill = NA, color = NA)
    )
  }

  if (show_border) {
    base_theme <- base_theme + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = axis_line_color, fill = NA, linewidth = 0.5)
    )
  }

  if (!show_grid) {
    base_theme <- base_theme + ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  } else {
    base_theme <- base_theme + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_line(color = "grey95", linewidth = 0.25)
    )
  }

  p <- p + base_theme

  # ============================================================================
  # 13. DISPLAY AND SAVE
  # ============================================================================

  print(p)

  if (!is.null(save_plot)) {
    if (is.character(save_plot)) {
      filename <- save_plot
    } else if (is.logical(save_plot) && save_plot) {
      indices_str <- paste(selected_indices, collapse = "_")
      filename <- paste0("multiple_compounds_", indices_str, ".png")
    } else {
      stop("save_plot must be either a file path or TRUE for auto-naming")
    }

    plot_dir <- dirname(filename)
    if (plot_dir != "." && !dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (transparent_background && grepl("\\.jpe?g$", filename, ignore.case = TRUE))
      warning("transparent_background = TRUE has no effect when saving as JPEG: ",
              "JPEG does not support transparency. Use PNG, PDF, or SVG instead.")

    ggplot2::ggsave(filename, plot = p, width = plot_width,
                    height = plot_height, dpi = plot_dpi,
                    bg = if (transparent_background) "transparent" else "white")

    if (verbose) message("Plot saved as: ", normalizePath(filename))
  }

  # ============================================================================
  # 14. RETURN METADATA
  # ============================================================================

  metadata <- list(
    selected_compounds = compound_labels,
    smart_legend_names = smart_legend_names,
    n_compounds = n_valid_compounds,
    selected_indices = selected_indices,
    unique_targets = unique_targets,
    unique_compounds = unique_compounds,
    match_type = match_type,
    input_used = input_used,
    is_normalized = is_normalized,
    plot_title = plot_title_final,
    y_axis_title = y_axis_title_final,
    point_shapes = point_shapes[1:n_valid_compounds],
    colors = colors_final[1:n_valid_compounds],
    color_palette_used = if (!is.null(color_palette)) color_palette else "hue",
    point_size = point_size,
    legend_position = legend_position,
    legend_ncol = legend_ncol,
    legend_title = legend_title_final,
    legend_text_size = legend_text_size,
    legend_title_size = legend_title_size,
    wrapped_labels = wrapped_labels,
    x_limits = x_limits,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    file_saved = if (!is.null(save_plot)) filename else NULL,
    shape_by_compound = shape_by_compound,
    show_grid = show_grid,
    available_palettes = list_available_palettes()
  )

  attr(p, "metadata") <- metadata

  if (verbose && !is.null(color_palette)) {
    message("\nTip: To see all available palettes, run: names(attr(p, 'metadata')$available_palettes)")
  }

  return(p)
}
