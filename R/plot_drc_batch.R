#' Plot Dose-Response Curves (DRC) from Batch Analysis Results
#'
#' @description
#' This function generates publication-quality Dose-Response Curve (DRC) plots from a list of
#' batch processing results. It handles data extraction, statistical aggregation, and visualization
#' in a single step. It supports advanced features like scientific color palettes, faceting,
#' custom styling, and selective plotting of specific plates.
#'
#' @param batch_drc_results A named list containing the results of the DRC batch analysis.
#' Structure expected: \code{list(plate_name = list(drc_result = list(detailed_results = ...)))}.
#' Can also accept a wrapper object from \code{batch_drc_analysis()} which contains a \code{drc_results} element.
#' @param construct_compound Character string. The specific construct/compound to plot.
#' Supports two formats:
#' \itemize{
#' \item \strong{"Construct:Compound"}: Performs an exact match (e.g., "PAK3:IPA-3").
#' \item \strong{"SearchString"}: Performs a fuzzy search across constructs and compounds (e.g., "IPA-3").
#' }
#' Ignored if \code{position} is provided.
#' @param position Integer. Selects a compound by its numerical index in the list of unique results.
#' Useful for iterating through compounds blindly.
#' @param plates Character vector or numeric vector. Specifies which plates to include in the plot.
#' Options:
#' \itemize{
#' \item \code{NULL} (default): Include all plates
#' \item Character vector: Plate names to include (e.g., \code{c("plate_01", "plate_02")})
#' \item Numeric vector: Plate indices to include (e.g., \code{c(1, 3)} for first and third plates)
#' }
#' Useful for comparing specific replicates or excluding problematic plates.
#' @param y_limits Numeric vector of length 2 (e.g., \code{c(0, 150)}). Sets the Y-axis limits.
#' Applied using \code{coord_cartesian()} to maintain all data points (no filtering).
#' If \code{NULL} (default), limits are automatic.
#' \emph{Note:} Ignored if \code{facet_scales} contains "free".
#' @param colors Character vector or string. Controls the line/point colors. Accepts:
#' \itemize{
#' \item \strong{Manual colors}: A vector of color codes/names (e.g., \code{c("red", "blue")}).
#' \item \strong{Single color name}: Any valid R color name (e.g., "red", "blue", "darkgreen").
#' \item \strong{Viridis palettes}: "viridis", "magma", "plasma", "inferno", "cividis", "turbo", "rocket", "mako".
#' Requires \code{viridisLite} package.
#' \item \strong{RColorBrewer palettes}: Names like "Set1", "Set2", "Set3", "Pastel1", "Pastel2",
#' "Paired", "Dark2", "Accent". Requires \code{RColorBrewer} package.
#' \item \strong{ggsci palettes}: "jama", "nature" (npg), "lancet", "nejm", "aaas", "d3", "locuszoom",
#' "igv", "uchicago", "startrek", "tron", "futurama", "rickandmorty", "simpsons", "gsea".
#' Requires \code{ggsci} package.
#' }
#' Palettes are automatically interpolated if more colors are needed than the palette provides.
#' @param point_shapes Integer vector. Custom shapes for points (0-25). Used only if \code{shape_by_duplicate = TRUE}.
#' @param show_error_bars Logical. If \code{TRUE} (default), adds error bars (mean +/- SD) to the points.
#' @param legend_position Character. Position of the legend ("right", "bottom", "top", "none").
#' @param show_legend Logical. Whether to display the legend.
#' @param shape_by_duplicate Logical. If \code{TRUE} and plotting a single compound group,
#' different replicates (plates) will have different point shapes.
#' @param show_grid Logical. If \code{TRUE}, displays grid lines. If \code{FALSE} (default),
#' produces a clean "classic" scientific style (axis lines only).
#' @param font_family Character. Font family for text elements.
#' Common options: "sans" (Arial/Helvetica), "serif" (Times New Roman), "mono" (Courier).
#' For custom fonts, use the \code{extrafont} package.
#' @param facet_by Character. Variable to split the plot into panels.
#' Options: \code{"compound"}, \code{"construct"}, \code{"plate"}. Default is \code{NULL} (single plot).
#' @param facet_ncol Integer. Number of columns for the faceted grid.
#' @param facet_scales Character. Scaling of axes in facets.
#' Options: \code{"fixed"} (default), \code{"free"}, \code{"free_y"}, \code{"free_x"}.
#' @param save_plot Character string. File path to save the plot (e.g., "plot.png", "plot.pdf").
#' If \code{NULL}, the plot is not saved automatically.
#' @param plot_width Numeric. Width of the saved plot in inches (default: 12).
#' @param plot_height Numeric. Height of the saved plot in inches (default: 8).
#' @param plot_dpi Integer. Resolution of the saved plot (default: 600).
#' @param plot_title Character. Custom main title. If \code{NULL}, generates a smart title based on selection.
#' @param title_size Numeric; size of plot title.
#' @param legend_title Character. Custom legend title. If \code{NULL}, generates based on data.
#' @param legend_title_size Numeric; size of legend title.
#' @param legend_text_size Numeric; size of legend text.
#' @param y_axis_title Character. Custom Y-axis title.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages to the console,
#' including which plates are being filtered when using the \code{plates} parameter.
#' @param axis_text_color Character. Color of axis tick labels (default: "black").
#' @param axis_text_size Numeric. Size of axis tick labels (default: 12).
#' @param axis_title_color Character. Color of axis titles (default: "black").
#' @param axis_title_size Numeric. Size of axis titles (default: 14).
#'
#' @return A \code{ggplot2} object. This allows further modification using standard ggplot syntax
#' (e.g., \code{plot + theme_dark()}).
#'
#' @details
#' \strong{Data Extraction:} The function iterates through the nested list structure of \code{batch_drc_results},
#' extracting raw data points and model predictions for valid fits. Curve predictions are generated
#' only within the range of experimental data (no extrapolation).
#'
#' \strong{Curve Behavior:} Model curves are plotted only within the concentration range of the
#' experimental data for each curve, avoiding extrapolation beyond measured points.
#'
#' \strong{Error Bars:} When multiple replicates exist at the same concentration, error bars
#' represent mean ± standard deviation.
#'
#' \strong{Faceting vs. Limits:} If \code{facet_scales} is set to "free" or "free_y", the global \code{y_limits}
#' argument will be ignored to allow each panel to scale independently.
#'
#' \strong{Automatic Title Generation:} When \code{plot_title = NULL}, the function generates
#' appropriate titles:
#' - Single compound: Compound name
#' - Multiple compounds, single construct: Construct name
#' - Multiple compounds, single compound name: Compound name
#' - Otherwise: "Comparative DRC Analysis"
#' The title does not indicate which plates are plotted when using the \code{plates} parameter.
#'
#' \strong{Interoperability:} The function works seamlessly with results from \code{batch_drc_analysis()}
#' and can handle both normalized and raw BRET ratio data automatically.
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_pdb")
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
#' cc <- drc_res$drc_results$plate_01$drc_result$summary_table$Compound[1]
#' p  <- plot_drc_batch(drc_res, construct_compound = cc, verbose = FALSE)
#' inherits(p, "ggplot")
#' }
#' @seealso
#' \code{\link{batch_drc_analysis}} for generating batch DRC results
#' \code{\link{batch_save_all_drc_plots}} for automatically saving all DRC plots
#' \code{\link[ggplot2]{ggplot}} for underlying plotting functionality
#'
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridisLite viridis
#' @importFrom ggsci pal_jama pal_npg pal_lancet pal_nejm pal_aaas pal_d3
#' @importFrom ggsci pal_locuszoom pal_igv pal_uchicago pal_startrek pal_tron
#' @importFrom ggsci pal_futurama pal_rickandmorty pal_simpsons pal_gsea
#' @export


plot_drc_batch <- function(batch_drc_results,
                           construct_compound = NULL,
                           position = NULL,
                           plates = NULL,
                           y_limits = NULL,
                           colors = NULL,
                           point_shapes = NULL,
                           show_error_bars = TRUE,
                           legend_position = "right",
                           show_legend = TRUE,
                           shape_by_duplicate = FALSE,
                           show_grid = FALSE,
                           font_family = "sans",
                           facet_by = NULL,
                           facet_ncol = NULL,
                           facet_scales = "fixed",
                           save_plot = NULL,
                           plot_width = 12,
                           plot_height = 8,
                           plot_dpi = 600,
                           plot_title = NULL,
                           legend_title = NULL,
                           y_axis_title = NULL,
                           verbose = TRUE,
                           axis_text_color = "black",
                           axis_text_size = 12,
                           axis_title_color = "black",
                           axis_title_size = 14,
                           title_size = 16,
                           legend_title_size = 12,
                           legend_text_size = 10
) {

  # ============================================================================
  # 1. SETUP & NAMESPACE CHECKS
  # ============================================================================
  required_packages <- c("ggplot2", "scales", "dplyr")
  missing_packages <- sapply(required_packages, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  })

  if (any(missing_packages)) {
    stop("The following packages are required: ", paste(required_packages[missing_packages], collapse = ", "))
  }

  if (is.list(batch_drc_results) && "drc_results" %in% names(batch_drc_results)) {
    if (verbose) message("Detected analysis wrapper object. Extracting 'drc_results'...")
    batch_drc_results <- batch_drc_results$drc_results
  }

  if (length(batch_drc_results) == 0) stop("The DRC results list is empty.")
  if (verbose) message("Processing DRC batch results for plotting...")

  # ============================================================================
  # 2. HELPER FUNCTIONS
  # ============================================================================
  parse_compound_name <- function(name) {
    if (is.null(name)) return(list(construct = "Unknown", compound = "Unknown"))
    clean_name <- trimws(gsub("\\.\\d+$", "", name))

    if (grepl(" \\| ", clean_name)) {
      parts <- strsplit(clean_name, " \\| ")[[1]]
      clean_name <- parts[1]
    }

    if (grepl(":", clean_name)) {
      parts <- strsplit(clean_name, ":")[[1]]
      return(list(construct = trimws(parts[1]), compound = if(length(parts)>1) trimws(parts[2]) else trimws(parts[1])))
    } else {
      return(list(construct = clean_name, compound = clean_name))
    }
  }

  generate_palette <- function(user_input, n_needed) {
    if (is.null(user_input)) {
      return(scales::hue_pal()(n_needed))
    }

    if (length(user_input) > 1) {
      return(rep(user_input, length.out = n_needed))
    }

    pal_name <- tolower(user_input)

    if (pal_name %in% colors()) {
      return(rep(pal_name, n_needed))
    }

    viridis_palettes <- c("viridis", "magma", "plasma", "inferno", "cividis", "turbo", "rocket", "mako")
    if (pal_name %in% viridis_palettes) {
      if (requireNamespace("viridisLite", quietly = TRUE)) {
        return(viridisLite::viridis(n_needed, option = pal_name))
      } else {
        stop("Package 'viridisLite' is required for viridis palettes. Please install it: install.packages('viridisLite')")
      }
    }

    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      brewer_palettes <- tolower(rownames(RColorBrewer::brewer.pal.info))
      if (pal_name %in% brewer_palettes) {
        correct_name <- rownames(RColorBrewer::brewer.pal.info)[tolower(rownames(RColorBrewer::brewer.pal.info)) == pal_name]
        max_colors <- RColorBrewer::brewer.pal.info[correct_name, "maxcolors"]

        if (n_needed <= max_colors) {
          return(RColorBrewer::brewer.pal(n_needed, correct_name))
        } else {
          if (verbose) message("Note: Interpolating palette '", user_input,
                               "' (requested ", n_needed, " colors, palette has ", max_colors, ").")
          return(colorRampPalette(RColorBrewer::brewer.pal(max_colors, correct_name))(n_needed))
        }
      }
    } else if (pal_name %in% tolower(c("Set1", "Set2", "Set3", "Pastel1", "Pastel2", "Paired", "Dark2", "Accent"))) {
      stop("Package 'RColorBrewer' is required for RColorBrewer palettes. Please install it: install.packages('RColorBrewer')")
    }

    ggsci_mapping <- list(
      "jama" = list(pal_func = "pal_jama", default_n = 7),
      "nature" = list(pal_func = "pal_npg", default_n = 10),
      "lancet" = list(pal_func = "pal_lancet", default_n = 9),
      "nejm" = list(pal_func = "pal_nejm", default_n = 8),
      "aaas" = list(pal_func = "pal_aaas", default_n = 10),
      "d3" = list(pal_func = "pal_d3", default_n = 10),
      "locuszoom" = list(pal_func = "pal_locuszoom", default_n = 7),
      "igv" = list(pal_func = "pal_igv", default_n = 51),
      "uchicago" = list(pal_func = "pal_uchicago", default_n = 9),
      "startrek" = list(pal_func = "pal_startrek", default_n = 7),
      "tron" = list(pal_func = "pal_tron", default_n = 7),
      "futurama" = list(pal_func = "pal_futurama", default_n = 12),
      "rickandmorty" = list(pal_func = "pal_rickandmorty", default_n = 12),
      "simpsons" = list(pal_func = "pal_simpsons", default_n = 16),
      "gsea" = list(pal_func = "pal_gsea", default_n = 12)
    )

    if (pal_name %in% names(ggsci_mapping)) {
      if (!requireNamespace("ggsci", quietly = TRUE)) {
        stop("Package 'ggsci' is required for palette '", user_input,
             "'. Please install it: install.packages('ggsci')")
      }

      pal_info <- ggsci_mapping[[pal_name]]
      pal_func_name <- pal_info$pal_func
      default_n <- pal_info$default_n

      pal_func <- get(pal_func_name, asNamespace("ggsci"))

      if (n_needed > default_n) {
        if (verbose) message("Note: Palette '", user_input, "' has ", default_n,
                             " colors, but ", n_needed, " are needed. Using interpolated colors.")

        base_colors <- pal_func()(default_n)
        return(colorRampPalette(base_colors)(n_needed))
      } else {
        return(pal_func()(n_needed))
      }
    }

    stop("Palette '", user_input, "' not recognized. ",
         "Valid options are: \n",
         "  - Color names (e.g., 'red', 'blue')\n",
         "  - Viridis palettes: 'viridis', 'magma', 'plasma', 'inferno', 'cividis', 'turbo'\n",
         "  - RColorBrewer palettes (install.packages('RColorBrewer'))\n",
         "  - ggsci palettes (install.packages('ggsci')): 'jama', 'nature', 'lancet', 'nejm', 'aaas', etc.\n",
         "  - Custom vector of colors")
  }

  get_detailed_results <- function(plate_obj) {
    if (is.null(plate_obj$drc_result)) return(NULL)
    res <- plate_obj$drc_result$detailed_results
    if (is.null(res)) res <- plate_obj$drc_result$curve_results
    if (is.null(res)) res <- plate_obj$drc_result$fits
    if (is.list(res) && length(res) > 0) return(res)
    return(NULL)
  }

  # ============================================================================
  # 3. DATA EXTRACTION (COM FILTRO DE PLACAS)
  # ============================================================================
  raw_data_list <- list(); curve_data_list <- list(); counter <- 0

  normalized_info <- list()

  available_plates <- names(batch_drc_results)

  if (!is.null(plates)) {
    if (is.character(plates)) {
      selected_plates <- plates[plates %in% available_plates]
      if (length(selected_plates) == 0) {
        warning("None of the specified plates found. Available plates: ",
                paste(available_plates, collapse = ", "))
        selected_plates <- available_plates
      }
    }
    else if (is.numeric(plates)) {
      selected_plates <- available_plates[plates[plates <= length(available_plates)]]
      if (length(selected_plates) == 0) {
        warning("Invalid plate indices. Using all plates.")
        selected_plates <- available_plates
      }
    } else {
      warning("Invalid 'plates' argument. Using all plates.")
      selected_plates <- available_plates
    }

    if (verbose && length(selected_plates) < length(available_plates)) {
      message("Filtering plates: ", paste(selected_plates, collapse = ", "))
    }
  } else {
    selected_plates <- available_plates
  }

  for (plate_name in selected_plates) {
    plate_obj <- batch_drc_results[[plate_name]]
    plate_res_list <- get_detailed_results(plate_obj)
    if (is.null(plate_res_list)) next

    plate_normalized <- NULL
    if (!is.null(plate_obj$drc_result) && !is.null(plate_obj$drc_result$normalized)) {
      plate_normalized <- plate_obj$drc_result$normalized
    }

    for (i in seq_along(plate_res_list)) {
      res <- plate_res_list[[i]]
      has_success <- !is.null(res$success) && isTRUE(res$success)
      has_data <- !is.null(res$data) && nrow(res$data) > 0

      if (!has_success || !has_data) next

      comp_name <- res$compound
      if(is.null(comp_name)) comp_name <- paste0("Cmpd_", i)

      info <- parse_compound_name(comp_name)
      cc_key <- paste(info$construct, info$compound, sep = ":")
      unique_id <- paste(plate_name, i, sep = "_")
      counter <- counter + 1

      valid_data <- dplyr::filter(res$data, is.finite(log_inhibitor) & is.finite(response))
      valid_data <- dplyr::mutate(valid_data,
                                  plate = plate_name, construct = info$construct,
                                  compound = info$compound, construct_compound = cc_key,
                                  unique_id = unique_id)

      if (nrow(valid_data) >= 2) {
        raw_data_list[[counter]] <- valid_data

        normalized_info[[unique_id]] <- if (!is.null(plate_normalized)) {
          plate_normalized
        } else if (!is.null(res$normalized)) {
          res$normalized
        } else {
          FALSE
        }

        model_obj <- res$model
        if (!is.null(model_obj)) {
          try({
            x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)

            x_seq <- seq(x_range[1], x_range[2], length.out = 100)

            pred_y <- predict(model_obj, newdata = data.frame(log_inhibitor = x_seq))

            curve_data_list[[counter]] <- data.frame(
              log_inhibitor = x_seq, response = pred_y, plate = plate_name,
              construct = info$construct, compound = info$compound, construct_compound = cc_key, unique_id = unique_id
            )
          }, silent = TRUE)
        }
      }
    }
  }

  if (length(raw_data_list) == 0) stop("No valid data found.")

  df_raw_master <- dplyr::bind_rows(raw_data_list)
  df_curve_master <- dplyr::bind_rows(curve_data_list)

  # ============================================================================
  # 4. FILTERING (CONSTRUCT/COMPOUND E POSITION)
  # ============================================================================
  unique_combinations <- dplyr::select(df_raw_master, construct_compound, construct, compound, unique_id)
  unique_combinations <- dplyr::distinct(unique_combinations)

  selected_ids <- NULL; match_desc <- ""

  if (!is.null(position)) {
    unique_cc <- unique(unique_combinations$construct_compound)
    if (position > length(unique_cc)) stop("Invalid position index.")
    target_cc <- unique_cc[position]

    matches_pos <- dplyr::filter(unique_combinations, construct_compound == target_cc)
    selected_ids <- dplyr::pull(matches_pos, unique_id)
    match_desc <- paste("Pos:", position)

  } else if (!is.null(construct_compound)) {
    search_term <- construct_compound
    if (grepl(":", search_term)) {
      parts <- strsplit(search_term, ":")[[1]]
      construct_query <- trimws(parts[1]); compound_query <- if(length(parts) > 1) trimws(parts[2]) else ""

      matches <- dplyr::filter(unique_combinations, construct == construct_query, compound == compound_query)
      if (nrow(matches) == 0) {
        matches <- dplyr::filter(unique_combinations,
                                 tolower(construct) == tolower(construct_query),
                                 tolower(compound) == tolower(compound_query))
      }
    } else {
      matches <- dplyr::filter(unique_combinations,
                               grepl(search_term, compound, ignore.case=TRUE) |
                                 grepl(search_term, construct, ignore.case=TRUE))
    }
    if (nrow(matches) == 0) stop("No match found for: ", search_term)
    selected_ids <- matches$unique_id; match_desc <- search_term

  } else {
    first_cc <- unique_combinations$construct_compound[1]
    matches_first <- dplyr::filter(unique_combinations, construct_compound == first_cc)
    selected_ids <- dplyr::pull(matches_first, unique_id)
    match_desc <- "First compound"
  }

  plot_raw <- dplyr::filter(df_raw_master, unique_id %in% selected_ids)
  plot_curve <- dplyr::filter(df_curve_master, unique_id %in% selected_ids)

  selected_norm_flags <- sapply(selected_ids, function(id) {
    if (!is.null(normalized_info[[id]])) {
      normalized_info[[id]]
    } else {
      FALSE
    }
  }, USE.NAMES = FALSE)

  # ============================================================================
  # 5. DATA AGGREGATION & LEGEND LOGIC
  # ============================================================================
  plot_stats <- dplyr::group_by(plot_raw, plate, construct, compound, construct_compound, unique_id, log_inhibitor)
  plot_stats <- dplyr::summarise(plot_stats,
                                 mean_response = mean(response, na.rm = TRUE),
                                 sd_response = sd(response, na.rm = TRUE),
                                 .groups="drop")

  n_groups <- length(unique(plot_raw$construct_compound))
  n_distinct_constructs <- length(unique(plot_raw$construct))
  n_distinct_compounds <- length(unique(plot_raw$compound))

  plate_levels <- unique(plot_raw$plate)
  nums <- suppressWarnings(as.numeric(gsub("\\D", "", plate_levels)))
  if (!any(is.na(nums)) && length(nums) == length(plate_levels)) plate_levels <- plate_levels[order(nums)] else plate_levels <- sort(plate_levels)

  create_cols <- function(df) {
    dplyr::mutate(df,
                  plate = factor(plate, levels = plate_levels),
                  dup_lbl = ifelse(grepl("\\d", plate), paste("Plate", gsub("\\D", "", plate)), as.character(plate)),

                  base_label = if(n_groups > 1) {
                    if (n_distinct_constructs == 1) compound else construct_compound
                  } else {
                    dup_lbl
                  },
                  legend_group = factor(base_label)
    )
  }

  plot_curve <- create_cols(plot_curve)
  plot_stats <- create_cols(plot_stats)

  # ============================================================================
  # 6. GGPLOT CONSTRUCTION
  # ============================================================================

  if (is.null(plot_title)) {
    if (n_groups == 1) {
      final_title <- unique(plot_raw$compound)[1]
    } else {
      if (n_distinct_constructs == 1) {
        final_title <- unique(plot_raw$construct)[1]
      } else if (n_distinct_compounds == 1) {
        final_title <- unique(plot_raw$compound)[1]
      } else {
        final_title <- "Comparative DRC Analysis"
      }
    }
  } else {
    final_title <- plot_title
  }

  n_legend_items <- length(unique(plot_curve$legend_group))
  final_colors <- generate_palette(colors, n_legend_items)

  if (verbose && !is.null(colors) && length(colors) == 1) {
    message("Using palette: ", colors, " (", n_legend_items, " colors needed)")
  }

  use_shapes <- shape_by_duplicate && n_groups == 1
  if (use_shapes) {
    if (is.null(point_shapes)) final_shapes <- 15:(15 + n_legend_items - 1) else final_shapes <- rep(point_shapes, length.out = n_legend_items)
    names(final_shapes) <- levels(plot_curve$legend_group)
  }

  yt <- if (!is.null(y_axis_title)) {
    y_axis_title
  } else {
    if (any(selected_norm_flags, na.rm = TRUE)) {
      "Normalized BRET ratio [%]"
    } else {
      "BRET ratio"
    }
  }

  if(!is.null(legend_title)) {
    lt <- legend_title
  } else {
    if (n_groups > 1) {
      if (n_distinct_constructs == 1) {
        lt <- "Compound"
      } else if (n_distinct_compounds == 1) {
        lt <- "Construct"
      } else {
        lt <- "Condition"
      }
    } else {
      lt <- "Replicate"
    }
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = plot_curve,
                       ggplot2::aes(x = log_inhibitor, y = response, color = legend_group, group = unique_id),
                       linewidth = 1, alpha = 0.8) +
    ggplot2::geom_point(data = plot_stats,
                        ggplot2::aes(x = log_inhibitor, y = mean_response, color = legend_group,
                                     shape = if(use_shapes) legend_group else NULL),
                        size = 3)

  if (show_error_bars) {
    p <- p + ggplot2::geom_errorbar(data = plot_stats,
                                    ggplot2::aes(x = log_inhibitor, ymin = mean_response - sd_response, ymax = mean_response + sd_response, color = legend_group),
                                    width = 0.05, linewidth = 0.5, alpha = 0.6)
  }

  p <- p + ggplot2::scale_color_manual(values = final_colors, name = lt)

  if (use_shapes) p <- p + ggplot2::scale_shape_manual(values = final_shapes, name = lt)
  else p <- p + ggplot2::guides(shape = "none")

  if (!is.null(facet_by)) {
    if (!facet_by %in% names(plot_curve)) warning("Facet column '", facet_by, "' not found.")
    else p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_by)), ncol = facet_ncol, scales = facet_scales)
  }

  base_theme <- if(show_grid) ggplot2::theme_minimal(base_family = font_family) else ggplot2::theme_classic(base_family = font_family)

  p <- p +
    ggplot2::labs(title = final_title, x = expression(paste("Log"[10], " Concentration [M]")), y = yt) +
    base_theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.text = ggplot2::element_text(color = axis_text_color, size = axis_text_size),
      axis.title = ggplot2::element_text(color = axis_title_color, size = axis_title_size, face = "bold"),
      legend.position = if(show_legend) legend_position else "none",
      legend.title = ggplot2::element_text(face = "bold", size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      strip.background = ggplot2::element_rect(fill = "#f0f0f0", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 11)
    )

  if (!show_grid) p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (!is.null(y_limits) && is.numeric(y_limits) && length(y_limits) == 2) {
    if (is.null(facet_by) || !grepl("free", facet_scales)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    } else {
      if (verbose) message("Note: 'y_limits' ignored due to free facet scales.")
    }
  }

  # ============================================================================
  # 7. OUTPUT
  # ============================================================================
  print(p)

  if (!is.null(save_plot)) {
    ggplot2::ggsave(save_plot, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi)
    if(verbose) message("Plot saved to: ", save_plot)
  }

  return(invisible(p))
}
