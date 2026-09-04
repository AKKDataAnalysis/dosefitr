#' Plot Dose-Response Curves from \code{rout_outliers()} Output
#'
#' Generates a multi-panel figure showing fitted 3PL/4PL curves, replicate
#' data points, and ROUT outliers for every compound in a
#' \code{rout_outliers()} result.
#'
#' @param rout_output Return value of \code{\link{rout_outliers}}.
#'
#' @param title Optional character string for the overall plot title.
#'   \code{NULL} (default) omits the title.
#'   
#' @param subplot_title Character. Controls what text is used as the title of
#'   each compound sub-plot. One of \code{"full"} (default, e.g.
#'   \code{"KinaseA:Cpd1"}), \code{"compound"} (e.g. \code{"Cpd1"}), or
#'   \code{"construct"} (e.g. \code{"KinaseA"}).
#' @param panel_spacing Numeric. Spacing between sub-plots in the panel, in
#'   centimetres (default \code{0.5}). Increase for more breathing room between
#'   plots.
#'
#' @param ncol Integer.  Number of columns in the compound panel grid
#'   (default \code{4}).
#'
#' @param file Character string giving the output file path
#'   (e.g. \code{"curves.png"}).  If \code{NULL} (default), the combined
#'   \pkg{patchwork} object is returned invisibly without saving.
#'
#' @param width Panel width in inches. If \code{NULL} (default), computed as
#'   \code{ncol * plot_width}.
#'
#' @param height Panel height in inches. If \code{NULL} (default), computed as
#'   \code{ceiling(n_compounds / ncol) * plot_height + 0.6}.
#'
#' @param plot_width Width of each individual subplot in inches. Default: 3.2.
#'
#' @param plot_height Height of each individual subplot in inches. Default: 3.0.
#' @param label_sep Character separator used in display labels between
#'   construct and compound names.  Defaults to \code{":"}.  Change to
#'   e.g. \code{"/"} to show \code{"Kinase/Cpd1"} instead of
#'   \code{"Kinase:Cpd1"} in plot titles and legends.  The internal data
#'   always uses \code{":"}; this parameter only affects display.
#' @param panel_title Character string for the overall panel title. If
#'   \code{NULL} (default), uses \code{title}. Overrides \code{title} when
#'   both are provided.
#' @param axis_label_size Numeric. Font size for axis titles (default: 8).
#' @param axis_text_size Numeric. Font size for axis tick labels (default: 7).
#' @param plot_title_size Numeric. Font size for the overall plot title
#'   (default: 13).
#' @param subplot_title_size Numeric. Font size for per-compound subplot
#'   titles (default: 9).
#' @param subplot_subtitle_size Numeric. Font size for per-compound subplot
#'   subtitles (default: 7).
#' @param legend_text_size Numeric. Font size for legend text (default: 7).
#' @param caption_size Numeric. Font size for the plot caption (default: 8).
#' @param point_size Numeric. Size of data points (default: 2.2).
#' @param outlier_point_size Numeric. Size of outlier points (default: 3.5).
#' @param outlier_label_size Numeric. Font size for outlier labels (default: 2.6).
#' @param line_width Numeric. Line width for the fitted curve (default: 0.7).
#' @param axis_line_width Numeric. Line width for axis lines (default: 0.8).
#' @param legend_key_size Numeric. Size of legend keys in cm (default: 0.35).
#' @param dpi Integer. Resolution of saved PNG files (default: 150).
#' @param width Numeric. Width of the entire panel in inches. If \code{NULL}
#'   (default), computed as \code{ncol * plot_width}.
#' @param height Numeric. Height of the entire panel in inches. If \code{NULL}
#'   (default), computed as \code{ceiling(n_compounds / ncol) * plot_height + 0.6}.
#' @param plot_width Numeric. Width of each individual subplot in inches. Default: 3.2.
#' @param plot_height Numeric. Height of each individual subplot in inches. Default: 3.0.
#' @param show_ic50 Logical. Add vertical dashed line at IC50/EC50 (default: TRUE).
#' @param y_limits Numeric vector of length 2. Force y-axis limits (default: NULL,
#'   auto-computed from data).
#' @param x_limits Numeric vector of length 2. Force x-axis limits (default: NULL,
#'   auto-computed from data).
#' @param show_grid Logical. Add major grid lines (default: FALSE).
#' @param base_family Character. Font family for all text (default: "Liberation Sans").
#' @param outlier_alpha Numeric. Transparency of outlier points (default: 0.7).
#' @param curve_alpha Numeric. Transparency of fitted curve (default: 1.0).
#' @param show_n Logical. Show number of non-outlier replicates in subtitle
#'   (default: TRUE).
#' @param theme Character. Plot theme: "prism" (default, ggprism) or "bw"
#'   (classic black-and-white).
#'
#' @return Invisibly returns the combined \pkg{patchwork} ggplot object.
#'   Each panel shows:
#'   \itemize{
#'     \item Smooth fitted 3PL/4PL curve (grey line, 200 points).
#'     \item Replicates shown with distinct colours and shapes. Rep 1 uses
#'       blue circles, Rep 2 orange triangles, and Rep 3 green squares;
#'       additional replicates are assigned further palette entries
#'       automatically.
#'     \item ROUT outliers as red \eqn{\times} with standardised residual
#'       label (via \pkg{ggrepel}).
#'     \item Subtitle: model used, dynamic range \%, convergence warning
#'       if applicable.
#'   }
#'
#' @details
#' Requires \pkg{ggplot2}, \pkg{ggprism}, \pkg{ggrepel}, and
#' \pkg{patchwork}.  The function checks for these packages at call time
#' and stops with an informative message if any are missing.
#'
#' The x-axis uses \eqn{10^x} notation (e.g. \eqn{10^{-9}}) regardless
#' of whether \code{log_base = "log10"} or \code{"ln"} was used during
#' fitting.
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_poc")
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
#' rr <- rout_outliers(
#'   data      = ratio_res$plate_01$result$modified_ratio_table,
#'   Q         = 0.01,
#'   n_param   = 4L,
#'   direction = "inhibition",
#'   verbose   = FALSE
#' )
#'
#' p <- plot_outliers_curves(rout_output = rr, ncol = 4L)
#' inherits(p, "gg") || inherits(p, "patchwork")
#' }
#' @seealso \code{\link{rout_outliers}},
#'   \code{\link{plot_outliers_batch_curves}}
#'
#' @export

plot_outliers_curves <- function(rout_output,
                                 title         = NULL,
                                 subplot_title = "full",
                                 panel_spacing = 0.5,
                                 ncol          = 4L,
                                 file   = NULL,
                                 width  = NULL,
                                 height = NULL,
                                 plot_width     = NULL,
                                 plot_height    = NULL,
                                 label_sep     = ":",
                                 panel_title   = NULL,
                                 axis_label_size      = 8,
                                 axis_text_size       = 7,
                                 plot_title_size      = 13,
                                 subplot_title_size   = 9,
                                 subplot_subtitle_size = 7,
                                 legend_text_size     = 7,
                                 caption_size         = 8,
                                 point_size           = 2.2,
                                 outlier_point_size   = 3.5,
                                 outlier_label_size   = 2.6,
                                 line_width           = 0.7,
                                 axis_line_width      = 0.8,
                                 legend_key_size      = 0.35,
                                 dpi                  = 150,
                                 show_ic50            = TRUE,
                                 y_limits             = NULL,
                                 x_limits             = NULL,
                                 show_grid            = FALSE,
                                 base_family          = "Liberation Sans",
                                 outlier_alpha        = 0.7,
                                 curve_alpha          = 1.0,
                                 show_n               = TRUE,
                                 theme                = "prism") {
  
  subplot_title <- match.arg(subplot_title, c("full", "compound", "construct"))
  
  # Derive the sub-plot title from a raw compound string ("Construct:Compound").
  .make_subplot_label <- function(compound_string) {
    parts <- strsplit(compound_string, ":", fixed = TRUE)[[1L]]
    switch(subplot_title,
           full      = gsub(":", label_sep, compound_string, fixed = TRUE),
           compound  = if (length(parts) >= 2L) parts[[2L]] else compound_string,
           construct = if (length(parts) >= 2L) parts[[1L]] else compound_string
    )
  }
  
  .check_plot_packages <- function() {
    missing_pkgs <- character(0)
    for (pkg in c("ggplot2", "ggrepel", "patchwork", "ggprism")) {
      if (!requireNamespace(pkg, quietly = TRUE))
        missing_pkgs <- c(missing_pkgs, pkg)
    }
    if (length(missing_pkgs) > 0L)
      stop("The following packages are required for plotting but are not installed: ",
           paste(missing_pkgs, collapse = ", "),
           "\nInstall with: install.packages(c(",
           paste0('"', missing_pkgs, '"', collapse = ", "), "))",
           call. = FALSE)
    invisible(TRUE)
  }
  
  # Returns a 3PL model function with hill slope fixed to hill_fixed.
  # Signature matches OptimModel::hill_model but with only 3 free parameters.
  .make_hill_3p <- function(hill_fixed) {
    function(theta, x)
      OptimModel::hill_model(c(theta[1L], theta[2L], theta[3L], hill_fixed), x)
  }
  
  .check_plot_packages()
  
  res   <- rout_output$results
  Q_val <- rout_output$params$Q
  caption_txt <- sprintf("Red \u2717 = ROUT outlier%s. Label = standardised residual.",
                         if (!is.null(Q_val)) sprintf(" (Q=%.3f)", Q_val) else "")
  
  # Build replicate aesthetics dynamically. The previous implementation only
  # defined entries for replicates 1 and 2, causing ggplot2 to assign NA
  # colour/shape values to replicate 3 (and silently omit those points).
  # Keep the original colours/shapes for reps 1 and 2 for backward visual
  # compatibility, then use colour-blind-friendly entries for further reps.
  rep_ids <- as.character(sort(unique(res$replicate[!is.na(res$replicate)])))

  colour_pool <- c(
    "#0279EE", # Rep 1: blue
    "#FF9400", # Rep 2: orange
    "#009E73", # Rep 3: green
    "#CC79A7", # Rep 4: purple
    "#D55E00", # Rep 5: vermillion
    "#56B4E9", # Rep 6: light blue
    "#F0E442", # Rep 7: yellow
    "#000000"  # Rep 8: black
  )
  shape_pool <- c(16, 17, 15, 18, 8, 25, 24, 23)

  rep_colours <- stats::setNames(
    rep_len(colour_pool, length(rep_ids)),
    rep_ids
  )
  rep_shapes <- stats::setNames(
    rep_len(shape_pool, length(rep_ids)),
    rep_ids
  )
  rep_labels <- stats::setNames(paste("Rep", rep_ids), rep_ids)
  
  # Filter out compounds whose name is NA / "NA" / "NA_N" (exact match only;
  # names that merely contain "NA" as a substring are kept).
  .is_na_name_oc <- function(x) {
    if (is.null(x) || length(x) == 0L || is.na(x)) return(TRUE)
    toupper(sub("_\\d+$", "", trimws(x))) == "NA"
  }
  .cmpd_part <- function(s) {
    p <- strsplit(s, ":", fixed = TRUE)[[1L]]
    if (length(p) >= 2L) trimws(p[[2L]]) else trimws(s)
  }
  .cons_part <- function(s) {
    p <- strsplit(s, ":", fixed = TRUE)[[1L]]
    if (length(p) >= 2L) trimws(p[[1L]]) else NA_character_
  }
  
  all_compounds <- unique(res$compound)
  compounds <- all_compounds[!vapply(all_compounds, function(cmpd) {
    .is_na_name_oc(.cmpd_part(cmpd)) || .is_na_name_oc(.cons_part(cmpd))
  }, logical(1L))]
  
  nrow_grid  <- ceiling(length(compounds) / ncol)
  
  plot_list <- lapply(compounds, function(cmpd) {
    
    df <- res[res$compound == cmpd, ]
    
    # Smooth fitted curve (200 points).
    # Use the correct model function per compound:
    #   - 4PL: hill_model(c(bottom, top, log10_EC50, hill), x)
    #   - 3PL: .make_hill_3p(hill_slope)(c(bottom, top, log10_EC50), x)
    # par[3] is log10(EC50) directly  --  no log(10^x) round-trip needed.
    # Detect concentration column name from results ("log10_conc" or "ln_conc").
    conc_col_name <- intersect(c("log10_conc", "ln_conc"), names(df))[1L]
    # Safety check: ensure concentration column is numeric and has no NA
    conc_values <- df[[conc_col_name]]
    if (!is.numeric(conc_values) || any(!is.finite(conc_values))) {
      warning(sprintf("Skipping compound '%s': concentration column '%s' contains NA or non-numeric values", cmpd, conc_col_name))
      return(NULL)
    }
    x_smooth <- seq(min(conc_values), max(conc_values), length.out = 200L)
    # hill_model expects ln(EC50) for par[3]; $log10_EC50 is in log10 scale.
    # Convert: ln(EC50) = log10_EC50 * log(10)
    # Safety check: if any required parameter is NA or non-numeric, skip this compound
    required_params <- c("log10_EC50", "bottom", "top", "hill_slope")
    for (param in required_params) {
      if (!param %in% names(df) || is.null(df[[param]]) || length(df[[param]]) == 0 ||
          !is.numeric(df[[param]][1L]) || !is.finite(df[[param]][1L])) {
        warning(sprintf("Skipping compound '%s': %s is NA or non-numeric", cmpd, param))
        return(NULL)
      }
    }
    ln_ec50 <- df$log10_EC50[1L] * log(10)
    # x_smooth is in log10 or ln units; hill_model expects linear concentration.
    # Use the correct inverse transform based on log_base stored in params.
    log_base_used <- rout_output$params$log_base
    x_linear <- if (log_base_used == "log10") 10^x_smooth else exp(x_smooth)
    if (df$model_used[1L] == "4PL") {
      y_smooth <- OptimModel::hill_model(
        c(df$bottom[1L], df$top[1L], ln_ec50, df$hill_slope[1L]),
        x_linear)
    } else {
      f3 <- .make_hill_3p(df$hill_slope[1L])
      y_smooth <- f3(c(df$bottom[1L], df$top[1L], ln_ec50), x_linear)
    }
    curve_df <- data.frame(x_smooth = x_smooth, y = y_smooth)
    
    # Safety check: ensure model_used and dynamic_range_pct are valid
    if (is.null(df$model_used) || length(df$model_used) == 0 || is.na(df$model_used[1L])) {
      warning(sprintf("Skipping compound '%s': model_used is NA", cmpd))
      return(NULL)
    }
    if (is.null(df$dynamic_range_pct) || length(df$dynamic_range_pct) == 0 ||
        !is.numeric(df$dynamic_range_pct[1L]) || !is.finite(df$dynamic_range_pct[1L])) {
      warning(sprintf("Skipping compound '%s': dynamic_range_pct is NA or non-numeric", cmpd))
      return(NULL)
    }

    conv_label  <- if (!all(df$converged)) " \u26a0 no conv." else ""
    n_label     <- if (show_n) sprintf(" | n=%d", sum(!df$outlier_fdr)) else ""
    subtitle    <- sprintf("%s | DR: %.0f%%%s%s",
                           df$model_used[1L], df$dynamic_range_pct[1L], conv_label, n_label)

    # Safety check: ensure bret_ratio is numeric
    if (!is.numeric(df$bret_ratio) || all(!is.finite(df$bret_ratio))) {
      warning(sprintf("Skipping compound '%s': bret_ratio contains only NA or non-numeric values", cmpd))
      return(NULL)
    }

    y_all  <- c(df$bret_ratio, y_smooth)
    y_pad  <- diff(range(y_all, na.rm = TRUE)) * 0.12
    y_lims <- if (!is.null(y_limits)) y_limits else c(min(y_all, na.rm = TRUE) - y_pad, max(y_all, na.rm = TRUE) + y_pad)
    x_lims <- if (!is.null(x_limits)) x_limits else range(df[[conc_col_name]], na.rm = TRUE)
    
    # Rename curve_df x column to match the concentration column name in df
    # so aes() references are consistent regardless of log_base setting.
    names(curve_df)[1L] <- conc_col_name
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data[[conc_col_name]], y = y),
        colour = "grey40", linewidth = line_width, alpha = curve_alpha) +
      {if (show_ic50 && !is.null(df$log10_EC50) && length(df$log10_EC50) > 0 && is.numeric(df$log10_EC50[1L]) && is.finite(df$log10_EC50[1L]))
        ggplot2::geom_vline(
          xintercept = df$log10_EC50[1L],
          linetype = "dashed", colour = "grey50", linewidth = 0.5)
      } +
      ggplot2::geom_point(
        data = df[!df$outlier_fdr, ],
        ggplot2::aes(x = .data[[conc_col_name]], y = bret_ratio,
                     colour = as.character(replicate),
                     shape  = as.character(replicate)),
        size = point_size, stroke = 0.6) +
      ggplot2::geom_point(
        data = df[df$outlier_fdr, ],
        ggplot2::aes(x = .data[[conc_col_name]], y = bret_ratio),
        colour = "#D62728", shape = 4, size = outlier_point_size, stroke = 1.4,
        alpha = outlier_alpha) +
      {if (any(df$outlier_fdr))
        ggrepel::geom_text_repel(
          data = df[df$outlier_fdr, ],
          ggplot2::aes(x = .data[[conc_col_name]], y = bret_ratio,
                       label = sprintf("%.1f SD", abs(std_residual))),
          colour = "#D62728", size = outlier_label_size, fontface = "bold",
          box.padding = 0.4, point.padding = 0.3,
          min.segment.length = 0.2, max.overlaps = 20)
      } +
      ggplot2::scale_colour_manual(
        values = rep_colours,
        breaks = rep_ids,
        labels = rep_labels,
        name = NULL) +
      ggplot2::scale_shape_manual(
        values = rep_shapes,
        breaks = rep_ids,
        labels = rep_labels,
        name = NULL) +
      ggplot2::scale_x_continuous(
        breaks = pretty(x_lims, n = 5),
        labels = function(x) parse(text = paste0("10^{", x, "}")),
        expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = x_lims, ylim = y_lims, clip = "on") +
      ggplot2::labs(
        title    = .make_subplot_label(cmpd),
        subtitle = subtitle,
        x        = "Concentration (M)",
        y        = "BRET ratio") +
      {if (theme == "prism") ggprism::theme_prism(base_size = axis_text_size, base_family = base_family)
        else ggplot2::theme_bw(base_size = axis_text_size, base_family = base_family)} +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = subplot_title_size,  face = "bold",  hjust = 0.5),
        plot.subtitle   = ggplot2::element_text(size = subplot_subtitle_size,  colour = "grey50", hjust = 0.5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(legend_key_size, "cm"),
        legend.text     = ggplot2::element_text(size = legend_text_size),
        panel.grid.major = if (show_grid) ggplot2::element_line(colour = "grey90", linewidth = 0.3) else ggplot2::element_blank(),
        axis.title      = ggplot2::element_text(size = axis_label_size),
        axis.text       = ggplot2::element_text(size = axis_text_size),
        axis.line       = ggplot2::element_blank(),
        panel.border    = ggplot2::element_blank(),
        plot.margin     = ggplot2::margin(t = 10, r = 8, b = 4, l = 6, unit = "pt"))
    
    x_range_oc <- range(df[[conc_col_name]], na.rm = TRUE)
    axis_segs_oc <- data.frame(
      x    = c(x_range_oc[1],  x_range_oc[1]),
      xend = c(x_range_oc[2],  x_range_oc[1]),
      y    = c(y_lims[1],      y_lims[1]),
      yend = c(y_lims[1],      y_lims[2])
    )
    p <- p +
      ggplot2::geom_segment(
        data = axis_segs_oc,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        colour = "black", linewidth = axis_line_width,
        inherit.aes = FALSE)
    p
  })

  # Remove NULL elements (skipped compounds)
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1L))]

  if (length(plot_list) == 0L) {
    warning("No compounds could be plotted (all had NA or non-numeric log10_EC50). Returning empty plot.")
    # Return an empty plot with a message
    p_empty <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "No valid compounds to plot\n(all had NA or non-numeric log10_EC50)",
                        size = 5, colour = "grey50") +
      ggplot2::theme_void() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)

    if (!is.null(file)) {
      subplot_width  <- plot_width  %||% 3.2
      subplot_height <- plot_height %||% 3.0
      final_width    <- width  %||% (ncol * subplot_width)
      final_height   <- height %||% (nrow_grid * subplot_height + 0.6)
      ggplot2::ggsave(file, p_empty, width = final_width, height = final_height,
                      dpi = dpi, bg = "white")
      message(sprintf("Saved: %s", file))
    }

    return(invisible(p_empty))
  }

  # Resolve panel title: explicit panel_title > title > NULL
  final_title <- panel_title %||% title

  combined <- (patchwork::wrap_plots(plot_list, ncol = ncol) &
                 ggplot2::theme(plot.margin = ggplot2::margin(
                   t = panel_spacing * 0.5, r = panel_spacing * 0.5,
                   b = panel_spacing * 0.5, l = panel_spacing * 0.5,
                   unit = "cm"))) +
    patchwork::plot_annotation(
      title   = final_title,
      caption = caption_txt,
      theme   = ggplot2::theme(
        plot.title   = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5),
        plot.caption = ggplot2::element_text(size = caption_size,  colour = "grey50", hjust = 0)))
  
  # Resolve panel dimensions: explicit width/height > plot_width/plot_height * ncol/nrow > defaults
  # plot_width/plot_height control individual subplot size; width/height control the whole panel
  subplot_width  <- plot_width  %||% 3.2
  subplot_height <- plot_height %||% 3.0
  final_width    <- width  %||% (ncol * subplot_width)
  final_height   <- height %||% (nrow_grid * subplot_height + 0.6)

  if (!is.null(file)) {
    ggplot2::ggsave(file, combined, width = final_width, height = final_height,
                    dpi = dpi, bg = "white")
    message(sprintf("Saved: %s", file))
  }
  
  invisible(combined)
}
