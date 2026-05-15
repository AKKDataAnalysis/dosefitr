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
#' @param ncol Integer.  Number of columns in the compound panel grid
#'   (default \code{4}).
#'
#' @param file Character string giving the output file path
#'   (e.g. \code{"curves.png"}).  If \code{NULL} (default), the combined
#'   \pkg{patchwork} object is returned invisibly without saving.
#'
#' @param width Plot width in inches.  Default: \code{ncol * 3.2}.
#'
#' @param height Plot height in inches.  Default:
#'   \code{ceiling(n_compounds / ncol) * 3.0 + 0.6}.
#'
#' @return Invisibly returns the combined \pkg{patchwork} ggplot object.
#'   Each panel shows:
#'   \itemize{
#'     \item Smooth fitted 3PL/4PL curve (grey line, 200 points).
#'     \item Rep 1 (blue circles) and Rep 2 (orange triangles).
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
#' \dontrun{
#' out <- rout_outliers(dat, Q = 0.01)
#'
#' # Display in RStudio viewer
#' plot_outliers_curves(out, title = "Plate 01")
#'
#' # Save to PNG
#' plot_outliers_curves(out, title = "Plate 01", file = "plate_01_curves.png")
#' }
#'
#' @seealso \code{\link{rout_outliers}},
#'   \code{\link{plot_outliers_batch_curves}}
#'
#' @export

plot_outliers_curves <- function(rout_output,
                                 title  = NULL,
                                 ncol   = 4L,
                                 file   = NULL,
                                 width  = NULL,
                                 height = NULL) {

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

  # Colour-blind friendly: rep1 = blue, rep2 = orange
  rep_colours <- c("1" = "#0279EE", "2" = "#FF9400")

  compounds  <- unique(res$compound)
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
    x_smooth <- seq(min(df[[conc_col_name]]), max(df[[conc_col_name]]), length.out = 200L)
    # hill_model expects ln(EC50) for par[3]; $log10_EC50 is in log10 scale.
    # Convert: ln(EC50) = log10_EC50 * log(10)
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

    conv_label  <- if (!all(df$converged)) " \u26a0 no conv." else ""
    subtitle    <- sprintf("%s | DR: %.0f%%%s",
                           df$model_used[1L], df$dynamic_range_pct[1L], conv_label)

    y_all  <- c(df$bret_ratio, y_smooth)
    y_pad  <- diff(range(y_all, na.rm = TRUE)) * 0.12
    y_lims <- c(min(y_all, na.rm = TRUE) - y_pad, max(y_all, na.rm = TRUE) + y_pad)

    # Rename curve_df x column to match the concentration column name in df
    # so aes() references are consistent regardless of log_base setting.
    names(curve_df)[1L] <- conc_col_name

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data[[conc_col_name]], y = y),
        colour = "grey40", linewidth = 0.7) +
      ggplot2::geom_point(
        data = df[!df$outlier_fdr, ],
        ggplot2::aes(x = .data[[conc_col_name]], y = bret_ratio,
                     colour = as.character(replicate),
                     shape  = as.character(replicate)),
        size = 2.2, stroke = 0.6) +
      ggplot2::geom_point(
        data = df[df$outlier_fdr, ],
        ggplot2::aes(x = .data[[conc_col_name]], y = bret_ratio),
        colour = "#D62728", shape = 4, size = 3.5, stroke = 1.4) +
      {if (any(df$outlier_fdr))
        ggrepel::geom_text_repel(
          data = df[df$outlier_fdr, ],
          ggplot2::aes(x = .data[[conc_col_name]], y = bret_ratio,
                       label = sprintf("%.1f SD", abs(std_residual))),
          colour = "#D62728", size = 2.6, fontface = "bold",
          box.padding = 0.4, point.padding = 0.3,
          min.segment.length = 0.2, max.overlaps = 20)
      } +
      ggplot2::scale_colour_manual(
        values = rep_colours,
        labels = c("1" = "Rep 1", "2" = "Rep 2"), name = NULL) +
      ggplot2::scale_shape_manual(
        values = c("1" = 16, "2" = 17),
        labels = c("1" = "Rep 1", "2" = "Rep 2"), name = NULL) +
      ggplot2::scale_x_continuous(
        breaks = pretty(range(df[[conc_col_name]]), n = 5),
        labels = function(x) parse(text = paste0("10^{", x, "}"))) +
      ggplot2::coord_cartesian(ylim = y_lims) +
      ggplot2::labs(
        title    = cmpd,
        subtitle = subtitle,
        x        = "Concentration (M)",
        y        = "BRET ratio") +
      ggprism::theme_prism(base_size = 9) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 9,  face = "bold",  hjust = 0.5),
        plot.subtitle   = ggplot2::element_text(size = 7,  colour = "grey50", hjust = 0.5),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(0.35, "cm"),
        legend.text     = ggplot2::element_text(size = 7),
        axis.title      = ggplot2::element_text(size = 8),
        axis.text       = ggplot2::element_text(size = 7),
        plot.margin     = ggplot2::margin(6, 8, 4, 6))
    p
  })

  combined <- patchwork::wrap_plots(plot_list, ncol = ncol) +
    patchwork::plot_annotation(
      title   = title,
      caption = caption_txt,
      theme   = ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
        plot.caption = ggplot2::element_text(size = 8,  colour = "grey50", hjust = 0)))

  if (is.null(width))  width  <- ncol * 3.2
  if (is.null(height)) height <- nrow_grid * 3.0 + 0.6

  if (!is.null(file)) {
    ggplot2::ggsave(file, combined, width = width, height = height,
                    dpi = 150, bg = "white")
    message(sprintf("Saved: %s", file))
  }

  invisible(combined)
}
