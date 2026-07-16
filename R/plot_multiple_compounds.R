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
#' @param legend_position Position of the legend: one of `"right"`, `"left"`,
#'   `"top"`, `"bottom"`, or `"none"`. Default: `"right"`.
#' @param show_grid Logical indicating whether to show background grid lines
#' @param transparent_background Logical.  If \code{TRUE}, the plot, panel,
#'   and legend backgrounds are set to fully transparent (\code{element_rect(fill = NA)})
#'   instead of white. Default: \code{FALSE}.
#' @param show_legend Logical indicating whether to display the legend
#' @param save_plot Defines whether to save the plot: \code{NULL} (do not
#'   save, default), \code{TRUE} (automatically saves with an auto-generated
#'   name in the current working directory), or a file path.  When a file
#'   path is supplied and it already carries an extension (\code{.png},
#'   \code{.pdf}, \code{.jpeg}, \code{.tiff}, \code{.svg}, \code{.eps},
#'   ...), that extension wins and \code{format} is ignored (a message is
#'   emitted when \code{verbose = TRUE}).  When a file path is supplied
#'   without an extension, or when \code{save_plot = TRUE}, the extension
#'   is taken from \code{format}.
#' @param format Character string giving the file format used when saving
#'   the plot.  Passed through as the file extension and forwarded to
#'   \code{ggplot2::ggsave()}, which infers the graphics device from it.
#'   Common values are \code{"png"} (default), \code{"pdf"}, \code{"svg"},
#'   \code{"jpeg"}, \code{"tiff"}, \code{"eps"}, and \code{"bmp"}; any
#'   extension \code{ggsave} recognises will work.  Only takes effect
#'   when the plot is actually saved (\code{save_plot} not \code{NULL}),
#'   and only when \code{save_plot} does not already specify an extension
#'   (see \code{save_plot} for the precedence rule).  Mirrors the
#'   \code{format} argument of \code{\link{batch_save_all_drc_plots}}.
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
#' @param y_axis_title Character string for the y-axis label.  If \code{NULL}
#'   (default), auto-detected from the batch result: \code{"Cell Viability (\%)"}
#'   or \code{"Luminescence"} for viability assays, \code{"Normalized BRET ratio [\%]"}
#'   or \code{"BRET ratio"} for NanoBRET assays (depending on whether
#'   \code{normalize} was \code{TRUE} or \code{FALSE}).
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
#' @param legend_width Numeric, \code{"auto"}, or NULL. Target width (in cm) for the legend
#'   column, including the legend keys, labels, and any internal margins.
#'   When specified, the legend column is padded via
#'   \code{legend.box.margin} so it occupies exactly this width — guaranteeing
#'   that the data panel is identically sized across plots with different
#'   label lengths. Use \code{legend_width = "auto"} to measure and return
#'   the current legend width without padding (useful for determining the
#'   target width across a set of plots). Only affects right- and
#'   left-positioned legends. \code{NULL} (default) disables padding.
#' @param legend_title Title for the legend (displayed above symbols).
#' @param verbose Logical; if TRUE (default), prints informative messages about processing
#'   steps, compound matches, and color assignments.
#' @param point_size Numeric value for the point size
#' @param axis_line_width Numeric. Line width of the manually drawn x/y axis
#'   lines (default: 0.8).
#' @param axis_vjust Numeric or NULL. Vertical justification (\code{vjust}) of
#'   the axis titles. \code{NULL} (default) leaves the ggplot2 default unchanged.
#' @param tick_length Numeric. Axis tick length in centimetres (default: 0.15).
#' @param error_linewidth Numeric. Line width of the error bars (default: 0.5).
#' @param point_alpha Numeric between 0 and 1. Opacity of the data points
#'   (default: 1, fully opaque). Does not affect error-bar opacity.
#' @param legend_spacing Numeric or NULL. Spacing (in points) between legend
#'   items, applied via \code{legend.spacing}. \code{NULL} (default) leaves the
#'   theme default unchanged.
#' @param aspect_ratio Numeric or NULL. Panel aspect ratio (\code{aspect.ratio}
#'   in \code{theme()}). \code{NULL} (default) leaves it unset.
#' @param byrow Logical. Whether legend keys are filled by row (\code{TRUE}) or
#'   by column (\code{FALSE}, default). Applied via \code{theme(legend.byrow=)}
#'   (the modern ggplot2 mechanism; the old \code{guide_legend(byrow=)} argument
#'   was removed in ggplot2 4.0).
#' @param label_sep Character string. Separator used for DISPLAY purposes in
#'   titles, legends, and filenames. When \code{NULL} (default), auto-detected
#'   from \code{attr(results, "label_sep")}; falls back to \code{":"}.
#'   The data separator (used for parsing compound names) is always read from
#'   the attribute and is never overridden by this argument. This separation
#'   ensures that changing the display separator does not break title or legend
#'   logic, which was a bug in the previous single-separator design.
#' @param legend_label Character string controlling the text shown in the
#'   side-legend entries.  One of:
#'   \describe{
#'     \item{\code{"auto"} (default)}{Smart shortening: if every entry
#'       shares the same construct (target), show only the compound part;
#'       if every entry shares the same compound, show only the construct
#'       part; otherwise show the full \code{construct:compound} label.}
#'     \item{\code{"compound"}}{Always show only the compound part,
#'       regardless of whether constructs differ.}
#'     \item{\code{"construct"}}{Always show only the construct (target)
#'       part, regardless of whether compounds differ.}
#'     \item{\code{"full"}}{Always show the full
#'       \code{construct:compound} label, even when auto would
#'       shorten it.}
#'   }
#'   Labels that do not contain the display separator (bare compound
#'   names) are passed through unchanged in every mode.  Line-wrapping
#'   via \code{legend_label_wrap} is applied after this choice.
#' @param plot_margin Margin or NULL. Plot margin applied via
#'   \code{theme(plot.margin = )}. Accepts a \code{ggplot2::margin()} object.
#'   \code{NULL} (default) uses the built-in margin
#'   (t = 12, r = 8, b = 8, l = 8 pt).
#' @param ic50_linetype Character or integer. Line type for the IC50 vertical
#'   reference lines (default: \code{"dashed"}). See \code{?par} for valid values.
#' @param ic50_linewidth Numeric. Line width of the IC50 vertical reference lines
#'   (default: 0.5).
#' @param ic50_line_alpha Numeric between 0 and 1. Opacity of the IC50 vertical
#'   reference lines (default: 0.8).
#' @param error_alpha Numeric between 0 and 1. Opacity of the error bars
#'   (default: 0.6).
#' @param point_shape Integer or character. Default shape of data points when
#'   shape mapping is not active (default: 16, filled circle). See \code{?points}
#'   for valid values. Ignored when \code{point_shapes} is provided.
#' @param grid_color Character string. Colour of the major grid lines when
#'   \code{show_grid = TRUE} (default: \code{"grey90"}).
#' @param grid_minor_color Character string. Colour of the minor grid lines when
#'   \code{show_grid = TRUE} (default: \code{"grey95"}).
#' @param grid_linewidth Numeric. Line width of the major grid lines
#'   (default: 0.5).
#' @param axis_expand Numeric vector of length 2. Expansion constants for the
#'   continuous axis scales, passed to \code{scale_*_continuous(expand = )}.
#'   Default \code{c(0, 0)} removes all padding between data and axis edge.
#' @param legend_key_fill Character string. Fill colour of the legend key
#'   backgrounds (default: \code{"white"}).
#' @param legend_bg_fill Character string. Fill colour of the legend background
#'   (default: \code{"white"}).
#' @param title_hjust Numeric. Horizontal justification of the plot title
#'   (default: 0.5, centred). Use 0 for left-aligned, 1 for right-aligned.
#' @param legend_key_linewidth Numeric. Line width of the fitted curve shown
#'   in legend keys via \code{guide_legend(override.aes)} (default: 0.8).
#' @param legend_lineheight Numeric. Line height multiplier for legend text
#'   (default: 0.8). Controls spacing between wrapped lines.
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
#' @return A ggplot2 object with the generated plot, containing metadata attributes.
#'   Metadata includes \code{legend_width_cm} (measured legend width in cm).
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
#'   \item \strong{X-axis limits}: Automatically calculated from data with 5\% margin
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
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_pmc")
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
#' p <- plot_multiple_compounds(
#'   results          = drc_res,
#'   compound_indices = 1:3,
#'   plate            = "plate_01",
#'   verbose          = FALSE
#' )
#' inherits(p, "ggplot")
#' }
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
                                    show_grid = FALSE,
                                    show_border = FALSE,
                                    transparent_background = FALSE,
                                    save_plot = NULL,
                                    format = "png",
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
                                    axis_line_width = 0.8,
                                    axis_vjust = NULL,
                                    tick_length = 0.15,
                                    error_linewidth = 0.5,
                                    point_alpha = 1,
                                    legend_spacing = NULL,
                                    aspect_ratio = NULL,
                                    byrow = FALSE,
                                    label_sep = NULL,
                                    legend_label = "auto",
                                    legend_width = NULL,
                                    ic50_linetype = "dashed",
                                    ic50_linewidth = 0.5,
                                    ic50_line_alpha = 0.8,
                                    error_alpha = 0.6,
                                    point_shape = 16,
                                    grid_color = "grey90",
                                    grid_minor_color = "grey95",
                                    grid_linewidth = 0.5,
                                    axis_expand = c(0, 0),
                                    legend_key_fill = "white",
                                    legend_bg_fill = "white",
                                    title_hjust = 0.5,
                                    legend_key_linewidth = 0.8,
                                    legend_lineheight = 0.8,
                                    plot_margin = NULL,
                                    plate = NULL) {


  # Null-coalescing operator
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

  # Resolve data separator (for parsing) and display separator (for titles/legends)
  # data_sep: the separator used in the actual compound names stored in the data.
  # label_sep: the separator used for DISPLAY in titles, legends, and filenames.
  # Resolution priority: explicit label_sep arg > attr(results, "label_sep") > ":"
  if (is.null(label_sep)) {
    label_sep <- attr(results, "label_sep")
    if (is.null(label_sep) || !is.character(label_sep) ||
        length(label_sep) != 1L || is.na(label_sep) || nchar(label_sep) == 0L) {
      label_sep <- ":"
    }
  }
  data_sep <- attr(results, "label_sep") %||% ":"

  # Validate legend_label. Kept close to label_sep since they share
  # responsibility for how construct/compound is displayed in legends.
  LEGEND_LABEL_MODES <- c("auto", "compound", "construct", "full")
  if (length(legend_label) != 1L || !is.character(legend_label) ||
      is.na(legend_label) || !legend_label %in% LEGEND_LABEL_MODES) {
    stop("`legend_label` must be one of: ",
         paste0("\"", LEGEND_LABEL_MODES, "\"", collapse = ", "),
         ". Got: ", paste(deparse(legend_label), collapse = " "),
         call. = FALSE)
  }

  # Validate `format`.  Mirrors batch_save_all_drc_plots(): no whitelist,
  # any extension ggsave() knows will work.  We only enforce that it is
  # a single, non-empty character string with no leading dot / whitespace
  # so paste0(name, ".", format) always yields a sane filename.
  if (length(format) != 1L || !is.character(format) || is.na(format) ||
      !nzchar(trimws(format))) {
    stop("`format` must be a single non-empty character string ",
         "(e.g. \"png\", \"pdf\", \"svg\", \"jpeg\", \"tiff\"). Got: ",
         paste(deparse(format), collapse = " "),
         call. = FALSE)
  }
  format <- sub("^\\.", "", trimws(format))  # tolerate ".png" and "  png "

  # ============================================================================
  # 0. PLATE EXTRACTION (when batch_drc_analysis result is supplied)
  # ============================================================================
  # Detect whether `results` is a batch_drc_analysis output by checking for
  # the $drc_results slot. If so, extract the requested plate's drc_result.
  # When `plate` is NULL, the first plate is used automatically.
  if (!is.null(results$drc_results)) {
    # Capture batch-level metadata before overwriting results with the
    # single-plate drc_result (which has no $metadata slot).
    batch_metadata <- results$metadata
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
    # Re-attach metadata so downstream code can read assay_type and normalize.
    results$metadata <- batch_metadata
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
    if (grepl(data_sep, input, fixed = TRUE)) {
      parts <- strsplit(input, data_sep, fixed = TRUE)[[1]]
      if (length(parts) >= 2 && nchar(parts[1]) > 0 && nchar(parts[length(parts)]) > 0) {
        return("both")
      }
    }
    return("unknown")
  }

  find_matches <- function(pattern, compound_list, exact_first = TRUE) {
    all_names <- sapply(compound_list, function(x) {
      name <- strsplit(x$full_name, " \\| ")[[1]][1]
      return(name)
    })

    searchable <- list()
    for (i in seq_along(all_names)) {
      name <- all_names[i]

      if (grepl(data_sep, name, fixed = TRUE)) {
        parts <- strsplit(name, data_sep, fixed = TRUE)[[1]]
        target <- parts[1]
        compound <- ifelse(length(parts) >= 2, paste(parts[-1], collapse = data_sep), name)

        compound_clean <- trimws(gsub("\\.2$", "", compound))
        target_clean <- trimws(gsub("\\.2$", "", target))

        searchable[[i]] <- list(
          full = name,
          target = target_clean,
          compound = compound_clean,
          target_compound = paste(target_clean, compound_clean, sep = data_sep),
          display_name = if (target_clean == compound_clean) compound_name else paste(target_clean, compound_clean, sep = label_sep)
        )
      } else {
        searchable[[i]] <- list(
          full = name,
          target = name,
          compound = name,
          target_compound = name,
          display_name = name
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

    if (grepl(data_sep, pattern, fixed = TRUE)) {
      pattern_parts <- strsplit(pattern, data_sep, fixed = TRUE)[[1]]
      pattern_target <- pattern_parts[1]
      pattern_compound <- ifelse(length(pattern_parts) >= 2, paste(pattern_parts[-1], collapse = data_sep), "")

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
    # Wrap a single word that contains hyphens/underscores by splitting on them.
    # Used as a fallback for individual words that are still too long after
    # space-based wrapping. Returns the word unchanged if no good split exists.
    wrap_hyphenated_word <- function(word, w) {
      if (nchar(word) <= w) return(word)
      if (!grepl("[-_]", word)) return(word)
      parts <- strsplit(word, "[-_]")[[1]]
      if (length(parts) < 2) return(word)
      current_length <- 0
      break_points <- numeric(0)
      for (j in seq_along(parts)) {
        current_length <- current_length + nchar(parts[j]) + 1
        if (current_length <= w + 1) break_points <- c(break_points, j)
      }
      if (length(break_points) > 0 && max(break_points) < length(parts)) {
        bp <- max(break_points)
        line1 <- paste(parts[1:bp], collapse = "-")
        line2 <- paste(parts[(bp + 1):length(parts)], collapse = "-")
        return(paste(line1, line2, sep = "\n"))
      }
      word
    }

    sapply(labels, function(label) {
      if (is.na(label) || nchar(label) <= width) return(label)

      # Strategy 1: space-based wrapping (handles "PF-05236216 HYDROCHLORIDE"
      # correctly by treating "PF-05236216" as one token).
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

      # Strategy 2: no spaces — try hyphen/underscore split within the single token.
      if (grepl("[-_]", label)) {
        result <- wrap_hyphenated_word(label, width)
        if (result != label) return(result)
      }

      # Strategy 3: hard break at midpoint (last resort for long single tokens).
      if (nchar(label) > width) {
        break_point <- ceiling(nchar(label) / 2)
        part1 <- substr(label, 1, break_point)
        part2 <- substr(label, break_point + 1, nchar(label))
        return(paste(part1, part2, sep = "\n"))
      }

      return(label)
    }, USE.NAMES = FALSE)
  }


  # Measure the actual rendered width (in cm) of the legend in a ggplot object.
  # Returns 0 if no legend is present.
  measure_legend_width <- function(p) {
    g <- ggplot2::ggplotGrob(p)
    guide_idx <- grep("guide-box", g$layout$name)
    if (length(guide_idx) == 0) return(0)
    for (i in guide_idx) {
      l_col <- g$layout$l[i]
      r_col <- g$layout$r[i]
      widths_cm <- sapply(l_col:r_col, function(col) {
        w <- g$widths[col]
        tryCatch(as.numeric(grid::convertUnit(w, "cm")), error = function(e) NA)
      })
      total <- sum(widths_cm, na.rm = TRUE)
      if (total > 0.01) return(total)
    }
    return(0)
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

      if (grepl(data_sep, compound_name, fixed = TRUE)) {
        name_parts <- strsplit(compound_name, data_sep, fixed = TRUE)[[1]]
        target <- name_parts[1]
        compound <- ifelse(length(name_parts) >= 2, paste(name_parts[-1], collapse = data_sep), target)
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
        target_compound = paste(target_clean, compound_clean, sep = data_sep),
        display_name = if (target_clean == compound_clean) compound_name else paste(target_clean, compound_clean, sep = label_sep),
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
      message("Selected by position ", position, ": ", compound_list[[position]]$display_name)
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
          message("  Match ", i, ": ", selected_compounds[[i]]$display_name)
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
        warning("Compound ", comp_info$display_name, ": insufficient data - skipping")
        next
      }

      x_range <- range(valid_data$log_inhibitor, na.rm = TRUE)
      x_seq <- seq(x_range[1], x_range[2], length.out = 100)

      curve_df <- data.frame(
        log_inhibitor = x_seq,
        response = predict(result$model, newdata = data.frame(log_inhibitor = x_seq)),
        compound = comp_info$display_name,
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
            compound = comp_info$display_name,
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

  # Auto-detect y-axis title from assay type when not supplied by the user.
  # Reads assay_type and normalize from the batch_drc_analysis() metadata if
  # available; falls back to the legacy normalization-only heuristic otherwise.
  y_axis_title_final <- if (!is.null(y_axis_title)) {
    y_axis_title
  } else {
    assay_src  <- results$metadata$assay_type
    normalized <- if (!is.null(results$metadata$normalize)) results$metadata$normalize else is_normalized
    if (!is.null(assay_src) && assay_src == "viability") {
      if (normalized) "Cell Viability (%)" else "Luminescence"
    } else {
      # nanobret or unknown — legacy behaviour
      if (normalized) "Normalized BRET ratio [%]" else "BRET ratio"
    }
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
        if (grepl(data_sep, input_used, fixed = TRUE)) {
          parts <- strsplit(input_used, data_sep, fixed = TRUE)[[1]]
          if (length(parts) >= 2) {
            title <- paste(parts[-1], collapse = data_sep)
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
      if (grepl(label_sep, label, fixed = TRUE)) {
        parts <- strsplit(label, label_sep, fixed = TRUE)[[1]]
        target_part <- parts[1]
        compound_part <- if (length(parts) >= 2) paste(parts[-1], collapse = label_sep) else parts[1]
        return(list(target = target_part, compound = compound_part))
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

  # If the user picked an explicit mode, override the smart shortening.
  # Labels lacking `label_sep` are passed through unchanged (bare compound
  # names never had a construct half to strip).
  if (legend_label != "auto") {
    split_labels <- function(labels) {
      # Returns a list of length(labels), each element list(target=, compound=)
      lapply(labels, function(lbl) {
        if (grepl(label_sep, lbl, fixed = TRUE)) {
          parts <- strsplit(lbl, label_sep, fixed = TRUE)[[1L]]
          list(
            target   = parts[[1L]],
            compound = if (length(parts) >= 2L) paste(parts[-1L], collapse = label_sep) else parts[[1L]]
          )
        } else {
          list(target = lbl, compound = lbl)
        }
      })
    }
    parts_list <- split_labels(compound_labels)

    smart_legend_names <- switch(legend_label,
      compound  = vapply(parts_list, `[[`, character(1L), "compound"),
      construct = vapply(parts_list, `[[`, character(1L), "target"),
      full      = compound_labels
    )
    if (isTRUE(verbose)) {
      message("Legend strategy: legend_label = '", legend_label,
              "' -> overriding auto shortening")
    }
  }

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
  # point_shapes = TRUE  → default optimal shapes, one per compound
  # point_shapes = NULL  → all points use shape 16 (filled circle), no per-compound mapping
  # point_shapes = <vec> → custom shapes, recycled to cover all compounds
  optimal_shapes <- c(16, 17, 15, 18, 8, 1, 2, 0, 5, 6, 7, 10, 11, 12, 13, 14)

  if (isTRUE(point_shapes)) {
    use_shape_mapping <- TRUE
    point_shapes <- if (n_valid_compounds <= length(optimal_shapes)) {
      optimal_shapes[1:n_valid_compounds]
    } else {
      rep(optimal_shapes, length.out = n_valid_compounds)
    }
  } else if (is.null(point_shapes)) {
    use_shape_mapping <- FALSE
    point_shapes <- rep(16L, n_valid_compounds)
  } else {
    use_shape_mapping <- TRUE
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
  #   "log10"  — already log10 molar (default, no conversion)
  #   "molar"  — raw molar (e.g. 1e-9); converted via log10()
  #   "uM"     — micromolar (e.g. 0.001 to 25); converted via log10(x * 1e-6)
  #   "nM"     — nanomolar; converted via log10(x * 1e-9)
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
    ggplot2::scale_y_continuous(expand = axis_expand) +
    ggplot2::scale_x_continuous(expand = axis_expand) +
    ggplot2::coord_cartesian(xlim = x_limits_final, ylim = y_limits_final, clip = "on")

  p <- p +
    ggplot2::geom_line(data = plot_data$curves,
                       ggplot2::aes(x = log_inhibitor, y = response,
                                    group = compound, color = compound),
                       linewidth = curve_linewidth, alpha = curve_alpha)

  if (nrow(plot_data$points) > 0) {
    if (use_shape_mapping) {
      p <- p +
        ggplot2::geom_point(data = plot_data$points,
                            ggplot2::aes(x = log_inhibitor, y = mean_response,
                                         shape = compound, color = compound),
                            size = point_size, alpha = point_alpha)

      p <- p + ggplot2::scale_shape_manual(
        values = setNames(point_shapes[1:n_valid_compounds], compound_labels),
        guide = "none"
      )
    } else {
      p <- p +
        ggplot2::geom_point(data = plot_data$points,
                            ggplot2::aes(x = log_inhibitor, y = mean_response,
                                         color = compound),
                            size = point_size, shape = point_shape, alpha = point_alpha)
    }

    if (show_error_bars && "sd_response" %in% colnames(plot_data$points)) {
      p <- p +
        ggplot2::geom_errorbar(data = plot_data$points,
                               ggplot2::aes(x = log_inhibitor,
                                            ymin = mean_response - sd_response,
                                            ymax = mean_response + sd_response,
                                            color = compound),
                               width = error_bar_width,
                               linewidth = error_linewidth, alpha = error_alpha)
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
        compound = selected_compounds[[i]]$display_name,
        stringsAsFactors = FALSE
      )
    })
    ic50_data <- do.call(rbind, Filter(Negate(is.null), ic50_rows))
    if (!is.null(ic50_data) && nrow(ic50_data) > 0) {
      p <- p + ggplot2::geom_vline(
        data = ic50_data,
        ggplot2::aes(xintercept = .data$log_ic50, color = .data$compound),
        linetype = ic50_linetype, linewidth = ic50_linewidth, alpha = ic50_line_alpha
      )
    }
  }

  # Apply colors and smart legend labels
  p <- p + ggplot2::scale_color_manual(
    values = setNames(colors_final[1:n_valid_compounds], compound_labels),
    labels = setNames(wrapped_labels, compound_labels)
  )

  # Configure legend guide
  # NOTE: row-vs-column fill order (`byrow`) is applied via theme(legend.byrow=)
  # below, not here -- ggplot2 >= 4.0 removed the `byrow` argument from
  # guide_legend() and ignores it silently.
  guide_args <- list(
    ncol = legend_ncol,
    override.aes = list(
      size = point_size,
      linetype = 1,
      linewidth = legend_key_linewidth,
      fill = NA
    ),
    title = legend_title_final
  )

  # Add shapes to legend override when shape mapping is active
  if (nrow(plot_data$points) > 0 && use_shape_mapping) {
    guide_args$override.aes$shape <- point_shapes[1:n_valid_compounds]
  }

  p <- p + ggplot2::guides(color = do.call(ggplot2::guide_legend, guide_args))

  # ============================================================================
  # 12. CONFIGURE PLOT THEME
  # ============================================================================

  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = ifelse(show_legend, legend_position, "none"),
      plot.title = ggplot2::element_text(hjust = title_hjust, face = "bold",
                                         size = if (!is.null(plot_title_size)) plot_title_size else axis_title_size + 2),
      axis.text = ggplot2::element_text(color = axis_text_color, size = axis_text_size),
      axis.title = ggplot2::element_text(color = axis_title_color, size = axis_title_size, face = "bold",
                                         vjust = axis_vjust),
      axis.line.x.bottom = ggplot2::element_blank(),
      axis.line.y.left   = ggplot2::element_blank(),
      axis.line.x.top    = ggplot2::element_blank(),
      axis.line.y.right  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = axis_line_color, linewidth = axis_line_width),
      axis.ticks.length = ggplot2::unit(tick_length, "cm"),
      legend.text = ggplot2::element_text(size = legend_text_size, lineheight = legend_lineheight),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.byrow = byrow,
      legend.key = ggplot2::element_rect(fill = legend_key_fill, color = NA),
      legend.background = ggplot2::element_rect(fill = legend_bg_fill, color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 12, r = 8, b = 8, l = 8, unit = "pt")
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
      panel.grid.major = ggplot2::element_line(color = grid_color, linewidth = grid_linewidth),
      panel.grid.minor = ggplot2::element_line(color = grid_minor_color, linewidth = grid_linewidth * 0.5)
    )
  }

  p <- p + base_theme

  # Optional theme tweaks (only applied when explicitly set, so the defaults
  # leave the appearance unchanged).
  if (!is.null(aspect_ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect_ratio)
  }
  if (!is.null(legend_spacing)) {
    p <- p + ggplot2::theme(legend.key.spacing.y = ggplot2::unit(legend_spacing, "pt"))
  }
  if (!is.null(plot_margin)) {
    p <- p + ggplot2::theme(plot.margin = plot_margin)
  }

  # Draw axis lines manually so they stop exactly at the data limits.
  # ggplot2's axis.line elements always span the full panel edge regardless of
  # coord_cartesian / expand, so we blank them above and draw our own here.
  # geom_segment with explicit data is more robust than annotate() under
  # coord_cartesian.
  axis_segs_mc <- data.frame(
    x    = c(x_limits_final[1],   x_limits_final[1]),
    xend = c(x_limits_final[2],   x_limits_final[1]),
    y    = c(y_limits_final[1],   y_limits_final[1]),
    yend = c(y_limits_final[1],   y_limits_final[2])
  )
  p <- p +
    ggplot2::geom_segment(
      data = axis_segs_mc,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      colour = axis_line_color, linewidth = axis_line_width,
      inherit.aes = FALSE)

  # ============================================================================
  # 13. DISPLAY AND SAVE
  # ============================================================================

  # --- Legend width alignment ---
  # Measure the rendered legend width and pad to a target so that
  # multiple plots with different label lengths have identical panel sizes.
  # Only right- and left-positioned legends affect panel width;
  # bottom/top legends span the full panel width and need no padding.
  legend_width_cm <- measure_legend_width(p)

  if (is.character(legend_width) && legend_width == "auto") {
    # "auto" mode: don't pad, just record the measured width in metadata
    # so the caller can determine the max across a set of plots
    if (verbose) message("Measured legend width: ", round(legend_width_cm, 3), " cm")
  } else if (!is.null(legend_width) && is.numeric(legend_width)) {
    if (legend_width_cm == 0) {
      # No side-positioned legend detected (bottom/top legend or no legend)
      if (verbose) message("legend_width has no effect: no side-positioned legend detected")
    } else if (legend_width_cm < legend_width) {
      pad_cm <- legend_width - legend_width_cm
      # Pad on the side opposite the legend to extend the legend column
      if (legend_position == "left") {
        p <- p + ggplot2::theme(legend.box.margin = ggplot2::margin(l = pad_cm, unit = "cm"))
      } else {
        p <- p + ggplot2::theme(legend.box.margin = ggplot2::margin(r = pad_cm, unit = "cm"))
      }
      if (verbose) message("Padded legend from ", round(legend_width_cm, 3), " cm to ", legend_width, " cm")
      legend_width_cm <- legend_width
    } else if (verbose) {
      message("Legend width (", round(legend_width_cm, 3), " cm) already >= target (", legend_width, " cm); no padding applied")
    }
  }


  print(p)

  if (!is.null(save_plot)) {
    # Regex matches "name.ext" where <ext> is 2-5 word characters (covers
    # png/pdf/svg/jpg/jpeg/tiff/tif/eps/bmp/webp).  Anchored to the end,
    # case-insensitive.
    ext_pattern <- "\\.[[:alnum:]]{2,5}$"

    if (is.character(save_plot)) {
      # Explicit path.  Two sub-cases: has extension, or doesn't.
      has_ext <- grepl(ext_pattern, save_plot)
      if (has_ext) {
        filename <- save_plot
        path_ext <- tolower(sub("^.*\\.", "", save_plot))
        if (!identical(path_ext, tolower(format)) && isTRUE(verbose)) {
          message("save_plot has extension '.", path_ext,
                  "'; ignoring `format = \"", format, "\"`.")
        }
      } else {
        # No extension -> use `format`
        filename <- paste0(save_plot, ".", format)
      }
    } else if (is.logical(save_plot) && save_plot) {
      indices_str <- paste(selected_indices, collapse = "_")
      filename <- paste0("multiple_compounds_",
                         if (!is.null(plate)) paste0(plate, "_") else "",
                         indices_str, ".", format)
    } else {
      stop("save_plot must be either a file path (character), TRUE for ",
           "auto-naming, or NULL to skip saving.  Got: ",
           paste(deparse(save_plot), collapse = " "),
           call. = FALSE)
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
    legend_width_cm = legend_width_cm,
    wrapped_labels = wrapped_labels,
    x_limits = x_limits,
    plot_dimensions = c(width = plot_width, height = plot_height, dpi = plot_dpi),
    file_saved = if (!is.null(save_plot)) filename else NULL,
    show_grid = show_grid,
    available_palettes = list_available_palettes()
  )

  attr(p, "metadata") <- metadata

  if (verbose && !is.null(color_palette)) {
    message("\nTip: To see all available palettes, run: names(attr(p, 'metadata')$available_palettes)")
  }

  return(p)
}
