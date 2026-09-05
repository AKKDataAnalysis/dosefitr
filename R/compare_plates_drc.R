# Internal mode normaliser used by compare_plates_drc() and its unit tests.
.normalise_compare_plates_mode <- function(compare_by) {
  if (!is.character(compare_by) || length(compare_by) != 1L ||
      is.na(compare_by) || !nzchar(trimws(compare_by))) {
    stop("compare_by must be one non-empty character value.")
  }

  mode <- tolower(trimws(compare_by))
  mode <- gsub("[- ]", "_", mode)
  aliases <- c(
    pair = "construct_compound",
    pairs = "construct_compound",
    both = "construct_compound",
    constructcompound = "construct_compound"
  )
  if (mode %in% names(aliases)) mode <- unname(aliases[[mode]])

  valid <- c("construct_compound", "compound", "construct")
  if (!mode %in% valid) {
    stop(
      "compare_by must be one of 'construct_compound', 'compound', or 'construct'."
    )
  }
  mode
}

# Internal grouping key kept separate so grouping semantics can be unit tested
# without invoking the plotting stack.
.compare_plates_group_key <- function(entry, compare_by) {
  mode <- .normalise_compare_plates_mode(compare_by)
  switch(
    mode,
    construct_compound = paste0(entry$construct, ":", entry$compound),
    compound = entry$compound,
    construct = entry$construct
  )
}

# Internal legend-mode normaliser. "auto" preserves the mode-dependent labels
# used before legend_by was introduced.
.normalise_compare_plates_legend <- function(legend_by) {
  if (!is.character(legend_by) || length(legend_by) != 1L ||
      is.na(legend_by) || !nzchar(trimws(legend_by))) {
    stop("legend_by must be one non-empty character value.")
  }
  mode <- tolower(trimws(legend_by))
  valid <- c("auto", "construct", "compound")
  if (!mode %in% valid) {
    stop("legend_by must be one of 'auto', 'construct', or 'compound'.")
  }
  mode
}

#' Compare Dose-Response Curves Across Plates
#'
#' @description
#' `compare_plates_drc()` takes the output of [`batch_drc_analysis()`] and
#' generates overlay plots according to the grouping selected in `compare_by`,
#' showing the fitted dose-response curves from every eligible plate.
#' This is the primary tool for assessing inter-plate reproducibility.
#'
#' The legacy/default mode keeps each full `Construct:Compound` pair separate.
#' Alternatively, curves may be grouped by compound across constructs or by
#' construct across compounds. Mean +- SD data points and error bars are
#' included by default.
#'
#' @param batch_drc_result The list returned by [`batch_drc_analysis()`].
#'   Must contain a `drc_results` element (one entry per plate), each with
#'   a `drc_result` sub-list that holds `detailed_results`.
#' @param compare_by Character. Controls how fits are grouped across plates.
#'   \describe{
#'     \item{`"construct_compound"`}{Default and backward-compatible mode.
#'       Compares only identical `Construct:Compound` pairs. Legend labels show
#'       the plate name. Alias: `"pair"`.}
#'     \item{`"compound"`}{Compares every successful fit having the same
#'       compound name, even when constructs differ. Legend labels show both
#'       plate and construct.}
#'     \item{`"construct"`}{Compares every successful fit having the same
#'       construct / target, even when compounds differ. Legend labels show
#'       both plate and compound.}
#'   }
#' @param output_dir Character.  Directory where the PNG files are saved.
#'   Created automatically if it does not exist.  Defaults to
#'   `"plate_comparison_plots"` inside the current working directory.
#' @param plot_width,plot_height Numeric.  Plot dimensions in inches.
#'   Defaults: `10` x `8`.
#' @param plot_dpi Numeric.  Resolution for saved PNGs.  Default: `600`.
#' @param y_limits Numeric vector of length 2 for fixed y-axis limits, or
#'   `NULL` (default) to let each plot auto-scale.
#' @param y_axis_title Character.  Y-axis label.  If \code{NULL} (default),
#'   auto-detected from the batch result: \code{"Cell Viability (\%)"} or
#'   \code{"Luminescence"} for viability assays, \code{"Normalised BRET ratio [\%]"}
#'   or \code{"BRET ratio"} for NanoBRET assays (depending on whether
#'   \code{normalize} was \code{TRUE} or \code{FALSE}).
#' @param y_number_format Character.  Y-axis tick-label number format,
#'   forwarded to [`plot_multiple_compounds()`]. One of \code{"integer"}
#'   (default; plain integers, e.g. \code{100000}), \code{"scientific"}
#'   (e.g. \code{1e+05}), or \code{"si"} (SI short-scale suffixes, e.g.
#'   \code{100K}, \code{1.5M}).
#' @param color_palette Character.  Any palette name accepted by
#'   [`plot_multiple_compounds()`] (e.g. `"set1"`, `"okabe_ito"`,
#'   `"nature"`).  Default: `"set1"`.
#' @param show_error_bars Logical.  Show mean +- SD error bars on data
#'   points.  Default: `TRUE`.
#' @param show_grid Logical.  Show background grid lines.  Default: `FALSE`.
#' @param legend_position Character.  One of `"right"`, `"left"`, `"top"`,
#'   `"bottom"`, `"none"`.  Default: `"right"`.
#' @param legend_by Character. Controls the text used for each curve in the
#'   legend. `"auto"` (default) preserves the current behaviour: plate name in
#'   `compare_by = "construct_compound"`, `"plate (construct)"` in
#'   `compare_by = "compound"`, and `"plate (compound)"` in
#'   `compare_by = "construct"`. Use `"construct"` or `"compound"` to show
#'   only that component. Curve identities remain distinct internally even
#'   when two displayed legend labels are equal.
#' @param legend_title Character.  Title printed above the legend.
#'   `"auto"` (default) shows `"Plate"`, `"Construct"`, or `"Compound"`
#'   according to `legend_by`. Supply another string to override it.
#' @param legend_text_size Numeric. Font size for the legend item labels.
#'   Default: `11`.
#' @param legend_title_size Numeric.  Font size for the legend title.
#'   Default: `11`.
#' @param legend_ncol Integer.  Number of columns in the legend.  `NULL`
#'   (default) uses adaptive logic from [`plot_multiple_compounds()`].
#' @param legend_label_wrap Integer.  Maximum character width before legend
#'   labels wrap to a new line.  Default: `25`.
#' @param show_legend Logical.  Whether to display the legend.  Default: `TRUE`.
#' @param axis_title_size Numeric.  Font size for axis titles.  Default: `14`.
#' @param axis_title_color Character.  Colour of axis title text.  Default:
#'   `"black"`.
#' @param axis_text_color Character.  Colour of axis tick labels.  Default:
#'   `"black"`.
#' @param colors Either `NULL` (use `color_palette`), a character vector of
#'   hex/colour names (one per plate), or `TRUE` for automatic ggplot2 hues.
#' @param point_shapes Controls point shapes. `NULL` (default) uses shape 16
#'   (filled circle) for all plates. `TRUE` activates a set of default
#'   distinct shapes, one per plate. A numeric vector assigns specific shape
#'   codes (recycled if needed).
#' @param point_stroke Numeric or \code{NULL}. Outline (border) width of the
#'   plotting symbols, forwarded to [`plot_multiple_compounds()`]. Thickens the
#'   outline of OPEN/hollow shapes (e.g. codes 0,1,2,5,6 and asterisk 8) so they
#'   are easier to see; filled shapes (16,17,15,18) are unaffected. \code{NULL}
#'   (default) uses a built-in width of 0.8 (bolder than ggplot2's 0.5). Try
#'   \code{1.2} for even bolder outlines when mixing open and filled shapes.
#' @param error_bar_width Numeric.  Width of the error bar caps.  Default:
#'   `0.05`.
#' @param axis_text_size Numeric.  Font size for axis tick labels.
#'   Default: `12`.
#' @param point_size Numeric.  Override the automatic point size.  `NULL`
#'   (default) uses the same adaptive logic as [`plot_multiple_compounds()`].
#' @param selected_entities Character vector. Restrict the comparison to a
#'   subset of grouping keys: `Construct:Compound` values in
#'   `compare_by = "construct_compound"`, compound names in
#'   `compare_by = "compound"`, or construct names in
#'   `compare_by = "construct"`. `NULL` processes all entities.
#' @param min_plates Integer.  Minimum number of plates that must contain an
#'   entity for a plot to be generated.  Default: `2` (skip entities that
#'   appear on only one plate, as there is nothing to compare).
#' @param x_limits Numeric vector of length 2 specifying the x-axis limits.
#'   `NULL` (default) auto-calculates from the data.  Unit is controlled by
#'   `x_limits_scale`.  Passed directly to [`plot_multiple_compounds()`].
#' @param x_limits_scale Character.  Unit of `x_limits`: `"log10"` (default),
#'   `"molar"`, `"uM"`, or `"nM"`.  Ignored when `x_limits = NULL`.
#' @param x_axis_title Character.  Custom x-axis label.  `NULL` (default) uses
#'   the standard \eqn{Log_{10}} Concentration \[M\] label.
#' @param curve_linewidth Numeric.  Line width of the fitted curves.
#'   Default: `1`.
#' @param curve_alpha Numeric (0-1).  Opacity of the fitted curves.
#'   Default: `0.7`.
#' @param show_ic50_lines Logical.  Draw a vertical dashed line at each
#'   compound's IC50.  Default: `FALSE`.
#' @param plot_title_size Numeric.  Font size for the plot title.  `NULL`
#'   (default) uses `axis_title_size + 2`.
#' @param axis_line_color Character.  Colour of axis lines, ticks, and (when
#'   `show_border = TRUE`) the panel border.  Default: `"black"`.
#' @param show_border Logical.  Draw a rectangular border around the plot
#'   panel.  Default: `FALSE`.
#' @param transparent_background Logical.  If \code{TRUE}, the plot and panel
#'   backgrounds are set to transparent.  Passed directly to
#'   [`plot_multiple_compounds()`].  Default: \code{FALSE}.
#' @param label_sep Character separator used in display labels between
#'   construct and compound names.  Defaults to \code{":"}.  Change to
#'   e.g. \code{"/"} to show \code{"EPHA1/KK135"} instead of
#'   \code{"EPHA1:KK135"} in plot titles and legends.  The internal data
#'   always uses \code{":"}; this parameter only affects display.
#' @param legend_width Numeric, \code{"auto"}, or \code{NULL}.  Target width
#'   (in cm) for the legend column in each comparison plot.  When specified,
#'   the legend is padded via \code{legend.box.margin} so that the data panel
#'   is identically sized across all comparison plots, regardless of plate-name
#'   length.  \code{"auto"} performs a two-pass approach: first measures all
#'   legend widths, then re-renders every plot padded to the maximum.  Only
#'   affects right- and left-positioned legends.  \code{NULL} (default)
#'   disables padding.
#' @param verbose Logical.  Print progress messages.  Default: `TRUE`.
#'
#' @return Invisibly returns a named list with one entry per entity plotted.
#'   Each entry contains:
#'   \describe{
#'     \item{`entity`}{The grouping key selected by `compare_by`.}
#'     \item{`constructs`, `compounds`}{Character vectors represented in the plot.}
#'     \item{`compare_by`}{The normalized grouping mode.}
#'     \item{`legend_by`}{The normalized legend-label mode.}
#'     \item{`plates`}{Character vector of plate names included.}
#'     \item{`n_plates`}{Number of plates overlaid.}
#'     \item{`plot`}{The `ggplot2` object.}
#'     \item{`file`}{Full path to the saved PNG.}
#'   }
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_cpd")
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
#' out_dir <- file.path(tempdir(), "dosefitr_ex_cpd_out")
#' dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#'
#' cmp <- compare_plates_drc(
#'   batch_drc_result  = drc_res,
#'   compare_by        = "construct_compound",
#'   output_dir        = out_dir,
#'   # Only two bridging compounds ("construct:compound" keys) so the example
#'   # stays under 5s.
#'   selected_entities = c("KinaseA:Cpd1", "KinaseB:Cpd9"),
#'   min_plates        = 2,
#'   plot_dpi          = 72,
#'   verbose           = FALSE
#' )
#' length(cmp)
#' }
#' @seealso
#' * [`batch_drc_analysis()`] - upstream function whose output is consumed here.
#' * [`plot_multiple_compounds()`] - underlying plotting engine.
#'
#' @export

compare_plates_drc <- function(batch_drc_result,
                               compare_by        = "construct_compound",
                               output_dir        = "plate_comparison_plots",
                               plot_width        = 10,
                               plot_height       = 8,
                               plot_dpi          = 600,
                               y_limits          = NULL,
                               y_axis_title      = NULL,
                               y_number_format   = c("integer", "scientific", "si"),
                               color_palette     = "set1",
                               show_error_bars   = TRUE,
                               show_grid         = FALSE,
                               legend_position   = "right",
                               legend_by         = "auto",
                               legend_title      = "auto",
                               legend_text_size  = 11,
                               legend_title_size = 11,
                               legend_ncol       = NULL,
                               legend_label_wrap = 25,
                               show_legend       = TRUE,
                               axis_title_size   = 14,
                               axis_title_color  = "black",
                               axis_text_size    = 12,
                               axis_text_color   = "black",
                               colors            = NULL,
                               point_size        = NULL,
                               point_shapes      = NULL,
                               point_stroke      = NULL,
                               error_bar_width   = 0.05,
                               selected_entities      = NULL,
                               min_plates             = 2,
                               x_limits               = NULL,
                               x_limits_scale         = "log10",
                               x_axis_title           = NULL,
                               curve_linewidth        = 1,
                               curve_alpha            = 0.7,
                               show_ic50_lines        = FALSE,
                               plot_title_size        = NULL,
                               axis_line_color        = "black",
                               show_border            = FALSE,
                               transparent_background = FALSE,
                               label_sep              = ":",
                               legend_width           = NULL,
                               verbose                = TRUE) {
  
  # ============================================================================
  # 1. VALIDATION
  # ============================================================================
  
  compare_by <- .normalise_compare_plates_mode(compare_by)
  legend_by  <- .normalise_compare_plates_legend(legend_by)

  legend_title_final <- if (identical(legend_title, "auto")) {
    switch(
      legend_by,
      auto = "Plate",
      construct = "Construct",
      compound = "Compound"
    )
  } else {
    legend_title
  }
  
  if (!is.list(batch_drc_result) || is.null(batch_drc_result$drc_results))
    stop("batch_drc_result must be the list returned by batch_drc_analysis().")
  
  drc_results <- batch_drc_result$drc_results
  if (length(drc_results) == 0)
    stop("No plate results found in batch_drc_result$drc_results.")
  
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  
  # ============================================================================
  # 2. HELPERS
  # ============================================================================
  
  # Safe fallback operator
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  
  # Parse "Construct:Compound" or plain name into target + compound parts
  parse_name <- function(raw_name) {
    # Strip trailing replicate suffix (.2, .3, ...) added by split_replicates
    clean <- gsub("\\.(\\d+)$", "", trimws(raw_name))
    if (grepl(":", clean)) {
      parts <- strsplit(clean, ":")[[1]]
      list(construct = trimws(parts[1]),
           compound  = trimws(paste(parts[-1], collapse = ":")))
    } else {
      list(construct = clean, compound = clean)
    }
  }
  
  # Sanitise a string for use as a filename
  safe_filename <- function(x, max_len = 60) {
    x <- gsub("[^A-Za-z0-9_.-]", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    if (nchar(x) > max_len) x <- substr(x, 1, max_len)
    x
  }
  
  # ============================================================================
  # 3. COLLECT ALL ENTRIES ACROSS PLATES
  # ============================================================================
  # Build a flat list: each element = one successful compound fit on one plate.
  # Fields: plate_name, construct, compound, full_label, result (the fit object)
  
  all_entries <- list()
  
  for (plate_name in names(drc_results)) {
    plate_obj  <- drc_results[[plate_name]]
    drc_result <- plate_obj$drc_result
    
    if (is.null(drc_result) || is.null(drc_result$detailed_results)) next
    if (!is.list(drc_result$detailed_results))                         next
    
    # Read normalisation flag at the plate level (most reliable source)
    plate_is_norm <- isTRUE(drc_result$used_normalized_data) ||
      isTRUE(drc_result$normalized)
    
    for (res in drc_result$detailed_results) {
      if (is.null(res$success) || !isTRUE(res$success)) next
      if (is.null(res$model))                            next
      
      # Compound label as stored: "Construct:Compound | replicate" or plain
      raw_label <- strsplit(res$compound %||% "Unknown", " \\| ")[[1]][1]
      parsed    <- parse_name(raw_label)
      
      all_entries[[length(all_entries) + 1]] <- list(
        plate_name = plate_name,
        construct  = parsed$construct,
        compound   = parsed$compound,
        full_label = raw_label,
        is_norm    = plate_is_norm,
        result     = res
      )
    }
  }
  
  if (length(all_entries) == 0)
    stop("No successful compound fits found across any plate.")
  
  # -- Filter out NA placeholder entries --------------------------------------
  # Removes entries where construct or compound is exactly NA, NA_2, NA_3, ...
  # (unnamed columns from the data). Uses a strict anchored regex so that real
  # gene/compound names containing "NA" (e.g. NAGA, CANAL, NAT1) are kept.
  na_pattern <- "^NA(_\\d+)?$"
  n_before <- length(all_entries)
  all_entries <- Filter(function(e) {
    construct_ok <- !is.na(e$construct) && !grepl(na_pattern, e$construct)
    compound_ok  <- !is.na(e$compound)  && !grepl(na_pattern, e$compound)
    construct_ok && compound_ok
  }, all_entries)
  n_removed <- n_before - length(all_entries)
  
  if (length(all_entries) == 0)
    stop("No valid entries remain after removing NA placeholders.")
  
  if (verbose) {
    message(sprintf("Collected %d successful fits across %d plates.",
                    n_before, length(drc_results)))
    if (n_removed > 0)
      message(sprintf("  Removed %d NA placeholder entr%s (construct or compound matched ^NA(_\\d+)?$).",
                      n_removed, if (n_removed == 1) "y" else "ies"))
  }
  
  # ============================================================================
  # 4. GROUP BY CHOSEN DIMENSION
  # ============================================================================
  
  # Group according to the user-selected comparison scope.
  group_key <- function(entry) {
    .compare_plates_group_key(entry, compare_by)
  }
  
  entity_map <- list()   # grouping key -> list of entries
  for (entry in all_entries) {
    key <- group_key(entry)
    entity_map[[key]] <- c(entity_map[[key]], list(entry))
  }
  
  # Apply selected_entities filter
  if (!is.null(selected_entities)) {
    entity_map <- entity_map[names(entity_map) %in% selected_entities]
    if (length(entity_map) == 0)
      stop("None of the selected_entities were found in the data.")
  }
  
  # Apply min_plates filter - count distinct plates per entity
  n_plates_per_entity <- sapply(entity_map, function(entries) {
    length(unique(sapply(entries, `[[`, "plate_name")))
  })
  
  skipped <- names(n_plates_per_entity)[n_plates_per_entity < min_plates]
  entity_map <- entity_map[n_plates_per_entity >= min_plates]
  
  if (verbose && length(skipped) > 0)
    message(sprintf(
      "Skipping %d entit%s with fewer than %d plates: %s",
      length(skipped),
      if (length(skipped) == 1) "y" else "ies",
      min_plates,
      paste(skipped, collapse = ", ")
    ))
  
  if (length(entity_map) == 0)
    stop(sprintf(
      "No entities found on >= %d plates. Lower min_plates or check your data.",
      min_plates
    ))
  
  if (verbose)
    message(sprintf("Will generate %d comparison plot(s), grouped by %s.",
                    length(entity_map), compare_by))
  
  # ============================================================================
  # 5. BUILD A SYNTHETIC fit_drc RESULT FOR plot_multiple_compounds()
  # ============================================================================

  build_synthetic_result <- function(entries) {
    # Use the plate-level normalisation flag collected during data gathering
    is_norm <- any(sapply(entries, function(e) isTRUE(e$is_norm)))
    
    synthetic_detailed <- lapply(seq_along(entries), function(entry_idx) {
      e <- entries[[entry_idx]]
      r <- e$result
      legend_label <- if (legend_by == "auto") {
        switch(
          compare_by,
          construct_compound = e$plate_name,
          compound = paste0(e$plate_name, " (", e$construct, ")"),
          construct = paste0(e$plate_name, " (", e$compound, ")")
        )
      } else {
        switch(
          legend_by,
          construct = e$construct,
          compound = e$compound
        )
      }

      # The left-hand component is a unique internal curve identifier. The
      # right-hand component is the user-facing legend label selected above.
      # plot_multiple_compounds() is asked to display only that right-hand part.
      r$compound <- paste0(sprintf("curve_%03d", entry_idx), ":", legend_label)
      r
    })
    
    list(
      detailed_results = synthetic_detailed,
      normalized       = is_norm,
      legend_label     = "compound"
    )
  }
  
  # ============================================================================
  # ============================================================================
  # 6. PLOT LOOP
  # ============================================================================

  # -- Pre-compute per-entity metadata (title, filename, y-axis title) ----------
  entity_meta <- list()
  for (ei in seq_along(entity_map)) {
    entity_name <- names(entity_map)[ei]
    entries     <- entity_map[[entity_name]]

    constructs <- unique(vapply(entries, `[[`, character(1L), "construct"))
    compounds  <- unique(vapply(entries, `[[`, character(1L), "compound"))

    entity_construct <- if (length(constructs) == 1L) constructs[[1L]] else NA_character_
    entity_compound  <- if (length(compounds) == 1L) compounds[[1L]] else NA_character_

    plates_present <- unique(sapply(entries, `[[`, "plate_name"))
    n_pl           <- length(plates_present)

    synthetic <- build_synthetic_result(entries)

    y_title_final <- if (!is.null(y_axis_title)) {
      y_axis_title
    } else {
      assay_src  <- batch_drc_result$metadata$assay_type
      normalized <- if (!is.null(batch_drc_result$metadata$normalize))
        batch_drc_result$metadata$normalize
      else isTRUE(synthetic$normalized)
      if (!is.null(assay_src) && assay_src == "viability") {
        if (normalized) "Cell Viability (%)" else "Luminescence"
      } else {
        if (normalized) "Normalised BRET ratio [%]" else "BRET ratio"
      }
    }

    plot_title_str <- switch(
      compare_by,
      construct_compound = paste(constructs[[1L]], compounds[[1L]], sep = " - "),
      compound = compounds[[1L]],
      construct = constructs[[1L]]
    )

    file_stem <- switch(
      compare_by,
      construct_compound = paste0(
        safe_filename(constructs[[1L]]), "__", safe_filename(compounds[[1L]])
      ),
      compound = paste0("compound__", safe_filename(compounds[[1L]])),
      construct = paste0("construct__", safe_filename(constructs[[1L]]))
    )
    fname <- file.path(output_dir, paste0(file_stem, ".png"))

    display_label <- switch(
      compare_by,
      construct_compound = paste(constructs[[1L]], compounds[[1L]], sep = " - "),
      compound = compounds[[1L]],
      construct = constructs[[1L]]
    )

    entity_meta[[entity_name]] <- list(
      construct     = entity_construct,
      compound      = entity_compound,
      constructs    = constructs,
      compounds     = compounds,
      compare_by    = compare_by,
      display_label = display_label,
      plates        = plates_present,
      n_plates      = n_pl,
      synthetic     = synthetic,
      y_title       = y_title_final,
      plot_title    = plot_title_str,
      fname         = fname
    )
  }

  total <- length(entity_map)

  # -- Helper: call plot_multiple_compounds for one entity ----------------------
  plot_one_entity <- function(entity_name, lw, save_file) {
    meta      <- entity_meta[[entity_name]]
    synthetic <- meta$synthetic
    p         <- NULL

    tryCatch({
      capture.output(
        { p <- suppressMessages(
            plot_multiple_compounds(
              results         = synthetic,
              plot_title      = meta$plot_title,
              y_axis_title    = meta$y_title,
              y_limits        = y_limits,
              y_number_format = y_number_format,
              color_palette   = color_palette,
              show_error_bars = show_error_bars,
              show_grid       = show_grid,
              legend_position = legend_position,
              legend_title      = legend_title_final,
              legend_label      = synthetic$legend_label,
              legend_text_size  = legend_text_size,
              legend_title_size = legend_title_size,
              legend_ncol       = legend_ncol,
              legend_label_wrap = legend_label_wrap,
              show_legend       = show_legend,
              axis_title_size   = axis_title_size,
              axis_title_color  = axis_title_color,
              axis_text_size    = axis_text_size,
              axis_text_color   = axis_text_color,
              colors            = colors,
              point_size        = point_size,
              point_shapes      = point_shapes,
              point_stroke      = point_stroke,
              error_bar_width   = error_bar_width,
              x_limits               = x_limits,
              x_limits_scale         = x_limits_scale,
              x_axis_title           = x_axis_title,
              curve_linewidth        = curve_linewidth,
              curve_alpha            = curve_alpha,
              show_ic50_lines        = show_ic50_lines,
              plot_title_size        = plot_title_size,
              axis_line_color        = axis_line_color,
              show_border            = show_border,
              transparent_background = transparent_background,
              save_plot              = save_file,
              plot_width             = plot_width,
              plot_height            = plot_height,
              plot_dpi               = plot_dpi,
              label_sep              = label_sep,
              legend_width           = lw,
              verbose                = FALSE
            )
          )
        },
        type = "output"
      )
    }, error = function(e) {
      warning(sprintf("Failed to plot '%s': %s", entity_name, e$message))
    })
    p
  }

  # ============================================================================
  # 6a. PASS 1: measure legend widths (only when legend_width = "auto")
  # ============================================================================
  max_legend_cm <- 0

  if (is.character(legend_width) && legend_width == "auto") {
    if (verbose) message("\nMeasuring legend widths (pass 1 of 2)...")
    for (ei in seq_along(entity_map)) {
      entity_name <- names(entity_map)[ei]
      meta        <- entity_meta[[entity_name]]

      if (verbose)
        message(sprintf("  [%d/%d] Measuring '%s'",
                        ei, total, meta$display_label))

      p <- plot_one_entity(entity_name, lw = "auto", save_file = NULL)
      if (!is.null(p)) {
        w <- attr(p, "metadata")$legend_width_cm
        if (!is.null(w) && w > max_legend_cm) max_legend_cm <- w
      }
    }
    if (verbose)
      message(sprintf("  Maximum legend width: %.3f cm", max_legend_cm))
  }

  # ============================================================================
  # 6b. PASS 2: render and save all plots
  # ============================================================================
  output_list <- list()

  # Determine the effective legend_width for pass 2
  effective_lw <- if (is.character(legend_width) && legend_width == "auto") {
    if (max_legend_cm > 0) max_legend_cm else NULL
  } else {
    legend_width   # NULL or numeric — pass through unchanged
  }

  if (is.character(legend_width) && legend_width == "auto" && verbose)
    message("\nRendering plots with legend_width = ", 
            if (is.null(effective_lw)) "NULL" else sprintf("%.3f cm", effective_lw),
            " (pass 2 of 2)...")

  for (ei in seq_along(entity_map)) {
    entity_name <- names(entity_map)[ei]
    meta        <- entity_meta[[entity_name]]

    if (verbose)
      message(sprintf("\n[%d/%d] '%s'  (%d plate%s: %s)",
                      ei, total, meta$display_label,
                      meta$n_plates, if (meta$n_plates == 1) "" else "s",
                      paste(meta$plates, collapse = ", ")))

    p <- plot_one_entity(entity_name, lw = effective_lw, save_file = meta$fname)

    if (!is.null(p)) {
      if (verbose) message("  Saved: ", meta$fname)
      output_list[[entity_name]] <- list(
        entity    = entity_name,
        construct = meta$construct,
        compound  = meta$compound,
        constructs = meta$constructs,
        compounds  = meta$compounds,
        compare_by = meta$compare_by,
        legend_by  = legend_by,
        plates    = meta$plates,
        n_plates  = meta$n_plates,
        plot      = p,
        file      = meta$fname
      )
    }
  }

  # ============================================================================
  # 7. SUMMARY
  # ============================================================================

  n_ok   <- length(output_list)
  n_fail <- total - n_ok

  if (verbose) {
    message("\n", paste(rep("=", 55), collapse = ""))
    message(sprintf("PLATE COMPARISON COMPLETE - %d plot%s saved to: %s",
                    n_ok, if (n_ok == 1) "" else "s",
                    normalizePath(output_dir)))
    if (n_fail > 0)
      message(sprintf("  %d plot%s failed (see warnings above).",
                      n_fail, if (n_fail == 1) "" else "s"))
    message(paste(rep("=", 55), collapse = ""))
  }

  invisible(output_list)

}
