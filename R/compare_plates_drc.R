#' Compare Dose-Response Curves Across Plates
#'
#' @description
#' `compare_plates_drc()` takes the output of [`batch_drc_analysis()`] and
#' generates one overlay plot per unique **compound** or **construct**,
#' showing the fitted dose-response curves from every plate side-by-side.
#' This is the primary tool for assessing inter-plate reproducibility.
#'
#' Each plot overlays the curves from all plates that contain the selected
#' construct-compound pair, with each plate rendered in a distinct colour.
#' Plots are always separated by the full `Construct:Compound` combination, so
#' the same compound tested on different constructs always produces separate
#' images.  The plot title is always `"Construct - Compound"` and the filename
#' is `Construct__Compound.png`.  Mean +- SD data points and error bars are
#' included by default.
#'
#' @param batch_drc_result The list returned by [`batch_drc_analysis()`].
#'   Must contain a `drc_results` element (one entry per plate), each with
#'   a `drc_result` sub-list that holds `detailed_results`.
#' @param compare_by Character.  Controls the **plot title** shown inside
#'   each image.  Plots are always separated by the full `Construct:Compound`
#'   pair regardless of this setting, so the same compound on different
#'   constructs always produces separate files.
#'   \describe{
#'     \item{`"compound"`}{(default) Title shows the compound name.}
#'     \item{`"construct"`}{Title shows the construct / target name.}
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
#' @param color_palette Character.  Any palette name accepted by
#'   [`plot_multiple_compounds()`] (e.g. `"set1"`, `"okabe_ito"`,
#'   `"nature"`).  Default: `"set1"`.
#' @param show_error_bars Logical.  Show mean +- SD error bars on data
#'   points.  Default: `TRUE`.
#' @param show_grid Logical.  Show background grid lines.  Default: `FALSE`.
#' @param legend_position Character.  One of `"right"`, `"left"`, `"top"`,
#'   `"bottom"`, `"none"`.  Default: `"right"`.
#' @param legend_title Character.  Title printed above the legend.
#'   Default: `"Plate"`.
#' @param legend_text_size Numeric.  Font size for the legend item labels
#'   (plate names).  Default: `11`.
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
#' @param error_bar_width Numeric.  Width of the error bar caps.  Default:
#'   `0.05`.
#' @param axis_text_size Numeric.  Font size for axis tick labels.
#'   Default: `12`.
#' @param point_size Numeric.  Override the automatic point size.  `NULL`
#'   (default) uses the same adaptive logic as [`plot_multiple_compounds()`].
#' @param selected_entities Character vector.  Restrict the comparison to a
#'   subset of compound or construct names.  `NULL` (default) processes all
#'   entities found across plates.
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
#' @param verbose Logical.  Print progress messages.  Default: `TRUE`.
#'
#' @return Invisibly returns a named list with one entry per entity plotted.
#'   Each entry contains:
#'   \describe{
#'     \item{`entity`}{The compound or construct name.}
#'     \item{`plates`}{Character vector of plate names included.}
#'     \item{`n_plates`}{Number of plates overlaid.}
#'     \item{`plot`}{The `ggplot2` object.}
#'     \item{`file`}{Full path to the saved PNG.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Run batch ratio + DRC analysis first
#' ratio_res <- batch_ratio_analysis("data/")
#' drc_res   <- batch_drc_analysis(ratio_res)
#'
#' # Title shows compound name (default)
#' compare_plates_drc(drc_res, compare_by = "compound",
#'                    output_dir = "comparison_plots")
#'
#' # Title shows construct name instead
#' compare_plates_drc(drc_res, compare_by = "construct",
#'                    output_dir = "comparison_plots")
#'
#' # Only compare specific entities
#' compare_plates_drc(drc_res, compare_by = "compound",
#'                    selected_entities = c("Staurosporine", "Compound_A"))
#' }
#'
#' @seealso
#' * [`batch_drc_analysis()`] - upstream function whose output is consumed here.
#' * [`plot_multiple_compounds()`] - underlying plotting engine.
#'
#' @export

compare_plates_drc <- function(batch_drc_result,
                               compare_by        = "compound",
                               output_dir        = "plate_comparison_plots",
                               plot_width        = 10,
                               plot_height       = 8,
                               plot_dpi          = 600,
                               y_limits          = NULL,
                               y_axis_title      = NULL,
                               color_palette     = "set1",
                               show_error_bars   = TRUE,
                               show_grid         = FALSE,
                               legend_position   = "right",
                               legend_title      = "Plate",
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
                               verbose                = TRUE) {
  
  # ============================================================================
  # 1. VALIDATION
  # ============================================================================
  
  compare_by <- tolower(compare_by)
  if (!compare_by %in% c("compound", "construct"))
    stop("compare_by must be either 'compound' or 'construct'.")
  
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
  
  # Always group by the full construct:compound pair so that the same compound
  # tested on different constructs always produces separate plots.
  # The plot title is always 'Construct - Compound'.
  group_key <- function(entry) {
    paste0(entry$construct, ":", entry$compound)
  }
  
  entity_map <- list()   # "construct:compound" -> list of entries
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
    message(sprintf("Will generate %d comparison plot(s) (titles showing %s name).",
                    length(entity_map), compare_by))
  
  # ============================================================================
  # 5. BUILD A SYNTHETIC fit_drc RESULT FOR plot_multiple_compounds()
  # ============================================================================
  # plot_multiple_compounds() expects:
  #   results$detailed_results  - list of individual fit objects
  #   results$normalized        - logical
  #
  # We construct one such object per entity, where each element of
  # detailed_results is one plate's fit for that entity.  We rename the
  # compound label to the plate name so the legend reads "Plate X".
  
  build_synthetic_result <- function(entries) {
    # Use the plate-level normalisation flag collected during data gathering
    is_norm <- any(sapply(entries, function(e) isTRUE(e$is_norm)))
    
    synthetic_detailed <- lapply(entries, function(e) {
      r <- e$result
      # Relabel compound so the legend shows only the plate name.
      # plot_multiple_compounds() splits on ":" to get target vs compound;
      # using "plate:<plate_name>" ensures the legend reads just the plate name
      # (the smart-label logic drops the shared "plate" prefix automatically).
      r$compound <- paste0("plate:", e$plate_name)
      r
    })
    
    list(
      detailed_results = synthetic_detailed,
      normalized       = is_norm
    )
  }
  
  # ============================================================================
  # 6. PLOT LOOP
  # ============================================================================
  
  output_list <- list()
  total <- length(entity_map)
  
  for (ei in seq_along(entity_map)) {
    entity_name <- names(entity_map)[ei]
    entries     <- entity_map[[entity_name]]
    
    # Parse the group key back into construct + compound parts
    key_parts        <- strsplit(entity_name, ":", fixed = TRUE)[[1]]
    entity_construct <- key_parts[1]
    entity_compound  <- if (length(key_parts) > 1)
      paste(key_parts[-1], collapse = ":") else key_parts[1]
    
    plates_present <- unique(sapply(entries, `[[`, "plate_name"))
    n_pl           <- length(plates_present)
    
    if (verbose)
      message(sprintf("\n[%d/%d] '%s - %s'  (%d plate%s: %s)",
                      ei, total, entity_construct, entity_compound,
                      n_pl, if (n_pl == 1) "" else "s",
                      paste(plates_present, collapse = ", ")))
    
    # Build the synthetic result object
    synthetic <- build_synthetic_result(entries)
    
    # Auto y-axis title: use assay_type + normalize from batch metadata when
    # available; fall back to normalization-only heuristic otherwise.
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
    
    # Plot title: compound name when compare_by = "compound", construct otherwise
    plot_title_str <- if (compare_by == "compound") entity_compound else entity_construct
    
    # Filename: "Construct__Compound.png" (always encodes both to avoid collisions)
    fname <- file.path(
      output_dir,
      paste0(safe_filename(entity_construct), "__",
             safe_filename(entity_compound), ".png")
    )
    
    # Call plot_multiple_compounds - suppress its print() so we control output
    p <- tryCatch({
      suppressMessages(
        plot_multiple_compounds(
          results         = synthetic,
          plot_title      = plot_title_str,
          y_axis_title    = y_title_final,
          y_limits        = y_limits,
          color_palette   = color_palette,
          show_error_bars = show_error_bars,
          show_grid       = show_grid,
          legend_position = legend_position,
          legend_title      = legend_title,
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
          save_plot              = fname,
          plot_width             = plot_width,
          plot_height            = plot_height,
          plot_dpi               = plot_dpi,
          verbose                = FALSE
        )
      )
    }, error = function(e) {
      warning(sprintf("Failed to plot '%s': %s", entity_name, e$message))
      NULL
    })
    
    if (!is.null(p)) {
      if (verbose) message("  Saved: ", fname)
      output_list[[entity_name]] <- list(
        entity    = entity_name,
        construct = entity_construct,
        compound  = entity_compound,
        plates    = plates_present,
        n_plates  = n_pl,
        plot      = p,
        file      = fname
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
