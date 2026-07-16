#' Plot NanoBRET Dose-Response Curves for All Batch Plates
#'
#' Generates a multi-panel dose-response figure for every plate in a
#' batch result, with one panel per compound.  Outlier points (replaced
#' by \code{NA} in \code{\link{rout_outliers_batch}}) are
#' overlaid in a distinct colour so they remain visible but are clearly
#' flagged.  Each plate is saved as a separate PNG file.
#'
#' @param batch_rout_output Named list returned by
#'   \code{\link{rout_outliers_batch}}.  Each plate element must
#'   contain \code{$result$modified_ratio_table} (cleaned data) and,
#'   when outliers were detected, \code{$result$modified_ratio_table_original}
#'   (pre-cleaning data used to overlay the removed points).
#'
#' @param output_dir Character string.  Directory where PNG files are
#'   saved.  Defaults to the current working directory.  Created
#'   automatically if it does not exist.
#'
#' @param plates Character vector of plate names to plot (must match
#'   names in \code{batch_rout_output}).  \code{NULL} (default) plots all
#'   plates.
#'
#' @param ncol Integer.  Number of compound panels per row (default
#'   \code{4L}).  Passed to \code{patchwork::wrap_plots()}.
#'
#' @param width_per_col Numeric.  Width in inches allocated to each
#'   column of panels (default \code{3.2}).  Total figure width is
#'   \code{ncol * width_per_col}.
#'
#' @param height_per_row Numeric.  Height in inches allocated to each
#'   row of panels (default \code{3.0}).  Total figure height is
#'   \code{ceiling(n_compounds / ncol) * height_per_row}.
#'
#' @param dpi Integer.  Resolution of saved PNG files (default
#'   \code{150}).
#'
#' @param verbose Logical.  Print per-plate progress messages (default
#'   \code{TRUE}).
#'   
#' @param panel_spacing Numeric. Spacing between sub-plots in the panel, in
#'   centimetres (default \code{0.5}). Increase for more breathing room between
#'   plots.
#' @param subplot_title Character. Controls what text is used as the title of
#'   each compound sub-plot. One of \code{"auto"} (default), \code{"full"}
#'   (e.g. \code{"KinaseA:Cpd1"}), \code{"compound"} (e.g. \code{"Cpd1"}),
#'   or \code{"construct"} (e.g. \code{"KinaseA"}). \code{"auto"} shows only
#'   the compound name when all compounds share one construct, only the construct
#'   name when all share one compound, and the full string otherwise.
#' @param label_sep Character separator used in display labels between
#'   construct and compound names.  Defaults to \code{":"}.  Change to
#'   e.g. \code{"/"} to show \code{"EPHA1/KK135"} instead of
#'   \code{"EPHA1:KK135"} in plot titles and legends.  The internal data
#'   always uses \code{":"}; this parameter only affects display.
#'
#' @return Invisibly returns a character vector of the PNG file paths
#'   written to \code{output_dir} (one path per successfully processed
#'   plate).  The primary side-effect is writing those PNG files.
#'
#' @section Output files:
#' One PNG per plate, named \code{<plate_name>_curves.png}, is written to
#' \code{output_dir}.  File dimensions are computed automatically from
#' \code{ncol}, \code{width_per_col}, \code{height_per_row}, and the
#' number of compounds on the plate.
#'
#' @section Outlier overlay:
#' When \code{$result$modified_ratio_table_original} is present (set by
#' \code{\link{rout_outliers_batch}} whenever at least one
#' outlier was removed), the original values are plotted as open red
#' triangles on top of the fitted curve.  This lets you visually confirm
#' that the removed points were genuine outliers rather than biologically
#' meaningful observations.
#'
#' @section Replicate handling:
#' Columns ending in \code{".2"} are treated as the second technical
#' replicate of the corresponding \code{".1"} column.  Both replicates
#' are passed to \code{\link{plot_outliers_curves}}, which plots them as
#' separate point series with a shared fitted curve.
#'
#' @section Dependencies:
#' Requires \pkg{ggplot2}, \pkg{ggprism}, \pkg{ggrepel}, and
#' \pkg{patchwork}.  Also requires \code{\link{rout_outliers}}
#' and \code{\link{plot_outliers_curves}} to be available in the current
#' environment (source \code{rout_outliers.R} before calling
#' this function).
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' # Copy only one plate + info to keep the example under 5s.
#' extdata_dir <- system.file("extdata", package = "dosefitr")
#' work_dir    <- file.path(tempdir(), "dosefitr_ex_pobc")
#' dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
#' invisible(file.copy(
#'   file.path(extdata_dir,
#'             c("nanobret_info.xlsx", "nanobret_plate_01.xlsx")),
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
#' # Suppress benign nls "false convergence" and NaN warnings.
#' rout_res <- suppressWarnings(rout_outliers_batch(
#'   batch_results = ratio_res, Q = 0.01, n_param = 4L,
#'   direction = "inhibition", verbose = FALSE
#' ))
#'
#' # Trim the batch object to two compounds so plotting stays under 5s.
#' keep_cmpds <- c("KinaseA:Cpd1", "KinaseA:Cpd2")
#' rr         <- rout_res$plate_01$result$rout_results$results
#' rout_res$plate_01$result$rout_results$results <- rr[rr$compound %in% keep_cmpds, ]
#'
#' mrt        <- rout_res$plate_01$result$modified_ratio_table
#' keep_cols  <- c(1L, grep("KinaseA:Cpd1|KinaseA:Cpd2", colnames(mrt)))
#' rout_res$plate_01$result$modified_ratio_table <- mrt[, keep_cols, drop = FALSE]
#' if (!is.null(rout_res$plate_01$result$modified_ratio_table_original)) {
#'   mrto        <- rout_res$plate_01$result$modified_ratio_table_original
#'   keep_cols_o <- c(1L, grep("KinaseA:Cpd1|KinaseA:Cpd2", colnames(mrto)))
#'   rout_res$plate_01$result$modified_ratio_table_original <-
#'     mrto[, keep_cols_o, drop = FALSE]
#' }
#'
#' out_dir <- file.path(tempdir(), "dosefitr_ex_pobc_out")
#' dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#'
#' res <- plot_outliers_batch_curves(
#'   batch_rout_output = rout_res,
#'   output_dir        = out_dir,
#'   dpi               = 72,
#'   verbose           = FALSE
#' )
#' names(res)
#' }
#' @seealso
#' \code{\link{plot_outliers_curves}} for the single-plate plotting
#' engine.
#'
#' \code{\link{rout_outliers_batch}} for the upstream outlier
#' detection step.
#'
#' \code{\link{batch_ratio_analysis}} for the upstream plate-processing
#' step.
#'
#' @export

plot_outliers_batch_curves <- function(batch_rout_output,
                                       output_dir      = NULL,
                                       plates          = NULL,
                                       ncol            = 4L,
                                       width_per_col   = 3.2,
                                       height_per_row  = 3.0,
                                       dpi             = 150,
                                       verbose         = TRUE,
                                       panel_spacing   = 0.5,
                                       subplot_title   = "auto",
                                       label_sep       = ":") {
  
  # --------------------------------------------------------------------------
  # 1. Dependency checks
  # --------------------------------------------------------------------------
  pkgs <- c("ggplot2", "ggprism", "ggrepel", "patchwork")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing_pkgs) > 0L)
    stop(sprintf(
      "plot_outliers_batch_curves() requires: %s\n  Install with: install.packages(c(%s))",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf('"%s"', missing_pkgs), collapse = ", ")),
      call. = FALSE)
  
  if (!exists("rout_outliers", mode = "function"))
    stop(paste0("rout_outliers() not found. ",
                "Please source('rout_outliers.R') before calling this function."))
  
  if (!exists("plot_outliers_curves", mode = "function"))
    stop(paste0("plot_outliers_curves() not found. ",
                "Please source('rout_outliers.R') before calling this function."))
  
  subplot_title <- match.arg(subplot_title, c("auto", "full", "compound", "construct"))
  
  # --------------------------------------------------------------------------
  # 2. Input validation
  # --------------------------------------------------------------------------
  if (!is.list(batch_rout_output) || length(batch_rout_output) == 0L)
    stop("batch_rout_output must be the non-empty return value of rout_outliers_batch().")
  
  if (is.null(batch_rout_output$params))
    stop(paste0("batch_rout_output$params not found. ",
                "Pass the direct return value of rout_outliers_batch()."))
  
  # --------------------------------------------------------------------------
  # 3. Setup
  # --------------------------------------------------------------------------
  if (is.null(output_dir)) output_dir <- getwd()
  
  plot_dir <- file.path(output_dir, "ROUT_Plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    if (verbose) message(sprintf("Created output folder: %s", plot_dir))
  }
  
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a
  
  # Extract ROUT parameters used during batch processing
  params     <- batch_rout_output$params
  Q          <- params$Q          %||% 0.01
  n_param    <- params$n_param    %||% 4L
  direction  <- params$direction  %||% "inhibition"
  ntry_retry <- params$ntry_retry %||% 3L
  log_base   <- params$log_base   %||% "log10"
  
  # Identify plate names (exclude reserved summary elements)
  reserved    <- c("outlier_summary", "skipped_summary", "rescued_summary", "params")
  plate_names <- setdiff(names(batch_rout_output), reserved)
  
  if (length(plate_names) == 0L)
    stop("No plate entries found in batch_rout_output.")
  
  # Filter to requested subset if `plates` is specified
  if (!is.null(plates)) {
    unknown <- setdiff(plates, plate_names)
    if (length(unknown) > 0L)
      warning(sprintf("plates not found in batch_rout_output and will be ignored: %s",
                      paste(unknown, collapse = ", ")))
    plate_names <- intersect(plate_names, plates)
    if (length(plate_names) == 0L)
      stop("None of the requested plates were found in batch_rout_output.")
  }
  
  if (verbose) {
    cat(strrep("=", 60), "\n")
    cat("NANOBRET BATCH CURVE PLOTS\n")
    cat(strrep("=", 60), "\n")
    cat(sprintf("Plates to plot   : %d\n", length(plate_names)))
    cat(sprintf("Output folder    : %s\n\n", plot_dir))
  }
  
  # --------------------------------------------------------------------------
  # 4. Helper: remove NA/empty columns (mirrors batch function logic)
  # --------------------------------------------------------------------------
  .drop_na_cols <- function(tbl) {
    value_cols <- seq(2L, ncol(tbl))
    keep <- vapply(value_cols, function(j) {
      nm        <- colnames(tbl)[j]
      cmpd_part <- strsplit(nm, ":")[[1L]][1L]
      if (grepl("^NA", cmpd_part)) return(FALSE)
      vals <- tbl[[j]]
      na_frac <- if (is.numeric(vals)) mean(is.na(vals) | is.nan(vals))
      else                  mean(is.na(vals))
      na_frac <= 0.8
    }, logical(1L))
    tbl[, c(1L, value_cols[keep]), drop = FALSE]
  }
  
  # --------------------------------------------------------------------------
  # 4b. Auto-detect subplot_title mode from batch composition
  # --------------------------------------------------------------------------
  # Collect all unique compound/construct names across all plates to be plotted.
  .extract_construct <- function(s) {
    if (is.null(s) || is.na(s)) return(NA_character_)
    p <- strsplit(s, ":", fixed = TRUE)[[1L]]
    if (length(p) >= 2L) trimws(p[[1L]]) else NA_character_
  }
  .extract_compound <- function(s) {
    if (is.null(s) || is.na(s)) return(NA_character_)
    p <- strsplit(s, ":", fixed = TRUE)[[1L]]
    if (length(p) >= 2L) trimws(p[[2L]]) else trimws(s)
  }
  .is_na_name <- function(x) {
    if (is.null(x) || length(x) == 0L || is.na(x)) return(TRUE)
    toupper(sub("_\\d+$", "", trimws(x))) == "NA"
  }
  
  all_cmpd_names <- character(0)
  all_cons_names <- character(0)
  for (.pn in plate_names) {
    .res <- tryCatch(batch_rout_output[[.pn]]$result$rout_results$results,
                     error = function(e) NULL)
    if (is.null(.res)) next
    for (.cmpd in unique(.res$compound)) {
      if (.is_na_name(.extract_compound(.cmpd)) ||
          .is_na_name(.extract_construct(.cmpd))) next
      all_cmpd_names <- c(all_cmpd_names, .extract_compound(.cmpd))
      all_cons_names <- c(all_cons_names, .extract_construct(.cmpd))
    }
  }
  auto_mode <- if (length(unique(all_cons_names[!is.na(all_cons_names)])) <= 1L) {
    "compound"
  } else if (length(unique(all_cmpd_names[!is.na(all_cmpd_names)])) <= 1L) {
    "construct"
  } else {
    "full"
  }
  effective_subplot_mode <- if (subplot_title == "auto") auto_mode else subplot_title
  
  # --------------------------------------------------------------------------
  # 5. Per-plate plotting loop
  # --------------------------------------------------------------------------
  saved_files <- list()
  
  for (plate_name in plate_names) {
    
    if (verbose) cat(sprintf("Plotting %s ... ", plate_name))
    
    plate <- batch_rout_output[[plate_name]]
    
    # ---- 5a. Extract the pre-cleaning modified_ratio_table ----

    mrt_original <- tryCatch(plate$result$modified_ratio_table_original,
                             error = function(e) NULL)
    mrt_cleaned  <- tryCatch(plate$result$modified_ratio_table,
                             error = function(e) NULL)
    
    mrt <- if (!is.null(mrt_original) && is.data.frame(mrt_original) &&
               nrow(mrt_original) > 0L && ncol(mrt_original) >= 3L) {
      mrt_original
    } else {
      mrt_cleaned
    }
    
    if (is.null(mrt) || !is.data.frame(mrt) || nrow(mrt) == 0L || ncol(mrt) < 3L) {
      if (verbose) cat("SKIPPED (no valid ratio table)\n")
      next
    }
    
    # ---- 5b. Drop NA/empty columns ----
    mrt_clean <- .drop_na_cols(mrt)
    
    if (ncol(mrt_clean) < 3L) {
      if (verbose) cat("SKIPPED (no valid compound columns after NA removal)\n")
      next
    }
    
    # ---- 5c. Extract dose rows only (drop ALL NA-concentration rows) ----

    conc_vals   <- mrt_clean[[1L]]
    dose_rows   <- which(!is.na(conc_vals))
    
    if (length(dose_rows) < 2L) {
      if (verbose) cat("SKIPPED (fewer than 2 dose rows)\n")
      next
    }
    
    tbl_for_fit           <- mrt_clean[dose_rows, , drop = FALSE]
    rownames(tbl_for_fit) <- NULL
    
    # ---- 5d. Obtain ROUT results for this plate ----

    rout_out <- tryCatch(plate$result$rout_results, error = function(e) NULL)
    
    if (is.null(rout_out) || is.null(rout_out$results) ||
        nrow(rout_out$results) == 0L) {
      
      if (verbose) message(sprintf(
        "  [%s] rout_results not cached; re-running ROUT (fallback).", plate_name))
      
      rout_out <- tryCatch(
        rout_outliers(
          data              = tbl_for_fit,
          Q                 = Q,
          n_param           = n_param,
          conc_col          = 1L,
          log_base          = log_base,
          direction         = direction,
          ntry_retry        = ntry_retry,
          verbose           = FALSE
        ),
        error = function(e) {
          warning(sprintf("Re-fit failed for plate '%s': %s", plate_name, e$message))
          NULL
        }
      )
      
      if (is.null(rout_out) || nrow(rout_out$results) == 0L) {
        if (verbose) cat("SKIPPED (fit failed)\n")
        next
      }
      
      # Fallback path: still apply rescue flag-clearing so the plot is correct
      rescued_plate <- tryCatch(plate$result$rescued_cytotoxic, error = function(e) NULL)
      if (!is.null(rescued_plate) && nrow(rescued_plate) > 0L) {
        conc_col_name_fb <- intersect(c("log10_conc", "ln_conc"),
                                      names(rout_out$results))[1L]
        for (i in seq_len(nrow(rescued_plate))) {
          .mask <- rout_out$results$compound == rescued_plate$compound[i] &
            abs(rout_out$results[[conc_col_name_fb]] -
                  rescued_plate[[conc_col_name_fb]][i]) < 1e-9
          rout_out$results$outlier_fdr[.mask] <- FALSE
          rout_out$results$outlier_raw[.mask] <- FALSE
        }
      }
    }
    
    # ---- 5g. Compute plot dimensions ----
    n_compounds <- length(unique(rout_out$results$compound))
    n_rows_grid <- ceiling(n_compounds / ncol)
    plot_width  <- ncol * width_per_col
    plot_height <- n_rows_grid * height_per_row + 0.8   # +0.8 for title/caption
    
    # ---- 5f. Build plate title (include data_file if available) ----
    data_file  <- plate$data_file %||% ""
    plate_title <- if (nchar(data_file) > 0L) {
      sprintf("%s  |  %s", plate_name, data_file)
    } else {
      plate_name
    }
    
    # ---- 5g. Save PNG ----
    out_file <- file.path(plot_dir, sprintf("%s_curves.png", plate_name))
    
    tryCatch({
      plot_outliers_curves(
        rout_output    = rout_out,
        title          = plate_title,
        ncol           = ncol,
        file           = out_file,
        width          = plot_width,
        height         = plot_height,
        subplot_title  = effective_subplot_mode,
        panel_spacing  = panel_spacing,
        label_sep      = label_sep
      )
      saved_files[[plate_name]] <- out_file
      if (verbose) cat(sprintf("saved (%d compounds, %dx%d in)\n",
                               n_compounds,
                               round(plot_width), round(plot_height)))
    }, error = function(e) {
      warning(sprintf("Failed to save plot for plate '%s': %s", plate_name, e$message))
      if (verbose) cat(sprintf("FAILED (%s)\n", e$message))
    })
  }
  
  # --------------------------------------------------------------------------
  # 6. Summary
  # --------------------------------------------------------------------------
  if (verbose) {
    cat(strrep("=", 60), "\n")
    cat(sprintf("Plots saved: %d / %d plates\n", length(saved_files), length(plate_names)))
    if (length(saved_files) > 0L) {
      cat(sprintf("Location   : %s\n", plot_dir))
    }
    cat(strrep("=", 60), "\n")
  }
  
  invisible(saved_files)
}
