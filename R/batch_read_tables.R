# ============================================================================
# Internal helper (not exported, not documented in man/)
# ----------------------------------------------------------------------------
# TRUE for each path that is absolute on the current platform.  Handles
# POSIX ("/..."), tilde ("~/..."), Windows drives ("C:\..." / "C:/..."),
# and Windows UNC paths ("\\server\...").  Vectorised.
# ============================================================================
.file_path_is_absolute <- function(x) {
    x <- as.character(x)
    grepl("^(/|~|[A-Za-z]:[\\/]|\\\\\\\\)", x)
}


#' Import pre-computed ratio / percentage dose-response tables
#'
#' Reads a directory of Excel files, each containing a pre-computed
#' dose-response table (for example an output from Prism or from a legacy
#' analysis script), and assembles them into a list whose shape matches
#' the return value of \code{\link{batch_ratio_analysis}} or
#' \code{\link{batch_viability_analysis}}.  The resulting object can be
#' passed directly to \code{\link{batch_drc_analysis}} and to every
#' downstream plotter or reporter in the package.
#'
#' @description
#' Use this function when the ratio and normalisation steps were already
#' performed outside \pkg{dosefitr}.  Each Excel file is expected to
#' contain exactly one dose-response table on the first sheet with:
#' \itemize{
#'   \item a log-concentration column (by default the first column of
#'     the sheet, whatever it is called; otherwise the column named by
#'     \code{log_conc_col_name}) containing the log-transformed
#'     concentrations.  NAs are allowed in padding rows, for example the
#'     top / bottom anchors used by Prism;
#'   \item one or more compound columns whose header follows the
#'     \code{construct:compound} pattern, with an optional replicate
#'     suffix \code{.2}, \code{.3}, ... appended to the compound token.
#' }
#'
#' @param directory Character string.  Directory containing the Excel
#'   files.  Defaults to \code{getwd()}, i.e. the user's current working
#'   directory.  Must exist.
#'
#' @param assay_source Character string.  Either \code{"nanobret"}
#'   (default) or \code{"viability"}.  Sets the top-level attribute read
#'   by \code{\link{batch_drc_analysis}} to choose assay-appropriate
#'   plausibility limits and downstream defaults.
#'
#' @param file_pattern Regular expression used to identify data files in
#'   \code{directory} (default \code{"\\\\.xlsx$"}).  Files whose names
#'   begin with \code{~$} (temporary Excel lock files) are always
#'   excluded.
#'
#' @param log_conc_col_name Character string or \code{NULL}.
#'   \itemize{
#'     \item \code{NULL} (default) - the log-concentration column is
#'       taken to be the first column of each sheet, whatever it is
#'       named.  The column is validated as numeric.
#'     \item A character string - the log-concentration column is looked
#'       up by this exact name.  If found and it is not already column 1,
#'       the returned \code{modified_ratio_table} is reordered so the
#'       log-concentration column is placed first (preserving the
#'       original order of the remaining compound columns).  If the name
#'       is not found in the sheet, an error is raised.
#'   }
#'
#' @param sheet Integer or character, length 1.  Sheet to read from each
#'   file, passed straight to \code{\link[openxlsx]{read.xlsx}}.  Default
#'   \code{1L}.
#'
#' @param plate_names Optional character vector, one element per file
#'   found.  If \code{NULL} (default), plates are named
#'   \code{plate_01}, \code{plate_02}, ... in sorted-filename order.
#'
#' @param file_map Optional explicit file list.  When supplied, this
#'   completely replaces the automatic discovery step and
#'   \code{file_pattern} is ignored.  Accepted forms:
#'   \itemize{
#'     \item a character vector of filenames or paths, e.g.
#'       \code{c("run_A.xlsx", "run_B.xlsx")};
#'     \item a list of length-1 character strings, matching the shape
#'       used by \code{\link{batch_ratio_analysis}}, e.g.
#'       \code{list("run_A.xlsx", "run_B.xlsx")}.  If the list is named,
#'       the names are \emph{ignored} - plate identifiers come from
#'       \code{plate_names}, or from the default
#'       \code{plate_01}, \code{plate_02}, ... sequence, in the order
#'       given.
#'   }
#'   Relative paths are resolved against \code{directory}; absolute
#'   paths are accepted as-is.  Missing files raise an error.
#'   Default \code{NULL} means auto-discover via \code{file_pattern}.
#'
#' @param verbose Logical.  If \code{TRUE} (default), prints a one-line
#'   summary and a QC caveat after a successful import.
#'
#' @details
#' \strong{Compatibility with the rest of the package}
#'
#' The returned object is deliberately shaped to mirror the slots that
#' \code{\link{batch_drc_analysis}}, \code{\link{scarab_table}} and
#' \code{\link{scarab_viability}} actually read.  In particular, the
#' modified ratio table is placed at
#' \code{result$plate_XX$result$modified_ratio_table} and the assay type
#' is set on \code{attr(result, "assay_source")}.  A sentinel attribute
#' \code{attr(result, "qc_available") <- FALSE} advertises that QC
#' inputs (Z-prime, positive / background means, Luciferase signal
#' comment) are unavailable because the raw plate reads are not present.
#'
#' \tabular{lll}{
#'   \strong{Downstream function} \tab \strong{Status} \tab \strong{Notes} \cr
#'   \code{batch_drc_analysis}                      \tab full    \tab assay auto-detected from attribute \cr
#'   \code{batch_save_all_drc_plots}                \tab full    \tab consumes DRC output only \cr
#'   \code{plot_multiple_compounds}                 \tab full    \tab consumes DRC output only \cr
#'   \code{compare_plates_drc}                      \tab full    \tab consumes DRC output only \cr
#'   \code{rout_outliers_batch}                     \tab partial \tab ROUT is applied to normalised percentages, not raw counts \cr
#'   \code{scarab_table} (NanoBRET)                 \tab partial \tab QC cells (Positive / Background means, Luciferase comment) will be blank \cr
#'   \code{scarab_viability}                        \tab full    \tab does not read \code{interval_means} \cr
#'   \code{merge_plate_replicates}                  \tab full    \tab operates on shape only \cr
#' }
#'
#' \strong{Column header rules}
#'
#' The log-concentration column is validated as numeric (NAs are allowed
#' in padding rows).  Every other column name must match the regular
#' expression \code{"^[^:]+:[^:]+(\\\\.\\\\d+)?$"}, i.e. one non-empty
#' non-\code{:} construct token, a single \code{:}, one non-empty
#' non-\code{:} compound token, and an optional \code{.<digits>}
#' replicate suffix.  Base-name collisions (columns that share
#' \code{construct:compound} after stripping any \code{.\\d+} suffix)
#' are allowed and are treated as replicates.  Empty column names,
#' whitespace-only names, and names containing multiple colons are
#' rejected with an informative error.
#'
#' @return A named list of plate results.  \code{names(result)} are the
#'   plate identifiers.  Each element has the same top-level fields as a
#'   single-plate slice of \code{\link{batch_ratio_analysis}}'s output.
#'   The list carries two attributes:
#'   \itemize{
#'     \item \code{attr(result, "assay_source")}: the string passed as
#'       \code{assay_source}.
#'     \item \code{attr(result, "qc_available")}: always \code{FALSE};
#'       downstream QC-consuming code can key off this sentinel.
#'   }
#'   The class is \code{c("dosefitr_batch_result", "list")}.
#'
#' @examples
#' stopifnot(requireNamespace("dosefitr", quietly = TRUE))
#' \donttest{
#' # Round-trip on a bundled NanoBRET plate --------------------------------
#' # 1. Run the full ratio pipeline on the bundled fixture
#' extdata <- system.file("extdata", package = "dosefitr")
#' work    <- tempfile("br_"); dir.create(work)
#' invisible(file.copy(
#'   list.files(extdata, "^nanobret_", full.names = TRUE),
#'   work, overwrite = TRUE
#' ))
#'
#' ratio_res <- batch_ratio_analysis(
#'   directory        = work,
#'   info_file        = "nanobret_info.xlsx",
#'   data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
#'   control_0perc    = "1",
#'   control_100perc  = "24",
#'   selected_columns = 1:24,
#'   generate_reports = FALSE,
#'   output_dir       = tempfile("br_out_"),
#'   verbose          = FALSE
#' )
#'
#' # 2. Export the modified ratio table to xlsx (mimicking Prism output)
#' shipped <- tempfile("shipped_"); dir.create(shipped)
#' openxlsx::write.xlsx(
#'   ratio_res$plate_01$result$modified_ratio_table,
#'   file = file.path(shipped, "plate_01.xlsx")
#' )
#'
#' # 3. Re-import it via batch_read_tables()
#' #    (assay_source defaults to "nanobret"; log_conc_col_name defaults
#' #     to NULL, i.e. use column 1 whatever it is called)
#' imported <- batch_read_tables(directory = shipped, verbose = FALSE)
#'
#' # 4. Feed straight into batch_drc_analysis
#' drc_res <- batch_drc_analysis(imported, model = "3pl",
#'                               generate_reports = FALSE, verbose = FALSE)
#' }
#'
#' @seealso
#'   \code{\link{batch_ratio_analysis}} for the full NanoBRET pipeline,
#'   \code{\link{batch_viability_analysis}} for the cell-viability
#'   pipeline, and \code{\link{batch_drc_analysis}} for the downstream
#'   dose-response fitting step.
#'
#' @export
batch_read_tables <- function(directory         = getwd(),
                              assay_source      = "nanobret",
                              file_pattern      = "\\.xlsx$",
                              log_conc_col_name = NULL,
                              sheet             = 1L,
                              plate_names       = NULL,
                              file_map          = NULL,
                              verbose           = TRUE) {

    # --- Validate top-level arguments --------------------------------------
    if (length(directory) != 1L || !is.character(directory) ||
        !nzchar(directory)) {
        stop("`directory` must be a single non-empty character string.",
             call. = FALSE)
    }
    if (!dir.exists(directory)) {
        stop("`directory` does not exist: ", directory, call. = FALSE)
    }

    if (length(assay_source) != 1L || !is.character(assay_source) ||
        !assay_source %in% c("nanobret", "viability")) {
        stop("`assay_source` must be exactly \"nanobret\" or \"viability\", ",
             "got: ", paste(deparse(assay_source), collapse = " "),
             call. = FALSE)
    }

    if (length(file_pattern) != 1L || !is.character(file_pattern)) {
        stop("`file_pattern` must be a single regular-expression string.",
             call. = FALSE)
    }
    if (!is.null(log_conc_col_name)) {
        if (length(log_conc_col_name) != 1L ||
            !is.character(log_conc_col_name) ||
            !nzchar(log_conc_col_name)) {
            stop("`log_conc_col_name` must be NULL or a single non-empty ",
                 "character string.", call. = FALSE)
        }
    }
    if (length(sheet) != 1L ||
        !(is.numeric(sheet) || is.character(sheet))) {
        stop("`sheet` must be a single numeric or character value.",
             call. = FALSE)
    }
    if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
        stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
    }

    # --- Normalise file_map to a character vector (or NULL) ----------------
    if (!is.null(file_map)) {
        if (is.list(file_map)) {
            # Reject nested / non-scalar-string entries early with a
            # crisp message rather than letting them fall through.
            ok <- vapply(file_map,
                         function(x) is.character(x) && length(x) == 1L &&
                                     !is.na(x) && nzchar(x),
                         logical(1L))
            if (!all(ok)) {
                stop("`file_map` list entries must each be a single ",
                     "non-empty character string.", call. = FALSE)
            }
            file_map <- unlist(file_map, use.names = FALSE)
        }
        if (!is.character(file_map) || length(file_map) == 0L ||
            any(is.na(file_map)) || any(!nzchar(file_map))) {
            stop("`file_map` must be a non-empty character vector or a ",
                 "list of single non-empty character strings.",
                 call. = FALSE)
        }
        if (anyDuplicated(file_map) > 0L) {
            stop("`file_map` contains duplicate entries.", call. = FALSE)
        }
    }

    # --- Discover files ----------------------------------------------------
    was_file_map <- !is.null(file_map)
    if (is.null(file_map)) {
        all_files <- list.files(directory,
                                pattern    = file_pattern,
                                full.names = TRUE)
        # Exclude Excel temporary lock files
        all_files <- all_files[!grepl("^~\\$", basename(all_files))]
        all_files <- sort(all_files)

        if (length(all_files) == 0L) {
            stop("No files matching `file_pattern` (", file_pattern,
                 ") found in directory: ", directory, call. = FALSE)
        }
    } else {
        # file_map: absolute paths kept as-is; relative paths joined to
        # `directory`.  Order is preserved (NOT sorted) so callers keep
        # control of plate ordering.
        is_abs   <- .file_path_is_absolute(file_map)
        resolved <- ifelse(is_abs, file_map, file.path(directory, file_map))
        missing  <- resolved[!file.exists(resolved)]
        if (length(missing) > 0L) {
            stop("`file_map` refers to file(s) that do not exist:\n  ",
                 paste0("- ", missing, collapse = "\n  "),
                 call. = FALSE)
        }
        all_files <- resolved
    }

    # --- Resolve plate names ----------------------------------------------
    n_files <- length(all_files)
    if (is.null(plate_names)) {
        plate_names <- sprintf("plate_%02d", seq_len(n_files))
    } else {
        if (!is.character(plate_names) || length(plate_names) != n_files) {
            stop("`plate_names` must be a character vector with one entry ",
                 "per file found (found ", n_files, " file(s)).",
                 call. = FALSE)
        }
        if (anyDuplicated(plate_names) > 0L) {
            stop("`plate_names` entries must be unique.", call. = FALSE)
        }
    }

    # --- Constants ---------------------------------------------------------
    COMPOUND_COL_RX <- "^[^:]+:[^:]+(\\.\\d+)?$"

    # --- Read + validate each file ----------------------------------------
    plate_list <- vector("list", n_files)
    names(plate_list) <- plate_names

    for (i in seq_len(n_files)) {
        f    <- all_files[[i]]
        fbn  <- basename(f)

        tbl <- tryCatch(
            openxlsx::read.xlsx(f, sheet = sheet, check.names = FALSE),
            error = function(e) {
                stop("Could not read sheet ",
                     if (is.numeric(sheet)) sheet
                     else paste0("\"", sheet, "\""),
                     " of file '", fbn, "': ", conditionMessage(e),
                     call. = FALSE)
            }
        )
        if (!is.data.frame(tbl) || nrow(tbl) == 0L || ncol(tbl) < 2L) {
            stop("File '", fbn, "' does not contain a usable table on sheet ",
                 if (is.numeric(sheet)) sheet
                 else paste0("\"", sheet, "\""),
                 " (need >=1 dose column and >=1 compound column).",
                 call. = FALSE)
        }

        col_names <- colnames(tbl)

        # 1. Resolve the log-concentration column.
        #    - log_conc_col_name = NULL: use column 1 as-is.
        #    - log_conc_col_name = "<string>": look up by name; if it is
        #      not the first column, reorder so it becomes the first.
        if (is.null(log_conc_col_name)) {
            log_idx <- 1L
        } else {
            hits <- which(col_names == log_conc_col_name)
            if (length(hits) == 0L) {
                stop("File '", fbn, "' has no column named '",
                     log_conc_col_name, "'. Available columns: ",
                     paste0("\"", col_names, "\"", collapse = ", "), ".",
                     call. = FALSE)
            }
            if (length(hits) > 1L) {
                stop("File '", fbn, "' has multiple columns named '",
                     log_conc_col_name, "' (found ", length(hits),
                     " times); expected exactly one.",
                     call. = FALSE)
            }
            log_idx <- hits[[1L]]
            if (log_idx != 1L) {
                # Reorder so the log-conc column becomes column 1.
                new_order <- c(log_idx, setdiff(seq_along(col_names), log_idx))
                tbl <- tbl[, new_order, drop = FALSE]
                col_names <- colnames(tbl)
            }
        }
        log_col_actual <- col_names[[1L]]

        # 2. Log-concentration column must be numeric (NAs allowed).
        first_col <- tbl[[1L]]
        if (!is.numeric(first_col)) {
            stop("File '", fbn, "' has non-numeric values in the ",
                 "log-concentration column ('", log_col_actual, "'). ",
                 "All dose values must be numeric (NAs are allowed in ",
                 "padding rows).",
                 call. = FALSE)
        }

        # 3. All remaining columns must match COMPOUND_COL_RX.
        cmpd_names <- col_names[-1L]
        if (length(cmpd_names) == 0L) {
            stop("File '", fbn, "' has no compound columns; ",
                 "at least one 'construct:compound' column is required.",
                 call. = FALSE)
        }
        bad <- cmpd_names[!grepl(COMPOUND_COL_RX, cmpd_names)]
        if (length(bad) > 0L) {
            stop(
                "File '", fbn,
                "' has invalid column name(s): ",
                paste0("\"", bad, "\"", collapse = ", "),
                ".\n",
                "Column names must follow the 'construct:compound' pattern ",
                "with an\noptional replicate suffix. Examples:\n",
                "  - LRRK2:DDD02453312            (first replicate, no suffix)\n",
                "  - LRRK2:DDD02453312.2          (second replicate)\n",
                "  - LRRK2:DDD02453312.3          (third replicate)\n",
                "  - LRRK2:DDD02453312.4          (fourth replicate, etc.)",
                call. = FALSE
            )
        }

        # 4. No entirely-NA compound columns.
        na_cols <- cmpd_names[vapply(cmpd_names,
                                     function(nm) all(is.na(tbl[[nm]])),
                                     logical(1L))]
        if (length(na_cols) > 0L) {
            stop("File '", fbn, "' has compound column(s) that are ",
                 "entirely NA: ",
                 paste0("\"", na_cols, "\"", collapse = ", "), ".",
                 call. = FALSE)
        }

        # --- Assemble per-plate list -------------------------------------
        # `sheet` may be an integer or a sheet name.  Downstream code only
        # peeks at this field for messages, so store the raw value.
        selected_cols <- seq_len(ncol(tbl) - 1L)

        plate_list[[i]] <- list(
            data_file        = normalizePath(f, mustWork = TRUE),
            sheet            = sheet,
            function_version = "imported",
            control_0perc    = NA_character_,
            control_100perc  = NA_character_,
            selected_columns = selected_cols,
            result = list(
                original_ratio_table   = NULL,
                modified_ratio_table   = tbl,
                interval_means         = NULL,
                general_means          = NULL,
                construct_intervals    = NULL,
                selected_columns_info  = list(
                    original_columns      = NULL,
                    selected_columns      = NULL,
                    was_selection_applied = FALSE
                )
            )
        )
    }

    # --- Attach attributes / class + optional summary ---------------------
    attr(plate_list, "assay_source") <- assay_source
    attr(plate_list, "qc_available") <- FALSE
    class(plate_list) <- c("dosefitr_batch_result", "list")

    if (isTRUE(verbose)) {
        first_tbl <- plate_list[[1L]]$result$modified_ratio_table
        first_col <- first_tbl[[1L]]
        log_range <- suppressWarnings(range(first_col, na.rm = TRUE))
        n_cmpd    <- ncol(first_tbl) - 1L

        message(sprintf(
            "batch_read_tables(): imported %d plate(s) with %d compound(s) each.",
            n_files, n_cmpd
        ))
        message(sprintf("   assay_source    = '%s'", assay_source))
        message(sprintf("   discovery       = %s",
                        if (was_file_map) "file_map (explicit)"
                        else sprintf("file_pattern '%s'", file_pattern)))
        message(sprintf("   log_conc column = '%s'",
                        colnames(first_tbl)[[1L]]))
        message(sprintf("   plate names     = %s",
                        paste(names(plate_list), collapse = ", ")))
        message(sprintf("   log[M] range    = %.2f to %.2f",
                        log_range[[1L]], log_range[[2L]]))
        message("")
        message(
            "Note: Z-prime, signal-to-background, positive/background means, ",
            "and the\nLuciferase signal comment are NOT available for ",
            "pre-computed ratio tables.\nWhen you later call scarab_table(), ",
            "those cells will be blank.\nscarab_viability() and every ",
            "plotting/curve-fitting function are unaffected."
        )
    }

    plate_list
}
