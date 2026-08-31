#' Reorganize 384-well plates by compound
#'
#' Rebuilds raw single-channel 384-well plate-reader workbooks so that each
#' output workbook contains all available rows for one compound. The mapping
#' between source plate rows, targets, and compounds is read from a multi-sheet
#' information workbook. A matching information workbook for the reorganized
#' plates and a row-level provenance workbook are written alongside the plate
#' files.
#'
#' This is useful when an experiment was acquired as one plate per target (for
#' example, one cell line per plate), but downstream analysis should operate on
#' one plate per compound across targets.
#'
#' @param directory A single character string giving the directory containing
#'   the raw plate workbooks and \code{info_file}. Defaults to the current working
#'   directory.
#' @param info_file A single character string naming the multi-sheet Excel
#'   information workbook. Relative paths are resolved against \code{directory}.
#'   Each sheet must be named for its source plate (for example, \code{plate_01})
#'   and its first four columns must contain log concentration, plate row,
#'   target, and compound, in that order.
#' @param data_pattern A regular expression used to select raw plate files.
#'   Matching filenames must end in a numeric plate identifier immediately
#'   before \code{.xlsx}; this identifier is matched to information-sheet names.
#'   Ignored when \code{file_map} is supplied.
#' @param file_map An optional named character vector or named list mapping
#'   information-sheet names to raw workbook paths, for example
#'   \code{c(plate_01 = "A375.xlsx", plate_02 = "HT29.xlsx")}. Relative paths
#'   are resolved against \code{directory}. When supplied, it completely
#'   replaces automatic discovery through \code{data_pattern}.
#' @param data_sheet A single numeric or character value identifying the sheet
#'   to read from each raw plate workbook.
#' @param output_dir A single character string giving the output directory.
#'   Relative paths are resolved against \code{directory}.
#' @param overwrite Logical. If \code{FALSE}, the function stops
#'   before writing when any target file already exists. Set to \code{TRUE} (the default) to
#'   replace existing generated files.
#' @param verbose Logical. If \code{TRUE} (the default), report a short completion
#'   message.
#'
#' @details
#' The function expects 384-well plates with rows \code{A} through \code{P} and 24 data
#' columns. In each raw workbook, the first column identifies the plate row;
#' the next 24 columns are copied to the reorganized workbook. Instrument
#' metadata above the data block is allowed and is ignored.
#'
#' Each output plate starts at the first worksheet row: row 1 contains the
#' column numbers 1 through 24, and rows 2 through 17 contain plate rows
#' \code{A} through \code{P}. No empty metadata rows are added above the table.
#'
#' The generated \code{reorganization_provenance.xlsx} workbook contains one
#' row per transferred plate row. It records the source file, source plate,
#' source row, target, compound, output file, output plate, and output row.
#'
#' The first column of every information sheet must contain the same 16-value
#' log-concentration template. Information records are processed in workbook
#' sheet order and then row order. Because a 384-well plate has 16 rows, a
#' compound may have at most 16 source records.
#'
#' All inputs are validated before output files are created. In particular,
#' missing plate references, missing source rows, duplicated plate identifiers,
#' non-numeric measurements, unsafe filename collisions, and accidental input
#' overwrites produce errors.
#'
#' @return Invisibly, a data frame with one row per generated compound plate
#'   and columns \code{plate}, \code{compound}, \code{source_rows}, and
#'   \code{file}. The paths to the generated information and provenance
#'   workbooks are stored in the \code{info_file} and \code{provenance_file}
#'   attributes, respectively.
#'
#' @examples
#' \donttest{
#' extdata <- system.file("extdata", package = "dosefitr")
#' example_dir <- tempfile("dosefitr_reorganize_")
#' dir.create(example_dir)
#' fixture_files <- list.files(
#'   extdata,
#'   pattern = "^viability_(info|plate_[0-9]+)\\.xlsx$",
#'   full.names = TRUE
#' )
#' stopifnot(all(file.copy(fixture_files, example_dir)))
#'
#' manifest <- batch_reorganize_plates(
#'   directory = example_dir,
#'   info_file = "viability_info.xlsx",
#'   file_map = c(
#'     plate_01 = "viability_plate_01.xlsx",
#'     plate_02 = "viability_plate_02.xlsx"
#'   ),
#'   output_dir = "reorganized",
#'   verbose = FALSE
#' )
#'
#' utils::head(manifest, 2L)
#' attr(manifest, "info_file")
#' attr(manifest, "provenance_file")
#' unlink(example_dir, recursive = TRUE)
#' }
#'
#' @seealso \code{\link{batch_viability_analysis}} for processing the
#'   reorganized raw plates.
#' @export
batch_reorganize_plates <- function(
    directory = getwd(),
    info_file = "info_tables.xlsx",
    data_pattern = "_[0-9]{2}\\.xlsx$",
    file_map = NULL,
    data_sheet = 1L,
    output_dir = "reorganized_plates",
    overwrite = TRUE,
    verbose = TRUE) {

    .reorg_assert_string(directory, "directory")
    .reorg_assert_string(info_file, "info_file")
    .reorg_assert_string(data_pattern, "data_pattern")
    .reorg_assert_string(output_dir, "output_dir")
    .reorg_assert_flag(overwrite, "overwrite")
    .reorg_assert_flag(verbose, "verbose")

    if (length(data_sheet) != 1L || is.na(data_sheet) ||
        !(is.numeric(data_sheet) || is.character(data_sheet))) {
        stop("`data_sheet` must be one non-missing numeric or character value.",
             call. = FALSE)
    }
    if (!dir.exists(directory)) {
        stop("`directory` does not exist: ", directory, call. = FALSE)
    }

    directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
    info_path <- .reorg_resolve_path(info_file, directory)
    output_path <- .reorg_resolve_path(output_dir, directory)

    if (!file.exists(info_path)) {
        stop("`info_file` does not exist: ", info_path, call. = FALSE)
    }
    if (file.exists(output_path) && !dir.exists(output_path)) {
        stop("`output_dir` exists but is not a directory: ", output_path,
             call. = FALSE)
    }

    info_sheets <- tryCatch(
        openxlsx::getSheetNames(info_path),
        error = function(e) {
            stop("Could not inspect `info_file`: ", conditionMessage(e),
                 call. = FALSE)
        }
    )
    if (length(info_sheets) == 0L) {
        stop("`info_file` contains no worksheets.", call. = FALSE)
    }
    if (anyDuplicated(tolower(info_sheets)) > 0L) {
        stop("Worksheet names in `info_file` must be unique when case is ignored.",
             call. = FALSE)
    }

    info_tables <- lapply(info_sheets, function(sheet_name) {
        tryCatch(
            openxlsx::read.xlsx(
                info_path,
                sheet = sheet_name,
                check.names = FALSE,
                skipEmptyRows = FALSE,
                skipEmptyCols = FALSE
            ),
            error = function(e) {
                stop("Could not read sheet '", sheet_name, "' from `info_file`: ",
                     conditionMessage(e), call. = FALSE)
            }
        )
    })

    too_narrow <- which(vapply(info_tables, ncol, integer(1L)) < 4L)
    if (length(too_narrow) > 0L) {
        stop("Every information sheet must contain at least four columns; ",
             "invalid sheet(s): ",
             paste(info_sheets[too_narrow], collapse = ", "), ".",
             call. = FALSE)
    }

    log_template <- info_tables[[1L]][[1L]]
    if (length(log_template) != 16L) {
        stop("The first column of each information sheet must contain exactly ",
             "16 values (one for each plate row); found ",
             length(log_template), " in sheet '", info_sheets[[1L]], "'.",
             call. = FALSE)
    }
    for (i in seq_along(info_tables)) {
        candidate <- info_tables[[i]][[1L]]
        if (!isTRUE(all.equal(candidate, log_template,
                              check.attributes = FALSE))) {
            stop("The first-column log-concentration values differ between ",
                 "sheet '", info_sheets[[1L]], "' and sheet '",
                 info_sheets[[i]], "'.", call. = FALSE)
        }
    }

    info_records <- vector("list", length(info_tables))
    for (i in seq_along(info_tables)) {
        sheet_name <- info_sheets[[i]]
        tab <- info_tables[[i]][, seq_len(4L), drop = FALSE]
        names(tab) <- c("log_value", "plate_row", "target", "compound")
        tab$plate_row <- toupper(trimws(as.character(tab$plate_row)))
        tab$target <- as.character(tab$target)
        tab$compound <- trimws(as.character(tab$compound))

        has_compound <- !is.na(tab$compound) & nzchar(tab$compound)
        tab <- tab[has_compound, , drop = FALSE]
        if (nrow(tab) > 0L) {
            invalid_rows <- is.na(tab$plate_row) |
                !tab$plate_row %in% LETTERS[seq_len(16L)]
            if (any(invalid_rows)) {
                stop("Sheet '", sheet_name,
                     "' contains invalid `Plate_Row` value(s): ",
                     paste(unique(tab$plate_row[invalid_rows]), collapse = ", "),
                     ". Expected A through P.", call. = FALSE)
            }
        }
        tab$source_plate <- sheet_name
        tab$source_plate_key <- tolower(sheet_name)
        info_records[[i]] <- tab
    }
    compound_info <- do.call(rbind, info_records)
    rownames(compound_info) <- NULL
    if (nrow(compound_info) == 0L) {
        stop("No non-empty compound records were found in `info_file`.",
             call. = FALSE)
    }

    if (is.null(file_map)) {
        raw_files <- sort(list.files(
            directory,
            pattern = data_pattern,
            full.names = TRUE,
            ignore.case = TRUE
        ))
        raw_files <- raw_files[!grepl("^~\\$", basename(raw_files))]
        raw_files <- raw_files[
            .reorg_path_key(raw_files) != .reorg_path_key(info_path)
        ]
        if (length(raw_files) == 0L) {
            stop("No raw plate files matching `data_pattern` were found in: ",
                 directory, call. = FALSE)
        }

        plate_ids <- vapply(raw_files, .reorg_plate_id, character(1L))
        raw_plate_names <- paste0("plate_", plate_ids)
        raw_plate_keys <- tolower(raw_plate_names)
        if (anyDuplicated(raw_plate_keys) > 0L) {
            duplicated_ids <- unique(
                raw_plate_names[duplicated(raw_plate_keys)]
            )
            stop("More than one raw file resolves to plate identifier(s): ",
                 paste(duplicated_ids, collapse = ", "), ".", call. = FALSE)
        }
    } else {
        resolved_map <- .reorg_resolve_file_map(file_map, directory)
        raw_files <- unname(resolved_map)
        raw_plate_names <- names(resolved_map)
        raw_plate_keys <- tolower(raw_plate_names)

        unknown_keys <- setdiff(raw_plate_keys, tolower(info_sheets))
        if (length(unknown_keys) > 0L) {
            unknown_names <- raw_plate_names[raw_plate_keys %in% unknown_keys]
            stop("`file_map` contains plate name(s) not found in `info_file`: ",
                 paste(unknown_names, collapse = ", "), ".", call. = FALSE)
        }
    }

    raw_data <- vector("list", length(raw_files))
    names(raw_data) <- raw_plate_keys
    raw_file_map <- stats::setNames(raw_files, raw_plate_keys)
    for (i in seq_along(raw_files)) {
        raw_data[[i]] <- .reorg_read_raw_plate(
            raw_files[[i]], data_sheet, raw_plate_names[[i]]
        )
    }

    missing_plate_keys <- setdiff(unique(compound_info$source_plate_key),
                                  names(raw_data))
    if (length(missing_plate_keys) > 0L) {
        missing_names <- unique(compound_info$source_plate[
            compound_info$source_plate_key %in% missing_plate_keys
        ])
        stop("No matching raw workbook was found for information sheet(s): ",
             paste(missing_names, collapse = ", "), ".", call. = FALSE)
    }

    for (i in seq_len(nrow(compound_info))) {
        plate_key <- compound_info$source_plate_key[[i]]
        row_name <- compound_info$plate_row[[i]]
        if (!row_name %in% names(raw_data[[plate_key]])) {
            stop("Source row '", row_name, "' referenced by compound '",
                 compound_info$compound[[i]], "' is missing from raw plate '",
                 compound_info$source_plate[[i]], "'.", call. = FALSE)
        }
    }

    compounds <- unique(compound_info$compound)
    counts <- vapply(compounds, function(compound_name) {
        sum(compound_info$compound == compound_name)
    }, integer(1L))
    if (any(counts > 16L)) {
        bad <- paste0(compounds[counts > 16L], " (", counts[counts > 16L], ")")
        stop("A compound cannot occupy more than 16 output rows; found: ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
    }

    new_plate_ids <- sprintf("%02d", seq_along(compounds))
    output_plate_names <- paste0("plate_", new_plate_ids)
    safe_compounds <- vapply(compounds, .reorg_safe_filename, character(1L))
    output_basenames <- paste0(safe_compounds, "_", new_plate_ids, ".xlsx")
    if (anyDuplicated(tolower(output_basenames)) > 0L) {
        collision <- unique(output_basenames[duplicated(tolower(output_basenames))])
        stop("Compound names produce colliding output filenames: ",
             paste(collision, collapse = ", "), ".", call. = FALSE)
    }

    output_files <- file.path(output_path, output_basenames)
    output_info <- file.path(output_path, "info_tables.xlsx")
    output_provenance <- file.path(
        output_path, "reorganization_provenance.xlsx"
    )
    input_keys <- .reorg_path_key(c(raw_files, info_path))
    output_keys <- .reorg_path_key(
        c(output_files, output_info, output_provenance)
    )
    if (any(output_keys %in% input_keys)) {
        stop("`output_dir` would overwrite one or more input workbooks. ",
             "Choose a separate output directory.", call. = FALSE)
    }

    target_files <- c(output_files, output_info, output_provenance)
    existing <- target_files[file.exists(target_files)]
    if (!overwrite && length(existing) > 0L) {
        stop("Output file(s) already exist. Use `overwrite = TRUE` to replace ",
             "them:\n  ", paste(existing, collapse = "\n  "), call. = FALSE)
    }

    new_plates <- vector("list", length(compounds))
    new_info <- vector("list", length(compounds))
    names(new_info) <- output_plate_names
    provenance_rows <- vector("list", nrow(compound_info))
    provenance_index <- 0L

    for (i in seq_along(compounds)) {
        compound_name <- compounds[[i]]
        records <- compound_info[
            compound_info$compound == compound_name, , drop = FALSE
        ]

        new_plate <- as.data.frame(
            matrix(NA, nrow = 17L, ncol = 25L),
            stringsAsFactors = FALSE
        )
        new_plate[1L, 2:25] <- seq_len(24L)

        new_info_table <- data.frame(
            "log(inhibitor) [M]" = log_template,
            Plate_Row = LETTERS[seq_len(16L)],
            Target = rep(NA_character_, 16L),
            Compound = rep(compound_name, 16L),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        for (j in seq_len(nrow(records))) {
            destination_row <- LETTERS[[j]]
            source_values <- raw_data[[records$source_plate_key[[j]]]][[
                records$plate_row[[j]]
            ]]
            new_plate[1L + j, 1L] <- destination_row
            new_plate[1L + j, 2:25] <- source_values
            new_info_table$Target[[j]] <- records$target[[j]]

            provenance_index <- provenance_index + 1L
            provenance_rows[[provenance_index]] <- data.frame(
                Source_File = basename(raw_file_map[[
                    records$source_plate_key[[j]]
                ]]),
                Source_Plate = records$source_plate[[j]],
                Source_Row = records$plate_row[[j]],
                Target = records$target[[j]],
                Compound = compound_name,
                Output_File = output_basenames[[i]],
                Output_Plate = output_plate_names[[i]],
                Output_Row = destination_row,
                stringsAsFactors = FALSE
            )
        }

        for (j in seq_len(16L)) {
            if (is.na(new_plate[1L + j, 1L])) {
                new_plate[1L + j, 1L] <- LETTERS[[j]]
            }
        }

        new_plates[[i]] <- new_plate
        new_info[[i]] <- new_info_table
    }
    provenance <- do.call(rbind, provenance_rows)
    rownames(provenance) <- NULL

    if (!dir.exists(output_path)) {
        if (!dir.create(output_path, recursive = TRUE, showWarnings = FALSE)) {
            stop("Could not create `output_dir`: ", output_path, call. = FALSE)
        }
    }

    for (i in seq_along(new_plates)) {
        tryCatch(
            openxlsx::write.xlsx(
                new_plates[[i]],
                file = output_files[[i]],
                colNames = FALSE,
                rowNames = FALSE,
                overwrite = TRUE
            ),
            error = function(e) {
                stop("Could not write output plate '", output_files[[i]],
                     "': ", conditionMessage(e), call. = FALSE)
            }
        )
    }
    tryCatch(
        openxlsx::write.xlsx(
            new_info,
            file = output_info,
            colNames = TRUE,
            rowNames = FALSE,
            overwrite = TRUE
        ),
        error = function(e) {
            stop("Could not write the reorganized information workbook: ",
                 conditionMessage(e), call. = FALSE)
        }
    )
    tryCatch(
        openxlsx::write.xlsx(
            list(Provenance = provenance),
            file = output_provenance,
            asTable = TRUE,
            overwrite = TRUE
        ),
        error = function(e) {
            stop("Could not write the reorganization provenance workbook: ",
                 conditionMessage(e), call. = FALSE)
        }
    )

    manifest <- data.frame(
        plate = output_plate_names,
        compound = compounds,
        source_rows = unname(counts),
        file = normalizePath(output_files, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
    )
    attr(manifest, "info_file") <- normalizePath(
        output_info, winslash = "/", mustWork = TRUE
    )
    attr(manifest, "provenance_file") <- normalizePath(
        output_provenance, winslash = "/", mustWork = TRUE
    )

    if (verbose) {
        message("Reorganized ", length(compounds), " compound plate(s) in '",
                output_path, "'.")
    }
    invisible(manifest)
}


.reorg_assert_string <- function(x, argument) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop("`", argument, "` must be one non-empty character string.",
             call. = FALSE)
    }
}


.reorg_assert_flag <- function(x, argument) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
}


.reorg_is_absolute <- function(path) {
    grepl("^(/|~|[A-Za-z]:[\\\\/]|\\\\\\\\)", path)
}


.reorg_resolve_path <- function(path, directory) {
    resolved <- if (.reorg_is_absolute(path)) path.expand(path) else
        file.path(directory, path)
    normalizePath(resolved, winslash = "/", mustWork = FALSE)
}


.reorg_path_key <- function(path) {
    tolower(normalizePath(path, winslash = "/", mustWork = FALSE))
}


.reorg_resolve_file_map <- function(file_map, directory) {
    if (is.list(file_map)) {
        valid_entries <- vapply(file_map, function(x) {
            is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
        }, logical(1L))
        if (!all(valid_entries)) {
            stop("Every `file_map` list entry must be one non-empty character ",
                 "string.", call. = FALSE)
        }
        map_names <- names(file_map)
        file_map <- unlist(file_map, use.names = FALSE)
        names(file_map) <- map_names
    }

    if (!is.character(file_map) || length(file_map) == 0L ||
        any(is.na(file_map)) || any(!nzchar(file_map))) {
        stop("`file_map` must be a non-empty named character vector or named ",
             "list.", call. = FALSE)
    }
    map_names <- names(file_map)
    if (is.null(map_names) || any(is.na(map_names)) ||
        any(!nzchar(map_names))) {
        stop("Every `file_map` entry must have a non-empty plate name, for ",
             "example `plate_01`.", call. = FALSE)
    }
    if (anyDuplicated(tolower(map_names)) > 0L) {
        stop("`file_map` plate names must be unique when case is ignored.",
             call. = FALSE)
    }

    resolved <- vapply(file_map, .reorg_resolve_path, character(1L),
                       directory = directory)
    missing <- resolved[!file.exists(resolved)]
    if (length(missing) > 0L) {
        stop("`file_map` refers to file(s) that do not exist:\n  ",
             paste(missing, collapse = "\n  "), call. = FALSE)
    }
    if (anyDuplicated(.reorg_path_key(resolved)) > 0L) {
        stop("`file_map` cannot map the same raw workbook more than once.",
             call. = FALSE)
    }
    stats::setNames(unname(resolved), map_names)
}


.reorg_plate_id <- function(path) {
    match <- regmatches(
        basename(path),
        regexec("([0-9]+)(?=\\.xlsx$)", basename(path),
                ignore.case = TRUE, perl = TRUE)
    )[[1L]]
    if (length(match) < 2L) {
        stop("Could not extract a trailing numeric plate identifier from file: ",
             basename(path), call. = FALSE)
    }
    match[[2L]]
}


.reorg_read_raw_plate <- function(path, sheet, plate_name) {
    raw <- tryCatch(
        openxlsx::read.xlsx(
            path,
            sheet = sheet,
            colNames = FALSE,
            check.names = FALSE,
            skipEmptyRows = FALSE,
            skipEmptyCols = FALSE
        ),
        error = function(e) {
            stop("Could not read raw plate '", basename(path), "': ",
                 conditionMessage(e), call. = FALSE)
        }
    )
    if (!is.data.frame(raw) || ncol(raw) < 25L) {
        stop("Raw plate '", basename(path),
             "' must contain a row-label column followed by at least 24 data ",
             "columns.", call. = FALSE)
    }

    row_labels <- toupper(trimws(as.character(raw[[1L]])))
    keep <- !is.na(row_labels) & row_labels %in% LETTERS[seq_len(16L)]
    raw <- raw[keep, , drop = FALSE]
    row_labels <- row_labels[keep]
    if (length(row_labels) == 0L) {
        stop("Raw plate '", basename(path),
             "' contains no rows labelled A through P in its first column.",
             call. = FALSE)
    }
    if (anyDuplicated(row_labels) > 0L) {
        duplicate_rows <- unique(row_labels[duplicated(row_labels)])
        stop("Raw plate '", basename(path),
             "' contains duplicated data row(s): ",
             paste(duplicate_rows, collapse = ", "), ".", call. = FALSE)
    }

    values <- vector("list", length(row_labels))
    names(values) <- row_labels
    for (i in seq_along(row_labels)) {
        cells <- unlist(raw[i, 2:25, drop = FALSE], use.names = FALSE)
        numeric_cells <- suppressWarnings(as.numeric(cells))
        nonempty <- !is.na(cells) & nzchar(trimws(as.character(cells)))
        if (any(nonempty & is.na(numeric_cells))) {
            stop("Raw plate '", basename(path), "' row '", row_labels[[i]],
                 "' contains non-numeric measurements.", call. = FALSE)
        }
        if (any(!is.na(numeric_cells) & !is.finite(numeric_cells))) {
            stop("Raw plate '", basename(path), "' row '", row_labels[[i]],
                 "' contains non-finite measurements.", call. = FALSE)
        }
        values[[i]] <- numeric_cells
    }
    values
}


.reorg_safe_filename <- function(compound) {
    safe <- gsub("[/\\\\:*?\"<>|]", "_", compound)
    safe <- gsub("[[:cntrl:]]", "_", safe)
    safe <- sub("[ .]+$", "", safe)
    if (!nzchar(safe)) {
        safe <- "compound"
    }

    stem <- toupper(sub("\\..*$", "", safe))
    reserved <- c("CON", "PRN", "AUX", "NUL",
                  paste0("COM", seq_len(9L)), paste0("LPT", seq_len(9L)))
    if (stem %in% reserved) {
        safe <- paste0("_", safe)
    }
    safe
}
