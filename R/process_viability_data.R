#' Process viability data from multi-well plate experiments
#'
#' This function processes raw viability data (e.g., from 384-well or 96-well plates)
#' and converts it into a structured format suitable for downstream analysis.
#' It supports automatic table detection, control normalization, replicate splitting,
#' and annotation using an optional metadata table.
#'
#' @description
#' The function detects a plate-like structure with rows labeled A-P and columns 1-24,
#' extracts the viability data, optionally filters and normalizes values, applies
#' control-based corrections, and reshapes the data for analysis.
#'
#' It is particularly useful for dose-response experiments, high-throughput screening,
#' or viability assays using fluorescence or luminescence readouts.
#'
#' @param data A data.frame or matrix containing raw plate data.
#'   Must contain at least 16 rows corresponding to plate rows A-P.
#'
#' @param control_0perc Numeric or character. Column representing the 0% control
#'   (e.g., background). If numeric, refers to column index (1-24).
#'
#' @param control_100perc Numeric or character. Column representing the 100% control
#'   (e.g., untreated or positive control). If numeric, refers to column index (1-24).
#'
#' @param split_replicates Logical. If TRUE (default), splits experimental concentrations
#'   into two technical replicates.
#'
#' @param info_table Optional data.frame containing metadata. Must contain at least
#'   four columns:
#'   \describe{
#'     \item{Column 1}{log(inhibitor) values}
#'     \item{Column 2}{Plate row identifiers (A-P)}
#'     \item{Column 3}{Construct name}
#'     \item{Column 4}{Compound name}
#'   }
#'
#' @param selected_columns Optional numeric vector specifying which data columns (1-24)
#'   to include. If NULL, all columns are used.
#'
#' @param low_value_threshold Numeric. Values below this threshold are replaced with NA.
#'   Default is 0.
#'
#' @param verbose Logical. If TRUE (default), prints progress messages and warnings.
#'
#' @param apply_control_means Logical. If TRUE (default), replaces control values with
#'   construct-specific means based on \code{info_table}.
#'
#' @param auto_detect Logical. If TRUE (default), automatically detects the plate
#'   structure within the input data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Detects plate layout (rows A-P and columns 1-24)
#'   \item Extracts viability values
#'   \item Filters low or zero values
#'   \item Maps control columns (0% and 100%)
#'   \item Optionally applies construct-specific control averaging
#'   \item Reorders columns to place controls at beginning and end
#'   \item Transposes the dataset
#'   \item Optionally renames columns using \code{info_table}
#'   \item Splits technical replicates (if enabled)
#'   \item Adds log(inhibitor) values (if provided)
#' }
#'
#' Column indices (1-24) always refer to the original plate layout,
#' regardless of subsetting.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{original_table}{Numeric data.frame of extracted viability values}
#'   \item{modified_table}{Processed data.frame ready for analysis}
#'   \item{processing_info}{List containing intermediate data and metadata}
#'   \item{selected_columns_info}{Information about selected columns}
#'   \item{version}{Function version (4.2)}
#'   \item{data_type}{Type of data ("viability")}
#'   \item{auto_detect}{Logical indicating if auto-detection was used}
#'   \item{behavior_mode}{Processing mode ("original")}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- process_viability_data(data = raw_data)
#'
#' # With control columns
#' result <- process_viability_data(
#'   data = raw_data,
#'   control_0perc = 1,
#'   control_100perc = 24
#' )
#'
#' # With metadata and column selection
#' result <- process_viability_data(
#'   data = raw_data,
#'   control_0perc = 1,
#'   control_100perc = 24,
#'   selected_columns = 2:23,
#'   info_table = metadata
#' )
#' }
#'
#' @seealso
#' \code{\link{t}}, \code{\link{apply}}
#'
#'
#' @export

process_viability_data <- function(data,
                                   control_0perc = NULL, control_100perc = NULL,
                                   split_replicates = TRUE, info_table = NULL,
                                   selected_columns = NULL,
                                   low_value_threshold = 0,
                                   verbose = TRUE,
                                   apply_control_means = TRUE,
                                   auto_detect = TRUE) {

  # --- SPECIFIC PATTERN DETECTION FUNCTION ---
  detect_specific_pattern <- function(data) {
    if (verbose) {
      message("Looking for specific pattern: A-P rows and 1-24 columns...")
    }

    # Convert to data.frame if not already
    if (!is.data.frame(data)) {
      data <- as.data.frame(data, stringsAsFactors = FALSE)
    }

    # Expected patterns
    expected_rows <- LETTERS[1:16]  # A to P
    expected_cols <- as.character(1:24)  # 1 to 24

    # --- STEP 1: Find header with numbers 1-24 ---
    header_row <- NULL
    header_start_col <- NULL
    found_by_colnames <- FALSE

    # FIRST: Check if column names already have numbers 1-24
    current_colnames <- colnames(data)
    if (verbose) {
      message("Checking column names for pattern 1-24...")
    }

    # Look for sequence 1-24 in column names
    for (start_col in 1:(length(current_colnames) - 23)) {
      potential_header <- current_colnames[start_col:(start_col + 23)]

      # Check if they are numbers 1 to 24
      numeric_vals <- suppressWarnings(as.numeric(potential_header))
      valid_nums <- numeric_vals[!is.na(numeric_vals)]

      if (length(valid_nums) >= 20 && all(valid_nums >= 1 & valid_nums <= 24)) {
        header_row <- NA  # No separate header row
        header_start_col <- start_col
        found_by_colnames <- TRUE

        if (verbose) {
          message("  OK - Found pattern 1-24 in column names starting at column ", start_col)
          message("    Column names: ", paste(potential_header[1:6], collapse=", "),
                  ifelse(length(potential_header) > 6, "...", ""))
        }
        break
      }
    }

    # SECOND: If not found in column names, look in a row
    if (is.null(header_row)) {
      if (verbose) {
        message("Pattern not found in column names, searching in rows...")
      }

      for (i in 1:min(30, nrow(data))) {
        row_values <- as.character(data[i, ])

        # Look for sequence starting at each column
        for (start_col in 1:(length(row_values) - 23)) {
          potential_header <- row_values[start_col:(start_col + 23)]

          # Check if they are numbers 1 to 24
          numeric_vals <- suppressWarnings(as.numeric(potential_header))
          valid_nums <- numeric_vals[!is.na(numeric_vals)]

          if (length(valid_nums) >= 20 && all(valid_nums >= 1 & valid_nums <= 24)) {
            header_row <- i
            header_start_col <- start_col

            if (verbose) {
              message("  OK - Found header at row ", i, ", starting at column ", start_col)
              message("    Header values: ", paste(potential_header[1:6], collapse=", "),
                      ifelse(length(potential_header) > 6, "...", ""))
            }
            break
          }
        }
        if (!is.null(header_row)) break
      }
    }

    if (is.null(header_row)) {
      if (verbose) warning("  X Could not find pattern 1-24 in column names or rows")
      return(NULL)
    }

    # --- STEP 2: Find rows A-P ---
    data_start_row <- NULL
    data_end_row <- NULL

    # If we found by column names, start searching from row 1
    if (found_by_colnames) {
      search_start <- 1
    } else {
      search_start <- header_row + 1
    }

    # Look for rows A-P
    for (i in search_start:min(search_start + 50, nrow(data))) {
      # Check first column of the row (where letters should be)
      # Use header_start_col for the labels column
      if (header_start_col <= ncol(data)) {
        first_val <- as.character(data[i, header_start_col])

        # Check if it's one of the letters A-P
        if (first_val %in% expected_rows) {
          if (is.null(data_start_row)) {
            data_start_row <- i
            if (verbose) message("  OK - Found start of data at row ", i, " (row label: ", first_val, ")")
          }

          # Continue until something that's not A-P is found
          data_end_row <- i
        } else if (!is.null(data_start_row)) {
          # Already started finding data, but found something different
          break
        }
      }
    }

    if (is.null(data_start_row)) {
      if (verbose) warning("  X Could not find rows A-P")
      return(NULL)
    }

    # Check if we found all 16 rows
    found_rows_count <- data_end_row - data_start_row + 1
    if (found_rows_count != 16) {
      if (verbose) {
        warning("  X Found ", found_rows_count, " rows (expected 16)")
      }
    } else {
      if (verbose) message("  OK - Found all 16 rows (A-P)")
    }

    # Check which letters were found
    found_letters <- sapply(data_start_row:data_end_row, function(i) {
      as.character(data[i, header_start_col])
    })

    if (verbose) {
      message("  Found rows: ", paste(found_letters, collapse=", "))
    }

    # --- STEP 3: Validate data structure ---
    # Check if data is numeric
    if (header_start_col + 1 <= ncol(data)) {
      is_numeric <- FALSE
      check_end <- min(data_start_row + 2, data_end_row)
      for (.check_row in data_start_row:check_end) {
        sample_cell <- data[.check_row, header_start_col + 1]
        if (!is.na(suppressWarnings(as.numeric(as.character(sample_cell))))) {
          is_numeric <- TRUE
          break
        }
      }

      if (!is_numeric) {
        if (verbose) warning("  X Data doesn't appear to be numeric")
        return(NULL)
      }
    }

    # Define table boundaries
    result <- list(
      header_row = header_row,
      data_start_row = data_start_row,
      data_end_row = data_end_row,
      label_col = header_start_col,  # Column with letters A-P
      first_data_col = header_start_col + 1,  # First data column (after letters)
      last_data_col = header_start_col + 24,  # 24 data columns
      found_rows = found_letters,
      found_by_colnames = found_by_colnames,
      success = TRUE
    )

    if (verbose) {
      message("\n  Auto-detection successful:")
      message("    Found by: ", ifelse(found_by_colnames, "column names", "row header"))
      if (!found_by_colnames) message("    Header row: ", header_row)
      message("    Data rows: ", data_start_row, " to ", data_end_row, " (", found_rows_count, " rows)")
      message("    Data columns: ", result$first_data_col, " to ", result$last_data_col, " (24 columns)")
      message("    Row labels: ", paste(head(found_letters, 3), collapse=", "),
              ifelse(length(found_letters) > 3, "...", ""))
    }

    return(result)
  }

  # --- DATA VALIDATION ---
  if (nrow(data) < 16) {
    stop("Data must have at least 16 rows")
  }

  # --- DETECTION OR USE OF DEFAULT POSITIONS ---
  if (auto_detect) {
    if (verbose) message("\n=== AUTO-DETECTION MODE ===")

    positions <- detect_specific_pattern(data)

    if (is.null(positions)) {
      warning("Auto-detection failed. Using default positions.")
      # Use default positions
      header_row <- 9
      data_start_row <- 10
      data_end_row <- 25
      label_col <- 1
      first_data_col <- 2
      last_data_col <- 25
      found_by_colnames <- FALSE
    } else {
      header_row <- positions$header_row
      data_start_row <- positions$data_start_row
      data_end_row <- positions$data_end_row
      label_col <- positions$label_col
      first_data_col <- positions$first_data_col
      last_data_col <- positions$last_data_col
      found_letters <- positions$found_rows
      found_by_colnames <- positions$found_by_colnames

      if (verbose) message("OK - Using auto-detected table positions\n")
    }
  } else {
    # Use fixed positions (original behavior)
    header_row <- 9
    data_start_row <- 10
    data_end_row <- 25
    label_col <- 1
    first_data_col <- 2
    last_data_col <- 25
    found_by_colnames <- FALSE

    if (verbose) {
      message("Using fixed table positions (original behavior):")
      message("  Header row: ", header_row)
      message("  Data rows: ", data_start_row, " to ", data_end_row)
      message("  Data columns: ", first_data_col, " to ", last_data_col, "\n")
    }
  }

  # --- SET COLNAMES BASED ON DETECTION MODE ---
  if (found_by_colnames) {
    # If we found by column names, use existing column names
    # No need to redefine column names
    if (verbose) message("Using existing column names (pattern found in column names)")
  } else {
    # If we found by header row, set column names from that row
    colnames(data) <- as.character(data[header_row, ])
    if (verbose) message("Setting column names from row ", header_row)
  }

  # Extract the viability subtable
  viability_data <- data[data_start_row:data_end_row, label_col:last_data_col]

  # First column has row names (A-P)
  final_rownames <- as.character(viability_data[, 1])

  # Check if row names are A-P
  if (verbose) {
    message("Row names: ", paste(final_rownames, collapse=", "))

    if (all(final_rownames %in% LETTERS[1:16])) {
      message("OK - Valid row names (A-P)")
    } else {
      warning("X Row names don't match expected A-P pattern")
    }
  }

  # --- GET DATA COLUMN NAMES ---
  # Extract data column names (columns 2 onward)
  if (found_by_colnames) {
    # Column names are already correct
    data_colnames <- colnames(viability_data)[-1]  # Exclude label column
  } else {
    # Use header row values
    data_colnames <- as.character(data[header_row, first_data_col:last_data_col])
  }

  # Check if column names are numbers 1-24
  if (verbose) {
    message("Data column names: ", paste(data_colnames[1:6], collapse=", "),
            ifelse(length(data_colnames) > 6, "...", ""))

    numeric_colnames <- suppressWarnings(as.numeric(data_colnames))
    valid_numbers <- sum(!is.na(numeric_colnames) & numeric_colnames >= 1 & numeric_colnames <= 24)

    if (valid_numbers >= 20) {
      message("OK - Data columns have expected pattern (numbers 1-24)")
    } else {
      warning("X Data columns don't fully match expected 1-24 pattern")
    }
  }

  # --- COLUMN SELECTION LOGIC (ORIGINAL BEHAVIOR) ---
  if (!is.null(selected_columns)) {
    # Validate selected_columns are numeric indices
    if (!is.numeric(selected_columns)) {
      stop("selected_columns must be numeric indices (e.g., c(2:23))")
    }

    # ORIGINAL BEHAVIOR:
    # selected_columns refers to columns in the DATA TABLE (1-24)

    # Validate indices are between 1 and 24
    if (max(selected_columns) > 24) {
      stop("Selected column index ", max(selected_columns),
           " is out of bounds. Maximum allowed: 24")
    }
    if (min(selected_columns) < 1) {
      stop("Selected column indices must be >= 1")
    }

    # Validate number of columns is even (for split_replicates)
    if (split_replicates && length(selected_columns) %% 2 != 0) {
      warning("Number of selected data columns is not even (", length(selected_columns),
              " columns selected). This may cause issues with split_replicates.")
    }

    # Map to actual columns in viability_data table
    actual_columns_in_table <- selected_columns + 1  # +1 to skip label column

    # Check if columns exist
    if (max(actual_columns_in_table) > ncol(viability_data)) {
      stop("Selected columns exceed available columns in data table")
    }

    # Always keep column 1 (labels A-P) + selected columns
    columns_to_keep <- c(1, actual_columns_in_table)

    # Apply column selection
    viability_selected <- viability_data[, columns_to_keep, drop = FALSE]

    if (verbose) {
      message("Column selection applied (original behavior):")
      message("  Selected data columns: ", paste(selected_columns, collapse=", "))
      message("  Corresponding to table columns: ", paste(actual_columns_in_table, collapse=", "))
    }

  } else {
    # If no selection, use all columns
    viability_selected <- viability_data

    if (verbose) {
      message("No column selection - using all 24 data columns")
    }
  }

  # Convert to numeric data frame (all columns except first)
  viability_num <- as.data.frame(apply(viability_selected[, -1, drop = FALSE], 2, as.numeric))
  rownames(viability_num) <- final_rownames

  # Set column names correctly
  if (!is.null(selected_columns)) {
    # Use selected numbers as column names
    colnames(viability_num) <- as.character(selected_columns)
  } else {
    # Use data column names
    # If we have numeric names 1-24, use them, otherwise use 1-24
    if (length(data_colnames) == 24) {
      colnames(viability_num) <- data_colnames
    } else {
      colnames(viability_num) <- as.character(1:24)
    }
  }

  # --- CONTROL COLUMN MAPPING (ORIGINAL BEHAVIOR) ---
  map_control_column <- function(control_spec) {
    if (is.null(control_spec)) return(NULL)

    if (is.numeric(control_spec)) {
      # ORIGINAL BEHAVIOR: control_spec is index 1-24 referring to data

      # Validate it's in the correct range (1-24)
      if (control_spec < 1 || control_spec > 24) {
        stop("Control column index ", control_spec, " must be between 1 and 24")
      }

      # Check if this column is available
      col_name <- as.character(control_spec)

      if (!col_name %in% colnames(viability_num)) {
        # Try to find by original column name
        if (control_spec <= length(data_colnames)) {
          original_name <- data_colnames[control_spec]
          if (original_name %in% colnames(viability_num)) {
            col_name <- original_name
          } else {
            stop("Control column ", control_spec, " not available in selected columns")
          }
        } else {
          stop("Control column ", control_spec, " not available")
        }
      }

      if (verbose) {
        message("Control column mapping: index ", control_spec, " -> '", col_name, "'")
      }

      return(list(name = col_name,
                  user_index = control_spec,
                  is_relative = TRUE))
    } else {
      # User provided a column name directly
      # Try to find by name
      col_index <- which(colnames(viability_num) == control_spec)[1]

      if (is.na(col_index)) {
        stop("Control column '", control_spec, "' not found in data columns")
      }

      # Determine relative index (1-24)
      if (is.null(selected_columns)) {
        # If no selection, try to convert name to number
        relative_index <- suppressWarnings(as.numeric(control_spec))
        if (is.na(relative_index)) {
          # If not a number, look in data_colnames
          relative_index <- which(data_colnames == control_spec)[1]
        }
      } else {
        relative_index <- selected_columns[col_index]
      }

      return(list(name = control_spec,
                  user_index = NA,
                  relative_index = relative_index,
                  is_relative = FALSE))
    }
  }

  # Map controls
  control_0_info <- map_control_column(control_0perc)
  control_100_info <- map_control_column(control_100perc)

  # Check if controls are in available columns
  if (!is.null(control_0_info)) {
    if (!control_0_info$name %in% colnames(viability_num)) {
      stop("Control column '", control_0_info$name, "' (index ", control_0_info$user_index,
           ") not found in available columns.")
    }
  }

  if (!is.null(control_100_info)) {
    if (!control_100_info$name %in% colnames(viability_num)) {
      stop("Control column '", control_100_info$name, "' (index ", control_100_info$user_index,
           ") not found in available columns.")
    }
  }

  # --- LOW VALUE FILTERING ---
  replace_low_values <- function(df, threshold = low_value_threshold) {
    for (col in colnames(df)) {
      low_vals <- df[[col]] < threshold & !is.na(df[[col]])
      if (any(low_vals)) {
        n_replaced <- sum(low_vals, na.rm = TRUE)
        if (verbose) {
          warning("Replaced ", n_replaced, " value(s) < ", threshold, " with NA in column '", col, "'")
        }
        df[[col]][low_vals] <- NA
      }
    }
    return(df)
  }

  viability_num <- replace_low_values(viability_num)

  if (any(viability_num == 0, na.rm = TRUE)) {
    warning("Zero values detected in viability data - replacing with NA")
    viability_num[viability_num == 0] <- NA
  }

  # Create a modified copy for further processing
  viability_modified <- viability_num

  # --- INFO TABLE PROCESSING (for column renaming) ---
  if (!is.null(info_table)) {
    if (ncol(info_table) < 4) {
      stop("Info table must have at least 4 columns: log(inhibitor), Plate_Row, Construct, and Compound")
    }

    base_id_values <- paste(info_table[[3]], info_table[[4]], sep = ":")
    info_table$Base_ID <- base_id_values

    id_counts <- table(base_id_values)
    duplicate_ids <- names(id_counts)[id_counts > 1]

    # Distinguish biological replicates by adding suffixes
    if (length(duplicate_ids) > 0) {
      suffix_counter <- setNames(rep(1, length(duplicate_ids)), duplicate_ids)
      new_construct_values <- info_table[[3]]

      for (i in seq_along(base_id_values)) {
        current_base_id <- base_id_values[i]
        if (current_base_id %in% duplicate_ids) {
          if (suffix_counter[current_base_id] > 1) {
            new_construct_values[i] <- paste0(info_table[[3]][i], "_", suffix_counter[current_base_id])
          }
          suffix_counter[current_base_id] <- suffix_counter[current_base_id] + 1
        }
      }

      info_table$Construct_Modified <- new_construct_values

      if (verbose) {
        message("Found and automatically distinguished ", length(duplicate_ids),
                " biological replicate(s): ", paste(duplicate_ids, collapse = ", "))
      }
    } else {
      info_table$Construct_Modified <- info_table[[3]]
    }

    info_table$ID <- paste(info_table$Construct_Modified, info_table[[4]], sep = ":")
  }

  # --- APPLY CONTROL MEANS BY CONSTRUCT ---
  if (apply_control_means && !is.null(info_table) &&
      !is.null(control_0_info) && !is.null(control_100_info)) {

    # Extract relevant values from info_table
    plate_row_values <- info_table[[2]]
    construct_values <- info_table$Construct_Modified

    # Map plate rows to indices (A=1, B=2, etc.)
    plate_row_to_index <- setNames(1:16, LETTERS[1:16])

    # Group constructs by plate rows
    construct_groups <- list()
    unique_constructs <- unique(construct_values)

    for (construct in unique_constructs) {
      construct_indices <- which(construct_values == construct)
      construct_groups[[construct]] <- plate_row_values[construct_indices]
    }

    # Convert plate rows to numeric indices
    row_intervals <- lapply(construct_groups, function(plate_rows) {
      indices <- plate_row_to_index[plate_rows]
      as.numeric(indices[!is.na(indices)])
    })

    row_intervals <- row_intervals[sapply(row_intervals, length) > 0]

    # Control columns to analyze
    mean_columns <- c(control_0_info$name, control_100_info$name)
    existing_columns <- mean_columns[mean_columns %in% colnames(viability_modified)]

    if (length(existing_columns) > 0 && length(row_intervals) > 0) {
      if (verbose) {
        message("Applying construct-specific control means to modified viability table...")
      }

      # Apply means by construct to controls
      for (construct_name in names(row_intervals)) {
        interval_rows <- row_intervals[[construct_name]]
        valid_rows <- interval_rows[interval_rows >= 1 & interval_rows <= nrow(viability_modified)]

        if (length(valid_rows) > 0) {
          # Calculate means for this construct
          if (control_0_info$name %in% colnames(viability_modified)) {
            mean_background <- mean(viability_modified[valid_rows, control_0_info$name], na.rm = TRUE)
            # Replace individual values with mean
            viability_modified[valid_rows, control_0_info$name] <- mean_background
          }

          if (control_100_info$name %in% colnames(viability_modified)) {
            mean_positive <- mean(viability_modified[valid_rows, control_100_info$name], na.rm = TRUE)
            # Replace individual values with mean
            viability_modified[valid_rows, control_100_info$name] <- mean_positive
          }

          if (verbose) {
            message("  Construct '", construct_name, "' (rows ",
                    paste(LETTERS[range(valid_rows)], collapse = "-"), "):")
            if (control_0_info$name %in% colnames(viability_modified)) {
              message("    - Background control replaced with mean: ",
                      round(mean_background, 2))
            }
            if (control_100_info$name %in% colnames(viability_modified)) {
              message("    - Positive control replaced with mean: ",
                      round(mean_positive, 2))
            }
          }
        }
      }
    }
  }

  # --- COLUMN REORGANIZATION ---
  # Place control columns at beginning and end
  if (!is.null(control_0_info) && !is.null(control_100_info)) {
    control_cols <- c(control_0_info$name, control_100_info$name)

    missing_controls <- control_cols[!control_cols %in% colnames(viability_modified)]
    if (length(missing_controls) > 0) {
      stop("Control columns not found: ", paste(missing_controls, collapse = ", "))
    }

    other_columns <- setdiff(colnames(viability_modified), control_cols)
    viability_modified <- viability_modified[, c(control_0_info$name, other_columns, control_100_info$name)]
  }

  # Transpose the table
  viability_transposed <- as.data.frame(t(viability_modified))
  colnames(viability_transposed) <- rownames(viability_modified)

  # --- COLUMN RENAMING WITH INFO TABLE ---
  # Replace column names with IDs from info_table
  if (!is.null(info_table)) {
    mapping <- setNames(info_table$ID, info_table[[2]])
    new_colnames <- mapping[colnames(viability_transposed)]

    colnames(viability_transposed) <- ifelse(
      is.na(new_colnames),
      colnames(viability_transposed),
      new_colnames
    )
  }

  # --- TECHNICAL REPLICATE SPLITTING ---
  # Split experimental concentrations into two technical replicates
  final_table <- if (split_replicates) {
    split_replicates_func <- function(df) {
      n_rows <- nrow(df)
      if (n_rows < 3) {
        warning("Not enough rows to split replicates. Returning original table.")
        return(df)
      }

      # First and last rows are controls
      control_rows <- c(1, n_rows)
      exp_rows <- 2:(n_rows - 1)  # Middle rows are experimental concentrations
      split_pt <- floor(length(exp_rows) / 2)

      rep1_rows <- exp_rows[1:split_pt]
      rep2_rows <- exp_rows[(split_pt + 1):length(exp_rows)]

      new_table <- data.frame()

      for (col in colnames(df)) {
        rep1 <- df[c(control_rows[1], rep1_rows, control_rows[2]), col]
        rep2 <- df[c(control_rows[1], rep2_rows, control_rows[2]), col]

        if (ncol(new_table) == 0) {
          new_table <- data.frame(rep1, rep2)
          colnames(new_table) <- c(col, paste0(col, ".2"))
        } else {
          new_table[[col]] <- rep1
          new_table[[paste0(col, ".2")]] <- rep2
        }
      }

      new_rownames <- c(
        rownames(df)[control_rows[1]],
        rownames(df)[rep1_rows],
        rownames(df)[control_rows[2]]
      )

      rownames(new_table) <- new_rownames
      return(new_table)
    }
    split_replicates_func(viability_transposed)
  } else {
    viability_transposed
  }

  # --- ADD LOG(INHIBITOR) COLUMN ---
  if (!is.null(info_table)) {
    log_col <- info_table[[1]]
    n_needed <- nrow(final_table)

    if (length(log_col) > 0 && !is.na(log_col[1])) {
      if (verbose) {
        warning("First row of log(inhibitor) is not NA. Shifting values and setting first row to NA.")
      }

      adjusted_log_col <- c(NA, log_col[1:min(length(log_col), n_needed - 1)])

      if (length(adjusted_log_col) < n_needed) {
        adjusted_log_col <- c(adjusted_log_col, rep(NA, n_needed - length(adjusted_log_col)))
      }
    } else {
      adjusted_log_col <- if (length(log_col) >= n_needed) {
        log_col[1:n_needed]
      } else {
        c(log_col, rep(NA, n_needed - length(log_col)))
      }
    }

    final_table <- cbind(adjusted_log_col, final_table)
    colnames(final_table)[1] <- colnames(info_table)[1]
  }

  # --- FINAL RESULT ASSEMBLY ---
  result <- list()
  result$original_table <- viability_num
  result$modified_table <- final_table

  # Store additional information needed for quality control
  result$processing_info <- list(
    viability_data = viability_num,        # Raw viability signals
    viability_modified_with_means = viability_modified,  # Viability with means applied
    control_0_info = control_0_info,       # 0% control info
    control_100_info = control_100_info,   # 100% control info
    info_table = info_table,               # Info table with modifications
    final_rownames = final_rownames,       # Row names
    selected_columns = selected_columns,   # Selected columns info
    apply_control_means = apply_control_means,  # Control flag
    auto_detected = auto_detect,           # Whether auto-detection was used
    detection_method = if (auto_detect && exists("found_by_colnames")) {
      if (found_by_colnames) "column_names" else "row_header"
    } else "fixed_positions",
    table_positions = if (auto_detect) {
      list(header_row = header_row,
           data_start_row = data_start_row,
           data_end_row = data_end_row,
           label_col = label_col,
           first_data_col = first_data_col,
           last_data_col = last_data_col,
           found_by_colnames = found_by_colnames)
    } else {
      list(header_row = 9,
           data_start_row = 10,
           data_end_row = 25,
           label_col = 1,
           first_data_col = 2,
           last_data_col = 25,
           found_by_colnames = FALSE)
    }
  )

  # Add selected columns information
  result$selected_columns_info <- if (!is.null(selected_columns)) {
    list(
      user_indices = selected_columns,
      description = paste("Selected data columns: ", paste(selected_columns, collapse=", "),
                          " (referring to data columns 1-24)", sep=""),
      all_data_columns = 1:24,
      behavior = "original (indices 1-24 refer to data columns)"
    )
  } else {
    list(
      user_indices = 1:24,
      description = "All data columns (1-24)",
      all_data_columns = 1:24,
      behavior = "original (using all 24 data columns)"
    )
  }

  result$version <- 4.2
  result$data_type <- "viability"
  result$auto_detect <- auto_detect
  result$behavior_mode <- "original"

  return(result)
}
