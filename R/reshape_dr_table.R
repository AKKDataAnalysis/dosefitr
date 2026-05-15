#' Reshape Dose-Response Results Table to Standardized Format
#'
#' Transforms dose-response analysis results from wide to structured format,
#' organizing parameters into logical sections following industry standards.
#'
#' @param results_table Data frame containing dose-response analysis results
#'   from \code{\link{fit_drc_3pl}} function.
#' @param output_file Character string specifying the output Excel file path
#'   (default: NULL, no file saved).
#' @param decimal_comma Logical indicating whether to use comma as decimal separator
#'   instead of point (default: FALSE).
#'
#' @return A data frame with the following structure:
#' \itemize{
#'   \item Rows organized into logical sections (best-fit values, confidence intervals, etc.)
#'   \item Compounds as columns
#'   \item Parameters as rows with section headers
#'   \item Professional formatting suitable for publications
#' }
#'
#' @details
#' This function restructures dose-response analysis results into a standardized
#' format that follows industry conventions for pharmacological data presentation.
#' The output organizes parameters into logical sections that facilitate
#' interpretation and comparison across multiple compounds.
#'
#' \strong{Output Structure:}
#' The table is organized into the following sections:
#' \itemize{
#'   \item \strong{log(inhibitor) vs. response (three parameters)}: Section header
#'   \item \strong{Best-fit values}: Bottom, Top, LogIC50, IC50, Span
#'   \item \strong{Lower 95% conf. limit}: Profile likelihood confidence intervals
#'   \item \strong{Upper 95% conf. limit}: Profile likelihood confidence intervals  
#'   \item \strong{Goodness of Fit}: Degrees of Freedom, R-squared, Sum of Squares, Syx
#'   \item \strong{Additional metrics}: Max Slope, Curve Quality
#' }
#'
#' \strong{Input Requirements:}
#' The input table must contain these required columns:
#' \code{Compound, Bottom, Top, LogIC50, IC50, Span, R_squared, Syx, Max_Slope, 
#' Curve_Quality, Degrees_of_Freedom}
#'
#' @examples
#' \dontrun{
#' # Perform dose-response analysis
#' analysis_results <- fit_dose_response(my_data, normalize = TRUE)
#'
#' # Reshape the summary table
#' reshaped_table <- reshape_dr_table(analysis_results$summary_table)
#'
#' # View the structured results
#' print(reshaped_table)
#'
#' # Save to Excel with European decimal format
#' reshape_dr_table(
#'   analysis_results$summary_table,
#'   output_file = "dose_response_structured.xlsx",
#'   decimal_comma = TRUE
#' )
#'
#' # Use with custom results
#' custom_results <- data.frame(
#'   Compound = c("Drug_A", "Drug_B", "Drug_C"),
#'   Bottom = c(10.2, 15.6, 8.9),
#'   Top = c(95.8, 92.3, 97.1),
#'   LogIC50 = c(-7.23, -6.89, -7.45),
#'   IC50 = c(5.89e-08, 1.29e-07, 3.55e-08),
#'   Span = c(85.6, 76.7, 88.2),
#'   R_squared = c(0.987, 0.956, 0.992),
#'   Syx = c(2.34, 3.45, 1.89),
#'   Max_Slope = c(-18.5, -16.2, -19.8),
#'   Curve_Quality = c("Good curve", "Good curve", "Good curve"),
#'   Degrees_of_Freedom = c(8, 8, 8)
#' )
#'
#' structured_table <- reshape_dr_table(custom_results)
#' print(structured_table)
#' }
#'
#' @section Output Format:
#' The reshaped table has the following characteristics:
#' \itemize{
#'   \item \strong{Columns}: Compound names
#'   \item \strong{Rows}: Parameters organized in logical groups
#'   \item \strong{Section headers}: Descriptive headers for parameter groups
#'   \item \strong{Professional ordering}: Follows standard pharmacological reporting
#'   \item \strong{Excel-ready}: Directly exportable to formatted spreadsheets
#' }
#'
#' @section Decimal Formatting:
#' When \code{decimal_comma = TRUE}:
#' \itemize{
#'   \item Converts decimal points to commas (e.g., 1.234 -> 1,234)
#'   \item Suitable for European and many international formats
#'   \item Only affects Excel export, not the returned data frame
#' }
#'
#' @seealso
#' \code{\link{fit_drc_3pl}} for generating analysis results
#' \code{\link{save_multiple_sheets}} for advanced Excel export options
#' \code{\link{plot_multiple_compounds}} for visual comparison of results
#'
#' @export
#'
#' @references
#' For pharmacological data reporting standards:
#' \itemize{
#'   \item "The IUPHAR/BPS Guide to PHARMACOLOGY" data standards
#'   \item Journal of Pharmacology and Experimental Therapeutics guidelines
#'   \item Nature Scientific Data reporting standards
#' }




reshape_dr_table <- function(results_table, output_file = NULL, decimal_comma = FALSE) {
  
  # Check if input is final_summary_table (transposed) or summary_table (original)
  is_final_summary <- !"Compound" %in% colnames(results_table) && 
    all(rownames(results_table) %in% c("Bottom", "Top", "LogIC50", "IC50", "Span", 
                                       "Bottom_Lower_95CI", "Bottom_Upper_95CI",
                                       "Top_Lower_95CI", "Top_Upper_95CI",
                                       "LogIC50_Lower_95CI", "LogIC50_Upper_95CI",
                                       "IC50_Lower_95CI", "IC50_Upper_95CI",
                                       "R_squared", "Syx", "Sum_of_Squares", 
                                       "Degrees_of_Freedom", "Max_Slope", 
                                       "Ideal_Hill_Slope", "Curve_Quality"))
  
  if (is_final_summary) {
    # Working with final_summary_table (already transposed)
    transposed <- results_table
    compound_names <- colnames(transposed)
    
  } else {
    # Working with summary_table (original format) - maintain original logic
    required_cols <- c("Compound", "Bottom", "Top", "LogIC50", "IC50", "Span", 
                       "R_squared", "Syx", "Max_Slope", "Curve_Quality", "Degrees_of_Freedom")
    
    if (!all(required_cols %in% colnames(results_table))) {
      stop("Table does not have the expected dose-response analysis structure")
    }
    
    # Transpose table - compounds become columns
    transposed <- as.data.frame(t(results_table[, -1]))
    colnames(transposed) <- results_table$Compound
  }
  
  # Define all sections and parameters in order
  sections <- list(
    list(header = "log(inhibitor) vs. response (three parameters)", params = NULL),
    list(header = "Best-fit values", params = c("Bottom", "Top", "LogIC50", "IC50", "Span")),
    list(header = "Lower 95% conf. limit (profile likelihood)", 
         params = c("Bottom_Lower_95CI", "Top_Lower_95CI", "LogIC50_Lower_95CI", "IC50_Lower_95CI")),
    list(header = "Upper 95% conf. limit (profile likelihood)", 
         params = c("Bottom_Upper_95CI", "Top_Upper_95CI", "LogIC50_Upper_95CI", "IC50_Upper_95CI")),
    list(header = "Goodness of Fit", 
         params = c("Degrees_of_Freedom", "R_squared", "Sum_of_Squares", "Syx")),
    list(header = "Additional Parameters", 
         params = c("Max_Slope", "Ideal_Hill_Slope", "Curve_Quality"))
  )
  
  # Build final table structure
  final_table <- data.frame(matrix(NA, nrow = 0, ncol = ncol(transposed)))
  colnames(final_table) <- colnames(transposed)
  
  # Populate table with sections and parameters
  for (section in sections) {
    # Add section header if it exists
    if (!is.null(section$header)) {
      header_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(transposed)))
      colnames(header_row) <- colnames(transposed)
      rownames(header_row) <- section$header
      final_table <- rbind(final_table, header_row)
    }
    
    # Add section parameters
    for (param in section$params) {
      if (param %in% rownames(transposed)) {
        row_data <- transposed[param, , drop = FALSE]
        final_table <- rbind(final_table, row_data)
      }
    }
  }
  
  # Save to Excel if output_file specified
  if (!is.null(output_file)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("The 'openxlsx' package is required to save Excel files. Install using: install.packages('openxlsx')")
    }
    
    # Prepare table for Excel export
    excel_table <- final_table
    excel_table$Parameter <- rownames(excel_table)
    excel_table <- excel_table[, c("Parameter", colnames(final_table))]
    
    # Convert decimal points to commas if requested
    if (decimal_comma) {
      for (col in 2:ncol(excel_table)) {
        excel_table[[col]] <- as.character(excel_table[[col]])
        excel_table[[col]] <- gsub("\\.", ",", excel_table[[col]])
        excel_table[[col]][is.na(excel_table[[col]])] <- NA
      }
    }
    
    # Save to Excel
    openxlsx::write.xlsx(excel_table, output_file, rowNames = FALSE)
    cat("Table saved to:", output_file, "\n")
  }
  
  return(final_table)
}
