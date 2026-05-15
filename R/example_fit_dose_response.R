#' @examples
#' # Create example dose-response data
#' example_data <- data.frame(
#'   log_conc = rep(c(-9, -8, -7, -6, -5, -4), 2),
#'   Compound_A_rep1 = c(10, 15, 45, 80, 95, 98, 12, 18, 48, 82, 94, 97),
#'   Compound_A_rep2 = c(12, 17, 43, 78, 96, 99, 11, 16, 46, 81, 95, 98),
#'   Compound_B_rep1 = c(85, 70, 45, 25, 15, 10, 87, 72, 43, 23, 14, 9),
#'   Compound_B_rep2 = c(83, 68, 47, 27, 16, 11, 86, 71, 44, 24, 15, 10)
#' )
#'
#' # Analyze the data
#' results <- fit_dose_response(example_data, normalize = TRUE)
#'
#' # View summary results
#' print(results$summary_table)
#'
#' # Check success rate
#' cat("Success rate:", results$successful_fits / results$n_compounds * 100, "%\n")
#'
#' # Export results to Excel
#' \dontrun{
#' fit_dose_response(example_data, output_file = "dose_response_results.xlsx")
#' }