#' @examples
#' \dontrun{
#' # Example 1: Complete workflow with quality assessment
#' # Load raw data and sample information
#' raw_data <- read.csv("dose_response_raw.csv")
#' sample_info <- read.csv("sample_information.csv")
#'
#' # Process with comprehensive quality metrics
#' analysis_result <- ratio_dose_response(
#'   data = raw_data,
#'   control_0perc = "DMSO_Control",
#'   control_100perc = "Staurosporine_Control", 
#'   info_table = sample_info,
#'   save_to_excel = "complete_analysis.xlsx"
#' )
#'
#' # Review quality metrics
#' quality_table <- analysis_result$interval_means
#' print(quality_table)
#'
#' # Identify high-quality assays
#' high_quality <- quality_table[, quality_table["Overall_Quality", ] == "high"]
#' medium_quality <- quality_table[, quality_table["Overall_Quality", ] == "medium"]
#'
#' cat("High quality assays:", ncol(high_quality), "\n")
#' cat("Medium quality assays:", ncol(medium_quality), "\n")
#'
#' # Example 2: Handling biological replicates
#' # Sample info with replicates (same Target-Compound combination)
#' sample_info_with_replicates <- data.frame(
#'   log_inhibitor = c(-9, -8, -7, -6, -5, -4),
#'   Plate_Row = c("A", "B", "C", "D", "E", "F"),
#'   Target = c("EGFR", "EGFR", "BRAF", "BRAF", "PI3K", "PI3K"),
#'   Compound = c("Gefitinib", "Gefitinib", "Dabrafenib", "Dabrafenib", "LY294002", "LY294002")
#' )
#'
#' result <- ratio_dose_response(
#'   data = raw_data,
#'   control_0perc = "Control_0",
#'   control_100perc = "Control_100", 
#'   info_table = sample_info_with_replicates
#' )
#'
#' # Check how replicates were handled
#' print(colnames(result$modified_ratio_table))
#' # Output: "EGFR_2-Gefitinib", "BRAF_2-Dabrafenib", etc.
#'
#'
#' # Example 5: Batch processing multiple plates
#' process_plate <- function(plate_data, plate_info, output_suffix) {
#'   result <- ratio_dose_response(
#'     data = plate_data,
#'     control_0perc = "DMSO_Ctrl",
#'     control_100perc = "Positive_Ctrl",
#'     info_table = plate_info,
#'     save_to_excel = paste0("plate_", output_suffix, ".xlsx")
#'   )
#'   return(result)
#' }
#'
#' # Process multiple plates
#' plate1_result <- process_plate(plate1_data, plate1_info, "1")
#' plate2_result <- process_plate(plate2_data, plate2_info, "2")
#' plate3_result <- process_plate(plate3_data, plate3_info, "3")
#'
#' # Combine quality metrics
#' all_quality <- cbind(
#'   Plate1 = plate1_result$interval_means,
#'   Plate2 = plate2_result$interval_means, 
#'   Plate3 = plate3_result$interval_means
#' )
#'
#' write.csv(all_quality, "combined_quality_metrics.csv")
#' }
#' 
#' 
#' 
#' @section Raw Data Preparation Guide:
#' 
#' \strong{Before using this function, ensure your raw data is properly prepared:}
#' 
#' \strong{1. Data Export from Instruments:}
#' \itemize{
#'   \item Export data in CSV or Excel format with minimal formatting
#'   \item Preserve original row structure and column headers
#'   \item Include all measurement data without manual editing
#' }
#' 
#' \strong{2. Expected Raw Data Format:}
#' \itemize{
#'   \item \strong{Rows 1-8}: May contain instrument metadata, plate maps, or experiment notes
#'   \item \strong{Row 9}: Must contain column identifiers (well positions or sample names)
#'   \item \strong{Rows 10-25}: First measurement set (typically 16 rows A-P)
#'   \item \strong{Rows 28-43}: Second measurement set (typically 16 rows A-P)
#'   \item \strong{Other rows}: Ignored but preserved in structure
#' }
#' 
#' \strong{3. Control Well Identification:}
#' \itemize{
#'   \item Control columns must be explicitly named in row 9
#'   \item Use consistent naming across experiments (e.g., "DMSO_Ctrl", "Max_Inhibition")
#'   \item Controls should be present in both measurement subtables
#' }
#' 
#' \strong{4. Common Raw Data Issues:}
#' \itemize{
#'   \item \strong{Mixed data types}: Ensure all data cells contain numeric values only
#'   \item \strong{Missing values}: Use NA or blank cells for failed measurements
#'   \item \strong{Format changes}: Avoid manual formatting that alters row structure
#'   \item \strong{Header modifications}: Do not modify row 9 column names
#' }
#' 
#' \strong{5. Example Raw Data Structure:}
#' \preformatted{
#'        [,1]       [,2]          [,3]          [,4]          [,5]
#' [1,]  "Experiment: DRG_001"    ""            ""            ""           ""
#' [2,]  "Plate: 1"              ""            ""            ""           ""
#' [3,]  "Date: 2024-01-15"      ""            ""            ""           ""
#' [4,]  ""                      ""            ""            ""           ""
#' [5,]  "Measurement 1: Donor"  ""            ""            ""           ""
#' [6,]  ""                      ""            ""            ""           ""
#' [7,]  "Temperature: 37C"      ""            ""            ""           ""
#' [8,]  ""                      ""            ""            ""           ""
#' [9,]  ""           "DMSO_Ctrl"  "Compound_1"  "Compound_2"  "Stauro_Ctrl"
#' [10,] "A"          15020       16240         15890         15210
#' [11,] "B"          14980       15560         15230         14890
#' ...    ...         ...         ...           ...           ...
#' [28,] "Measurement 2: Acceptor" ""           ""            ""           ""
#' [29,] "A"          12015       11580         11240         9850
#' [30,] "B"          11890       11620         11380         9920
#' ...    ...         ...         ...           ...           ...
#' }
#' 
#' @section Raw Data Quality Checks:
#' The function automatically performs these quality checks on raw data:
#' \itemize{
#'   \item \strong{Signal intensity}: Flags measurements <1000 as potentially failed
#'   \item \strong{Data consistency}: Checks for numeric values in measurement areas
#'   \item \strong{Control validity}: Verifies control columns exist and contain data
#'   \item \strong{Structure integrity}: Confirms required row ranges contain data
#' }