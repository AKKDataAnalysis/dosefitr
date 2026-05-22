# Download the dosefitr package
# install.packages("remotes")       ## first install this package
# remotes::install_github("AKKDataAnalysis/dosefitr")

## NanoBRET analysis
library(dosefitr)

# Set your working directory to the folder containing your Excel data files
# before running this script, e.g.:
# setwd("C:/Users/YourName/Documents/my_experiment")

#### Create the Bret Ratio Table v1
results <- batch_ratio_analysis(
  control_0perc = 24,
  control_100perc = 12,
  output_dir = "./drc_quality",
  verbose = TRUE,
  low_value_threshold = 3000
)


#### Create the Bret Ratio Table v2
# results <- batch_ratio_analysis(
#   control_0perc = 16,
#   control_100perc = c(12,24),
#   output_dir = "./drc_quality",
#   verbose = TRUE,
#   low_value_threshold = 3000,
#   function_version = "v2"
# )


## If you want to use for Cell Viability
# via_results <- batch_viability_analysis(
#   control_0perc    = 13,
#   control_100perc  = 12,
#   selected_columns = c(2:23)   # passed straight through to process_viability_data
# )

## Detect and remove outliers
## (if you want to use the outliers table, put results_clean in the
## batch_results argument of batch_drc_analysis)
results_clean <- rout_outliers_batch(results, Q = 0.01)

outliers <- results_clean$outlier_summary
original_ratio_table <- results_clean$plate_01$result$modified_ratio_table_original  # choose which plate to inspect
outlier_table        <- results_clean$plate_01$result$modified_ratio_table            # choose which plate to inspect

## Plot the outlier curves
plot_outliers_batch_curves(results_clean)


## Merging plates for replicates (optional)
# merged <- merge_plate_replicates(results)


## Three-parameter logistic (3PL) dose-response model
drc_results <- batch_drc_analysis(
  batch_results = results_clean, ## or results (without outliers removed)
  normalize = TRUE,              ## normalization TRUE or FALSE
  output_dir = "./drc_results",
  verbose = TRUE
)

## Save all individual DRC plots (also generates a panel image per plate)
batch_save_all_drc_plots(
  batch_drc_results = drc_results,
  verbose = TRUE
)


## Plot multiple construct:compound curves together
plot_multiple_compounds(
  drc_results,
  compound_indices = 1:6,        # choose which compounds to include
  legend_text_size = 20,
  color_palette = "colorblind",
  plot_title = "",
  save_plot = TRUE,
  axis_text_size = 18,
  axis_title_size = 20,
  plate = "plate_01"             # select which plate to use
)


## Create SCARAB table
table <- scarab_table(
  results,
  drc_results,
  plate_name = "plate_01",       # choose plate
  date = "260303",
  experimenter_abbrev = "JD",    # abbreviation of the experimenter
  nLuc_orientation = "C",        # orientation of the nLuc
  tracer_kd_app = -7.5,          # value used for tracer_kd_app
  tracer_concentration_used = -7.5,
  tracer = "Tracer K10",
  decimal_separator = ","        # "," or "." for decimals
)


## Compare plates
compare_plates_drc(drc_results,
                   compare_by = "compound")
