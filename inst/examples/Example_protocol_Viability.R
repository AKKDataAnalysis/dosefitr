# Download the dosefitr package
# install.packages("remotes")       ## first install this package
# remotes::install_github("AKKDataAnalysis/dosefitr")

## NanoBRET analysis
library(dosefitr)

# Set your working directory to the folder containing your Excel data files
# before running this script, e.g.:
# setwd("C:/Users/YourName/Documents/my_experiment")


## If you want to use for Cell Viability
via_results <- batch_viability_analysis(
 control_0perc    = 13,
 control_100perc  = 12,
 selected_columns = c(2:23)   # passed straight through to process_viability_data
)

## Detect and remove outliers
## (if you want to use the outliers table, put results_clean in the
## batch_results argument of batch_drc_analysis)
results_clean <- rout_outliers_batch(via_results, Q = 0.01, keep_cytotoxic = TRUE)

outliers <- results_clean$outlier_summary
original_ratio_table <- results_clean$plate_01$result$modified_ratio_table_original  # choose which plate to inspect
outlier_table        <- results_clean$plate_01$result$modified_ratio_table            # choose which plate to inspect

## Plot the outlier curves
plot_outliers_batch_curves(results_clean)


## Merging plates for replicates (optional)
# merged <- merge_plate_replicates(via_results)


## Three-parameter logistic (3PL) dose-response model
drc_results <- batch_drc_analysis(
  batch_results = results_clean, ## or via_results (without outliers removed)
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



# Viability Scarab
table <- scarab_viability(
 via_results, 
 drc_results,
 cell_line = "AGP-01",
 cell_type = "Metastatic Gastric Adenocarcinoma",
 date = "260521",
 plate_name = "plate_02",
 experimenter_abbrev = "TL",
 decimal_separator = ","
)


## Compare plates
compare_plates_drc(drc_results,
                   compare_by = "compound")