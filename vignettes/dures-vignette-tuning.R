## ----setup, include = FALSE---------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org/"))
#this code block will
# Clean your folder_path every time the vignette is run
# 
# Prevent file conflicts and stale output errors
# 
# Ensure reproducibility for every re-render
# Suppress title check if necessary
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# === Custom folder cleanup ===

# Define the base folder path
folder_path <- "/data1/sbanerjee/test_data_MTBLS606/"  # Update this to your actual working path to the test data 

# Files and folders to keep
keep_files <- c("feature_list.txt", "Sample-and-feature-wise-RT-tolerance.txt", "mzml_files")

# List all contents in folder_path
all_contents <- list.files(folder_path, full.names = TRUE)

# Delete everything except essentials
to_delete <- all_contents[!basename(all_contents) %in% keep_files]
unlink(to_delete, recursive = TRUE, force = TRUE)

message("✅ Cleaned folder_path — kept only required files/folders.")

## ----initial-setup, eval=FALSE------------------------------------------------
# install_and_load_dures_dependencies()

## ----initial-setup-contd, eval=TRUE-------------------------------------------
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", quiet = TRUE)
}
devtools::install_github("BiosystemEngineeringLab-IITB/dures", auth_token = NULL)
library(dures)

## ----execute first step of the workflow---------------------------------------
# Define the path to the  folder (replace with the actual path on your system)
folder_path <- "/data1/sbanerjee/test_data_MTBLS606/"

# Check if the directory exists
if (!dir.exists(folder_path)) {
  stop("The directory test_1 does not exist. Please modify the path to the actual location of the test_1 folder")
}

l1 <- preprocess(folder_path = folder_path, tol_mz = 25, tol_rt = 0.1)

## ----execute intermediate step of workflow------------------------------------

print(paste0("Length of stats file after removing features from which MS2 spectra couldn't be extracted: ", as.character(dim(l1$stats_file_ms2_only)[1])))

## ----Subsetting the feature list to reduce computation------------------------
# Set seed for reproducibility
set.seed(123)

# Step 1: Get all IDs from the dataframe
all_ids <- l1$stats_file_ms2_only$ID

# Step 2: Sample 50 IDs
selected_ids <- sample(all_ids, 2)

# Step 3: Subset dataframe and reorder by selected_ids
subset_stats <- l1$stats_file_ms2_only[l1$stats_file_ms2_only$ID %in% selected_ids, ]
subset_stats <- subset_stats[match(selected_ids, subset_stats$ID), ]

# Step 4: Subset and reorder spectra list to match the same order
subset_spectra <- l1$spectra_ms2_only[selected_ids]  # list is named by ID

# Step 5: Store both in a new list
l1_subset <- list(
  stats_file_ms2_only = subset_stats,
  spectra_ms2_only = subset_spectra
)

## ----second step--------------------------------------------------------------
l2 = extract_raw_spectra(folder_path = folder_path, l1_subset, 0.05, 0.8)

## ----check_and_assign, echo=TRUE, error=TRUE----------------------------------
try({
# Check if 'sps_top80_tic_2' exists before assigning
if (!exists("sps_top80_tic_2", envir = .dures_env)) {
  stop("Object 'sps_top80_tic_2' not found in .dures_env.")
} else {
  assign("sps_top80_tic_2", get("sps_top80_tic_2", envir = .dures_env), envir = .dures_env)
}

})

## ----third step - most expensive----------------------------------------------
l3 = call_aggregate(l2$sps_top_tic_2, 0.05, folder_path)

## ----label_spec---------------------------------------------------------------
l4 = label_individual_spectrum(l3, folder_path, 0.05)

