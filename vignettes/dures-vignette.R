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
folder_path <- "/data1/sbanerjee/test_1/"  # Update this to your actual working path to the test data 

# Files and folders to keep
keep_files <- c("feature_list.txt", "Sample-and-feature-wise-RT-tolerance.txt", "mzml_files")

# List all contents in folder_path
all_contents <- list.files(folder_path, full.names = TRUE)

# Delete everything except essentials
to_delete <- all_contents[!basename(all_contents) %in% keep_files]
unlink(to_delete, recursive = TRUE, force = TRUE)

message("✅ Cleaned folder_path — kept only required files/folders.")

## ----initial-setup, eval=FALSE------------------------------------------------
# # Dependencies
# install.packages("BiocManager", quiet = TRUE)
# BiocManager::install(c("Spectra", "S4Vectors", "mzR"), quiet = TRUE)
# 
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools", quiet = TRUE)
# }

## ----initial-setup-contd, eval=TRUE-------------------------------------------
devtools::install_github("BiosystemEngineeringLab-IITB/dures")
library(dures)
library(Spectra)
library(mzR)

## ----execute first step of the workflow---------------------------------------
# Define the path to the test_1 folder (replace with the actual path on your system)
folder_path <- "/data1/sbanerjee/test_1//"

# Check if the directory exists
if (!dir.exists(folder_path)) {
  stop("The directory test_1 does not exist. Please modify the path to the actual location of the test_1 folder")
}

l1 <- preprocess(folder_path = folder_path, tol_mz = 5, tol_rt = 0.1)

## ----first element------------------------------------------------------------
cat("name of first feature\n")
names(l1$spectra_ms2_only)[1] #name of the first feature
l1$spectra_ms2_only[[1]]

## ----second-------------------------------------------------------------------
cat("name of third feature\n")
names(l1$spectra_ms2_only)[3] #name of the first feature
l1$spectra_ms2_only[[3]]

## ----peaksdata----------------------------------------------------------------
peaksData(l1$spectra_ms2_only[[1]]$aggregate)[[1]] #first replicate spectrum
peaksData(l1$spectra_ms2_only[[2]]$aggregate)[[1]] #second replicate spectrum

## ----spectrum identity--------------------------------------------------------
l1$spectra_ms2_only[[1]]$aggregate$spectrumId[1:10] #will be a vector of length 83


## ----second step--------------------------------------------------------------
l2 = extract_raw_spectra(folder_path = folder_path, l1, 0.05, 0.8)

## ----check_and_assign, echo=TRUE, error=TRUE----------------------------------
try({
# Check if 'sps_top80_tic_2' exists before assigning
if (!exists("sps_top80_tic_2", envir = .dures_env)) {
  stop("Object 'sps_top80_tic_2' not found in .dures_env.")
} else {
  assign("sps_top80_tic_2", get("sps_top80_tic_2", envir = .dures_env), envir = .dures_env)
}

})

## ----intra--------------------------------------------------------------------
l2$df[c(1,3),]
print(l2$sps_top_tic_2[[1]]) # From 83 before to 66 spectra after after top x% TIC
print(l2$sps_top_tic_2[[3]]) # From 43 before to 34 spectra remain after top x% TIC

## ----compare_spectra----------------------------------------------------------
print(l1$spectra_ms2_only[[1]]$aggregate$spectrumId[1:10]) #spectra identities before step 2
print(l2$sps_top_tic_2[[1]]$spectrumId[1:10]) #spectra identitied after x% TIC cutoff implementation in step 2

## ----peaksdata_before_intra---------------------------------------------------
peaksData(l2$sps_top_tic_2[[1]][3])[[1]]

## ----after_grouping-----------------------------------------------------------
spectrum_size_after_grouping = read.delim("/data1/sbanerjee/test_1/MS2_scans_before_denoising/1982/20092020_CovidMild_P294_P_2.mzML_scan_1587.txt")
print(dim(spectrum_size_after_grouping)[1])

## ----another example----------------------------------------------------------
print(head(peaksData(l2$sps_top_tic_2[[3]][1])[[1]]))
print(dim(peaksData(l2$sps_top_tic_2[[3]][1])[[1]])[1])

## ----after_grouping_1---------------------------------------------------------
spectrum_size_after_grouping = read.delim("/data1/sbanerjee/test_1//MS2_scans_before_denoising/872/13092020_CovidMild_P282_P_2.mzML_scan_1120.txt")
print(dim(spectrum_size_after_grouping)[1])

## ----third step---------------------------------------------------------------
l3 = call_aggregate(l2$sps_top_tic_2, 0.05, folder_path)

## ----inspect-l3---------------------------------------------------------------
print(head(l3[[1]]$Df))
print(dim(l3[[1]]$Df)[1])
l3[[1]]$Mean #will return one consensus spectrum


## ----recur--------------------------------------------------------------------
l = l3[[1]]$Df
print(head(l[order(l$Frequency, decreasing = TRUE),],10))

## ----label_spec---------------------------------------------------------------
l4 = label_individual_spectrum(l3, folder_path, 0.05)

## ----check step4--------------------------------------------------------------
print(length(l4[[1]]))
print(names(l4[[1]][3]))
print(l4[[1]][3])

## ----last_step----------------------------------------------------------------
l5 = generate_denoised_spectra(l4, folder_path, 0.12, ion_mode = "pos") 

## ----eg-----------------------------------------------------------------------
#20092020_CovidMild_P294_P_2.mzML_scan_1587
backend = Spectra::MsBackendMzR()
sps <- Spectra::Spectra("/data1/sbanerjee/test_1/Denoised_spectra_mzML/20092020_CovidMild_P294_P_2.mzML", source = backend)
print(peaksData(l2$sps_top_tic_2[[1]][3])[[1]])

## ----test 5-------------------------------------------------------------------
print(l4[[1]][3])


## ----test 6-------------------------------------------------------------------
print(peaksData(sps)[[1]])

