#' precursor_matching
#'
#' Matches precursor m/z with library within a given tolerance (As outlined in tol_mz param in preprocess()) and generates a library of reference spectra that can be used for MS/MS-level fragment matching. Also groups fragments in the reference spectra that are close-by within a given tolerance tol
#'
#' @param l1 output from the preprocess()
#' @param folder_path folder containing the required input directory (mzml) and feature list file in .txt format.
#' @param  ionization_mode positive or negative. Based on this input, the library data is loaded and stored as a global variable
#' @param tol tolerance to group reference spectra fragments, defaults to 0.01
#' @return A dataframe with only those precursors, their mz, RT, etc, which had a corresponding match with the reference library
#' @examples
#' # Example usage of the function
#' precursor_matching(l1, folder_path, "positive", 0.05)
#' @export



precursor_matching <- function(l1, folder_path, ionization_mode, tol = 0.01) {
  path <- file.path(folder_path, "Reference/")

  # Create output directory if it doesn't exist
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  # Load the library only once and store it globally
  if (!exists(".positive_lib", envir = .GlobalEnv) && ionization_mode == "positive") {
    file_path_pos_lib <- system.file("extdata", "library_positive.rds", package = "dures")
    .GlobalEnv$.positive_lib <- readRDS(file_path_pos_lib)
    cat("Loaded positive mode library from file\n")
  }

  if (!exists(".negative_lib", envir = .GlobalEnv) && ionization_mode == "negative") {
    file_path_neg_lib <- system.file("extdata", "library_negative.rds", package = "dures")
    .GlobalEnv$.negative_lib <- readRDS(file_path_neg_lib)
    cat("Loaded negative mode library from file\n")
  }

  # Use the preloaded library
  lib <- if (ionization_mode == "positive") .GlobalEnv$.positive_lib else .GlobalEnv$.negative_lib
  cat("Matching precursor with library in", ionization_mode, "mode\n")

  # Process each entry
  for (j in seq_len(nrow(l1$stats_file_ms2_only))) {
    message("Processing metabolite: ", l1$stats_file_ms2_only$ID[j])
    name <- l1$stats_file_ms2_only$ID[j]
    p <- file.path(path, name)
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)

    # Prepare references for matching
    sps_ref <- reference_preparation_for_matching(j, l1$stats_file_ms2_only, lib, p)
  }

  ### **Step 1: Identify and Remove Empty Folders**
  ### **Step 1: Identify and Remove Empty Folders**
  all_folders <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  empty_folders <- all_folders[sapply(all_folders, function(f) length(list.files(f, recursive = TRUE)) == 0)]

  # Remove empty folders
  if (length(empty_folders) > 0) {
    message("Removing empty folders: ", paste(basename(empty_folders), collapse = ", "))
    unlink(empty_folders, recursive = TRUE, force = TRUE)  # Delete empty directories
  }

  empty_folder_names <- basename(empty_folders)  # Store removed folder names


  ### **Step 2: Update `l1$stats_file_ms2_only`**
  new_l1_stats <- l1$stats_file_ms2_only[!l1$stats_file_ms2_only$ID %in% empty_folder_names, ]

  ### **Step 3: Update `sps_top80_tic_3`**
  remaining_names <- setdiff(list.files(path), empty_folder_names)  # Names of non-empty folders
  sps_top80_tic_2 <- get("sps_top80_tic_2", envir = .dures_env)

  # Find metabolites that were removed (those in sps_top80_tic_2 but not in remaining_names)
  removed_metabolites <- empty_folder_names

  sps_top80_tic_3 <- sps_top80_tic_2[names(sps_top80_tic_2) %in% remaining_names]

  if (length(removed_metabolites) > 0) {
    message("Precursor matching with the library did not yield spectra for the following metabolites: ",
            paste(removed_metabolites, collapse = ", "))
  }


  ### **Step 4: Assign `sps_top80_tic_3` to `.dures_env`**
  assign("sps_top80_tic_3", sps_top80_tic_3, envir = .dures_env)

  ### **Step 5: Perform Intrascan Grouping and Overwrite Files in `Reference/`**
  for (i in seq_len(nrow(new_l1_stats))) {

    metabolite <- new_l1_stats$ID[i]
    metabolite_folder <- file.path(path, metabolite)

    files <- list.files(metabolite_folder, full.names = TRUE)
    if (length(files) == 0) next  # Skip empty directories

    for (file in files) {
      scan <- read.delim(file)
      scan_name <- basename(file)

      # Sorting m/z values in decreasing order
      scan <- scan[order(scan$fragments, decreasing = TRUE), ]
      scan$unique_id <- paste("frag", scan$fragments, sep = "_")

      # Perform intrascan grouping
      scan <- scan %>%
        dplyr::arrange(fragments) %>%
        mutate(group = cumsum(c(1, diff(fragments) > tol)))

      agg_df <- scan %>%
        dplyr::group_by(group) %>%
        dplyr::summarize(fragments = mean(fragments), intensity = sum(intensity), .groups = "drop")

      # Save new grouped scan data (overwrite the original file)
      new_scan <- data.frame(fragments = agg_df$fragments, intensity = agg_df$intensity)
      #print(paste("Before grouping:", nrow(scan)))
      #print(paste("After grouping:", nrow(new_scan)))

      write.table(new_scan, file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }


  ### **Step 6: Return Updated Dataframe**
  return(new_l1_stats)
}










