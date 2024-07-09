#' extract_raw_spectra
#'
#' Extracts MS2 spectra for every feature, groups fragments within a given mass tolerance and store only the top x% TIC spectra for further analysis. It is required to run preprocess() before running this.
#'
#' @param folder_path: same folder path used in preprocess(). The extracted fragment-grouped spectra in .txt format will be stored here
#' @param list_of_spectra_files: output from preprocess()
#' @param intrascan_grouping_tolerance: if fragments are found within this mass tolerance (in Da), they will be grouped into a single fragment by taking average mz value and sum of intensities. Default value is 0.05 Da.
#' @param top_percent_tic: this parameter will tell the user that only top x% TIC spectra will be extracted and stored in the folder in .txt format. Default value is 0.8 (or 80%)
#' @return a list containing all top x% TIC scans in spectra format, a dataframe showing the number of spectra remaining before and after applying the top x% cutoff, and a vector containing features with no fragments
#' @examples
#' # Example usage of the function
#' extract_raw_spectra(folder_path, l, 0.05, 0.7)
extract_raw_spectra <- function(folder_path, list_of_spectra_files, intrascan_grouping_tolerance = 0.05, top_percent_tic = 0.8) {
  stats_ms2 <- list_of_spectra_files$stats_file_ms2_only; spec_null = c()
  spectra_ms2 <- list_of_spectra_files$spectra_ms2_only
  sps_aggregate_all_mets <- spectra_ms2

  if (!dir.exists(file.path(folder_path, "MS2_scans_before_denoising"))) {
    dir.create(file.path(folder_path, "MS2_scans_before_denoising"))
  }

  # Initialize lists to store results
  sps_top80_tic_2 <- vector("list", length(sps_aggregate_all_mets))
  df <- data.frame(matrix(ncol = 3, nrow = length(sps_aggregate_all_mets)))
  colnames(df) <- c("Feature_ID", "Num_spectra_before_intrascan_grouping", "Num_spectra_before_intrascan_grouping")

  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(sps_aggregate_all_mets), style = 3)

  # Define the base folder path
  base_folder_path <- file.path(folder_path, "MS2_scans_before_denoising")

  # Use lapply with progress bar
  results <- lapply(seq_along(sps_aggregate_all_mets), function(i) {
    sps_j <- sps_aggregate_all_mets[[i]]
    name <- names(sps_aggregate_all_mets)[i]

    if (length(sps_j) != 0) {
      cat("Processing:", name, "\n")
      f_path <- file.path(base_folder_path, name)
      dir.create(f_path, recursive = TRUE, showWarnings = FALSE)

      # Process the scans without changing the working directory
      l_1 <- extract_features_from_scans_raw_data(sps_j, name, f_path, top_percent_tic, intrascan_grouping_tolerance)
      sps_top80_tic_2[[name]] <<- l_1[[1]]# Use `<<-` to assign to sps_top80_tic_2
      df[i, ] <<- c(name, l_1[[2]], l_1[[3]])  # Use `<<-` to assign to df

      # Update progress bar
      setTxtProgressBar(pb, i)

      return(NULL)  # No need to return anything specific
    } else {
      # Update progress bar even if the length is 0
      spec_null = c(spec_null, name)
      setTxtProgressBar(pb, i)
      return(name)
    }
  })

  # Close the progress bar
  close(pb)
  names(sps_top80_tic_2) =
  assign("sps_top_tic_2",sps_top80_tic_2, envir = .dures_env)
  # Return results as a list
  names(sps_top80_tic_2) <- df[,1]
  return(list(sps_top_tic_2 = sps_top80_tic_2, df = df, num_feats_no_spectra  = spec_null))
}
