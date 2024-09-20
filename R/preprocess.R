#' preprocess
#'
#' Reads mzml files from the predefined folder path and extracts spectra within the given mz and RT range.
#'
#' @param folder_path folder containing the required input directory (mzml) and feature list file in .txt format.
#' @param tol_mz mass tolerance
#' @param tol_rt RT tolerance
#' @return A list containing the MS2 spectra extracted from the mzml files using the predefined tolerances and the updated feature list file with MS2 features only
#' @examples
#' # Example usage of the function
#' preprocess(folder_path, 10, 0.1)
#' @export
preprocess <- function(folder_path, tol_mz = 5, tol_rt = 0.1667) {
  cat("*******************Sanity check for input files***********************\n")
  output_list = read_files(folder_path = folder_path, tol_mz, tol_rt)
  #set_user_path(folder_path)
  # Function to set user-defined path
  # set_user_path <- function(path) {
  #   .pkg_env$user_path <- path
  # }
  # #output_list[[1]] = contain spectral files
  #output_list[[2]] = contain ID, mz_up, mz_down, RT_median and sample-wise RT_min and RT_max (total = 6 + 2*number of mzml files with MS2 specta)
  #output_list[[3]] = contain names of mzml files used to generate output_list[[1]]
  sps_aggregate = concatenate_spectra(output_list)
  no_ms2_id = sps_aggregate$ids_with_no_ms2
  cat(paste("MS2 spectra couldn't be extracted from a total of ", length(no_ms2_id), " features out of ",length(sps_aggregate$sps_aggregate_all_mets), " features. They are as follows:",sep=""))
  cat(no_ms2_id, sep=",")
  stats_file = output_list[[2]]
  stats_file_ms2_triggered = stats_file[-c(which(stats_file$ID %in% no_ms2_id)), ]
  #update sps_aggregate
  #null_entries <- sapply(sps_aggregate$sps_aggregate_all_mets, is.null)
  sps_aggregate_ms2_triggered = Filter(Negate(is.null), sps_aggregate$sps_aggregate_all_mets)

  return(list(spectra_ms2_only = sps_aggregate_ms2_triggered, stats_file_ms2_only = stats_file_ms2_triggered))
 }
