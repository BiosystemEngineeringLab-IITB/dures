#' derive_aggregate_spectra
#'
#' Helper code to run derive_aggregate_spectra
#' @param spectra_list: list of spectra objects for all features, each containing top x% TIC spectra concatenated
#' @param mz_tol: mass tolerance (in Da) required when grouping fragments across multiple spectra for a given feature
#' @param folder_path: the path to the folder where you want to store the data reduction dataframe
#' @return a list containing the aggregate spectra object and a dataframe of the aggregate spectra for all features
#' @examples
#' # Example usage of the function
#' call_aggregate(sps_top_tic, 0.05)
#' @export
call_aggregate <- function(spectra_list, mz_tol, folder_path){
  df = data.frame()
  print(paste("Creating aggregate spectra for ", length(spectra_list), " features. This is an expensive operation. Please wait..", sep=""))
  result_list <- pbapply::pblapply(spectra_list, derive_aggregate_spectra, mz_tol)
  dims_df <- data.frame(
    num_fragments_before_grouping = sapply(result_list, function(res) res[[3]]),
    num_fragments_after_grouping = sapply(result_list, function(res) res[[4]])
  )
  write.csv(dims_df, paste(folder_path, "num_fragments_before_after_grouping.csv", sep=""), row.names = FALSE)

  return(result_list)
}
