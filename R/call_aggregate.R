#' call_aggregate
#'
#' Helper function to run derive_aggregate_spectra.
#'
#' @param spectra_list A list of spectra objects for all features, each containing the top x% TIC spectra concatenated.
#' @param mz_tol Mass tolerance (in Da) required when grouping fragments across multiple spectra for a given feature.
#' @param folder_path Path to the folder where the data reduction dataframe will be stored.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{aggregate_spectra_object}}: Aggregate spectra object.
#'   \item{\code{aggregate_dataframe}}: Dataframe of aggregate spectra for all features.
#' }
#'
#' @examples
#' folder_path <- "path/to/folder"
#' sp <- preprocess(folder_path, 5, 0.01)
#' l <- sp[[1]]
#' sp1 <- extract_raw_spectra(folder_path, l, 0.05, 0.8)
#' spectra_list <- sp1[[1]]
#' call_aggregate(spectra_list, 0.05, folder_path)
#'
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
