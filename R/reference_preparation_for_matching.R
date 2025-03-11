#' reference_preparation_for_matching
#'
#' Prepares the reference library for matching with current experimental spectra. This neednot be run independently. This is a companion function to precursor_matching.R
#'
#' @param j index
#' @param annotation_file stats file with ID, mz, RT etc
#' @param  lib positive or negative library
#' @param path folder_path to store the spectra as txt files
#' @return reference spectra in Spectra format
#' @examples
#' # Example usage of the function
#' reference_preparation_for_matching(j, l1$stats_file_ms2_only, lib, p)
#' @export

reference_preparation_for_matching <- function(j, annotation_file, lib, path) {
  # Extract relevant information
  name <- annotation_file$ID[j]
  left <- annotation_file$mz_down[j]
  right <- annotation_file$mz_up[j]

  # Find matches based on m/z tolerance
  ind <- which(dplyr::between(lib$mz, left, right))

  if (length(ind) == 0) {
    message("No reference match using the current m/z tolerance for index: ", j)
    return(NULL)
  } else {
    # Extract relevant library data
    library_names <- lib$library_name[ind]
    splash_names <- lib$splash_keys[ind]
    mz <- lib$mz[ind]

    # Prepare references
    sps_ref_new <- vector("list", length(ind))
    for (i in seq_along(ind)) {
      sps_ref_new[[i]] <- prepare_reference(ind[i], paste(mz[i], library_names[i], splash_names[i], sep = "_"), lib, path)
    }
    return(sps_ref_new)
  }
}
