#' concat_ms2_spec_from_stats_file
#'
#' Helper function for concatenate_spectra(). Need not be called separately. Gets auto-invoked once we call concatenate_spectra()
#'
#' @param j jth feature from the stats file.
#' @param stats_file edited stats file with RT_min and RT_max prefixes added to every sample.
#' @param fl list of sample names.
#' @return A list with all spectra concatenated for a particular feature
#' @examples
#' # Example usage of the function (need not be called separately)
#' concat_ms2_spec_from_stats_file(1,stats_file,fl)
#' @export
concat_ms2_spec_from_stats_file <- function(j, stats_file, fl) {
  # This function extracts all MS2 scans within a given RT range and m/z tolerance window
  # for a given metabolite index 'j' from a stats_file.
  # It supports both mzML and mzXML input formats.

  # Initialize list to store per-file spectra
  sps_comb <- list()

  # Clean file names (remove extensions like .mzML or .mzXML)
  fl_clean <- sub("\\.(mzML|mzXML)$", "", fl)

  # Get all spectral files loaded into environment
  spectral_files_f <- get_spectral_files()

  for (i in seq_along(fl_clean)) {
    sample_name <- fl_clean[i]

    # Find matching RT columns
    pattern <- paste0("RT_(min|max)_", sample_name, "\\b")
    matching_columns <- grep(pattern, names(stats_file), value = TRUE)
    column_indices <- which(names(stats_file) %in% matching_columns)

    # Skip if matching RT columns are not found
    if (length(column_indices) != 2) next

    rtime <- c(stats_file[j, column_indices[1]], stats_file[j, column_indices[2]])

    # Skip if both RTs are 0
    if (all(rtime == 0)) next

    # Try to find matching mzML or mzXML in loaded spectral_files
    spec_file_key <- NULL
    if (!is.null(spectral_files_f[[paste0(sample_name, ".mzML")]])) {
      spec_file_key <- paste0(sample_name, ".mzML")
    } else if (!is.null(spectral_files_f[[paste0(sample_name, ".mzXML")]])) {
      spec_file_key <- paste0(sample_name, ".mzXML")
    } else {
      next  # No matching spectral file
    }

    # Filter by RT and then by precursor m/z range
    sps <- Spectra::filterRt(spectral_files_f[[spec_file_key]], rt = rtime)
    sps_mz <- Spectra::filterPrecursorMzRange(sps, mz = c(stats_file$mz_down[j], stats_file$mz_up[j]))

    if (length(sps_mz) > 0) {
      sps_comb[[length(sps_comb) + 1]] <- sps_mz
    }
  }

  # Combine all matched spectra across files
  if (length(sps_comb) == 0) return(NULL)

  sps_aggregate <- Spectra::concatenateSpectra(sps_comb)

  if (length(sps_aggregate) > 0) {
    return(list(Aggregate = sps_aggregate))
  } else {
    return(NULL)
  }
}


