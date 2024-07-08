concatenate_spectra <- function(output_list) {
  stats_file <- output_list[[2]]
  id_no_ms2 <- c()

  # Initialize the list to store results
  sps_aggregate_all_mets <- vector("list", dim(stats_file)[1])

  # Use pbapply's pblapply for parallel processing with progress bar
  sps_aggregate_all_mets <- pbapply::pblapply(seq_len(dim(stats_file)[1]), function(i) {
    print(paste("Iteration", i))

    # Concatenate spectra
    l <- dures::concat_ms2_spec_from_stats_file(i, stats_file, output_list[[3]])

    if (is.null(l)) {
      id_no_ms2 <<- c(id_no_ms2, stats_file$ID[i])
      print(paste("No MS2 spectra could be extracted for feature ", stats_file$ID[i], " for the given mz and RT range", sep = ""))
      return(NULL)
    }

    # l is a list containing individual scans and combined (union) of all scans
    aggregate <- l$Aggregate
    aggregate$spectrumId <- paste("scan", seq_along(aggregate), sep = "_")

    gc()

    return(list(aggregate = aggregate))
  })

  names(sps_aggregate_all_mets) = stats_file$ID

  return(list(sps_aggregate_all_mets = sps_aggregate_all_mets, ids_with_no_ms2 = id_no_ms2))
}
