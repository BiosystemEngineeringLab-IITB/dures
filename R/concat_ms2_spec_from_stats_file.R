#' concat_ms2_spec_from_stats_file
#'
#' Helper function for concatenate_spectra(). Need not be called separately. Gets auto-invoked once we call concatenate_spectra()
#'
#' @param j: jth feature from the stats file
#' @param stats_file: edited stats file with RT_min and RT_max prefixes added to every sample
#' @param fl: list of sample names
#' @return a list with all spectra concatenated for a particular feature
#' @examples
#' # Example usage of the function (need not be called separately)
#' concat_ms2_spec_from_stats_file(1,stats_file,fl)
#' @export

concat_ms2_spec_from_stats_file <- function(j, stats_file, fl){

  #this function extracts all scans within a given RT range
  #and mz tolerance (decided by ppm error = Dalton error/mz * 10^6)
  #i = index of metabolite in the annotation_file_final dataframe
  #fl = names(spectral_files)
  sps_comb = list()
  fl = unlist(lapply(fl, function(x) strsplit(x, ".mzML")[[1]][1]))
  for (i in 1:length(fl)) {
    # Create the prefix pattern
    pattern <- paste0("RT_(min|max)_", fl[i], "\\b")

    # Find column names that match the pattern
    matching_columns <- grep(pattern, names(stats_file), value = TRUE)

    # Find the column indices
    column_indices <- which(names(stats_file) %in% matching_columns)

    rtime = c(stats_file[, column_indices[1]][j], stats_file[, column_indices[2]][j])

    if(rtime[1] == 0 & rtime[2] ==0){
      next
    }
    else{
      spectral_files_f <- get_spectral_files()
      sps = Spectra::filterRt(spectral_files_f[[paste(fl[i],".mzML",sep="")]], rt =  c(rtime[1], rtime[2]))
      sps_mz= Spectra::filterPrecursorMzRange(sps, mz=c(stats_file$mz_down[j],
                                                       stats_file$mz_up[j]))
      sps_comb[[i]] = sps_mz
    }
  }

  sps_aggregate = Spectra::concatenateSpectra(sps_comb)

  flag =0

  #just return a list with all scans for a particular metabolite
  if(length(sps_aggregate)!=0){
    #sps_union = create_consensus_table(sps_aggregate, "IDA_1", annotation_file_final$Inchi[j])
    sps_return = list(Aggregate = sps_aggregate)
    flag = 1

  }

  if(flag!=0){
    return(sps_return)
  }
}


