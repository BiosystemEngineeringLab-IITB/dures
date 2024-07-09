#' generate_denoised_spectrum
#'
#' Takes in a list of individual spectrum with fragments labeled with frequencies, denoises it using fixed cutoff and writes the resulting spectrum into a file
#' @param aggregate_list: output from label_individual_spectrum()
#' @param folder_path: same path as the one used in preprocess()
#' @param custom_threshold: must be between 0 and 1 (default = NULL). If the user wants to apply a threshold different from ours (0.12 (95%CI- 0.08-0.15)) they can do so using this parameter. The threshold when multiplied by the number of replicate spectra (for a given feature) will tell us that only fragments with frequencies above this value will be considered as signal fragments and retained in the spectrum
#' @return directly prints the spectra to the file in .txt and .mzML format
#' @examples
#' # Example usage of the function
#' generate_denoised_spectra(l3, folder_path, NULL)
generate_denoised_spectra <- function(aggregate_list, folder_path, custom_threshold = NULL){
  if (!dir.exists(file.path(folder_path, "Denoised_spectra"))) {
    dir.create(file.path(folder_path, "Denoised_spectra/"))
  }
  sps_top_tic_2 <- get_sps_top_tic_2()
  freq_df = aggregate_list
  for(j in 1:length(freq_df)){
    if (!dir.exists(file.path(paste(folder_path,names(freq_df)[j],"/",sep="")))) {
      dir.create(file.path(paste(folder_path,"Denoised_spectra/",names(freq_df)[j],"/",sep="")))
      path = file.path(paste(folder_path,"Denoised_spectra/", names(freq_df)[j],"/",sep=""))
    }
    idx = which(names(sps_top_tic_2) %in% names(freq_df)[j])
    #print(length(sps_top_tic_2[[idx]]))
    if(length(sps_top_tic_2[[idx]]) <= 25){
      threshold = (3/length(sps_top_tic_2[[idx]]))
    } else if(length(sps_top_tic_2[[idx]]) >= 416){
      threshold = (50/length(sps_top_tic_2[[idx]]))
    } else if(length(custom_threshold)!=0) {
      threshold = custom_threshold
    }else{
      threshold = 0.12
    }
    #file_names = paste(names(freq_df)[[j]], names(freq_df[[j]]), sep="/")
    file_names = names(freq_df[[j]])

    for(k in 1:length(freq_df[[j]])){
      f = freq_df[[j]][[k]]
      sps_tic_20_all_freq = tune_parameters(f, threshold, file_names[k], path)

    }
  }

}


