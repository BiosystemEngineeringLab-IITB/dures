#' tune_parameters
#'
#' Helper function for generate_denoise_spectra(). Need not be called separately
#' @param agg_df: spectrum dataframe with columns Fragment_ID, MZ, Intensity & Freq
#' @param threshold: see function description for generate_denoised_spectra
#' @param name: names of the spectra such as "scan_1.txt" "scan_2.txt" and so on
#' @param path: folder path to store the denoised spectra in .txt and .mzML format
#' @return directly prints the spectra to the file in .txt format and returns the spectrum-wise mz and intensity
#' @examples
#' # Example usage of the function
#' tune_parameters(agg_df, threshold, "scan_1.txt", path)
#' @export
tune_parameters <- function(agg_df, threshold, name, path){
  colnames(agg_df)[2:4] = c("Mean_MZ", "Mean_Intensity", "Frequency")
  agg_df_cutoff = subset(agg_df, agg_df$Frequency >= threshold)
  agg_df_cutoff = agg_df_cutoff[order(agg_df_cutoff$Mean_MZ, decreasing = F),]
  # #creating the composite spectra
  # spd <- DataFrame(
  #   msLevel = c(2L),
  #   polarity = c(1L),
  #   id = c(name),
  #   name = c(name))
  # spd$mz <- list(
  #   c(agg_df_cutoff$Mean_MZ))
  # spd$intensity <- list(
  #   c(agg_df_cutoff$Mean_Intensity))
  # sps <- Spectra(spd)

  #plotSpectra(sps, main = paste(sps$name, " (Freq_cutoff_", as.character(threshold),")", sep=""))
  # folder_create = dir.create(paste("~/metabolomics/WTC+DIELdata/results_ppm=5/Denoised_spectra/",name,"/", sep=""))
  #setwd(paste("~/metabolomics/Sportomics/Database_metabolites/Aggregate_top_5_files_SS_score/Files_of_interest_Sportomic/After_filtering_stage_2_spectra/", name,sep=""))
  name_exp= paste(path, name, sep="")
  #sps_df=data.frame(fragments = peaksData(sps)[[1]][,1], intensity =  peaksData(sps)[[1]][,2])
  sps_df = data.frame(fragments = agg_df_cutoff$Mean_MZ, intensity = agg_df_cutoff$Mean_Intensity)
  write.table(sps_df, name_exp, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  return(list(agg_df_cutoff$Mean_MZ, agg_df_cutoff$Mean_Intensity))
}
