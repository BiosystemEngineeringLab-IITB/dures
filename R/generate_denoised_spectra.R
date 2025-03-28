#' generate_denoised_spectrum
#'
#' Takes in a list of individual spectrum with fragments labeled with frequencies, denoises it using fixed cutoff and writes the resulting spectrum into a file
#' @param aggregate_list: output from label_individual_spectrum()
#' @importFrom S4Vectors DataFrame
#' @importFrom Spectra export
#' @importFrom Spectra MsBackendMzR
#' @param folder_path: same path as the one used in preprocess()
#' @param ion_mode: ionization mode = "pos" or "neg" (default = "pos")
#' @param custom_threshold: must be between 0 and 1 (default = NULL). If the user wants to apply a threshold different from ours (0.12 (95%CI- 0.08-0.15)) they can do so using this parameter. The threshold when multiplied by the number of replicate spectra (for a given feature) will tell us that only fragments with frequencies above this value will be considered as signal fragments and retained in the spectrum
#' @return directly prints the spectra to the file in .txt and .mzML format
#' @examples
#' # Example usage of the function
#' generate_denoised_spectra(l3, folder_path, NULL)
#' @export
generate_denoised_spectra <- function(aggregate_list, folder_path, custom_threshold = NULL, ion_mode = "pos"){
  if (!dir.exists(file.path(folder_path, "Denoised_spectra_mzML"))) {
    dir.create(file.path(folder_path, "Denoised_spectra_mzML/"))
  }
  if (!dir.exists(file.path(folder_path, "Denoised_spectra_txt"))) {
    dir.create(file.path(folder_path, "Denoised_spectra_txt/"))
  }
  path_mzML = file.path(folder_path, "Denoised_spectra_mzML/")
  path_txt = file.path(folder_path, "Denoised_spectra_txt/")
  sps_top_tic_2 <- get_sps_top_tic_2(); nm_s = list()
  freq_df = aggregate_list

  mz = list(); inten = list()

  for(j in 1:length(freq_df)){
    idx = which(names(sps_top_tic_2) %in% names(freq_df)[j])
    #   if(length(sps_top_tic_2[[idx]]) <= 25){
    #     threshold = (3/length(sps_top_tic_2[[idx]]))
    #   } else if(length(sps_top_tic_2[[idx]]) >= 416){
    #     threshold = (50/length(sps_top_tic_2[[idx]]))
    #   } else if(length(custom_threshold)!=0) {
    #     threshold = custom_threshold
    #   }else{
    #     threshold = 0.12
    #   }
      #for feature number j, I have m scans belonging to n samples
    name_samples = unlist(lapply(names(freq_df[[j]]), function(x) strsplit(x,"_scan")[[1]][1]))
    name_scans = unlist(lapply(names(freq_df[[j]]), function(x) strsplit(x,"_scan_")[[1]][2]))
    nm_s[[j]] = name_samples

    #sps_tic_20_all_freq = list();
    for(k in 1:length(freq_df[[j]])){
      f = freq_df[[j]][[k]]
      colnames(f)[2:4] = c("Mean_MZ", "Mean_Intensity", "Frequency")
      f_cutoff = subset(f, f$Frequency >= custom_threshold)
      f_cutoff = f_cutoff[order(f_cutoff$Mean_MZ, decreasing = F),]
      #name_exp= paste(path, name, sep="")
      #sps_df = data.frame(fragments = f_cutoff$Mean_MZ, intensity = f_cutoff$Mean_Intensity)
      #print(f_cutoff)
      #mz[[paste(as.character(k), "_", names(freq_df)[j],"_",name_samples[k],sep="")]] = f_cutoff$Mean_MZ
      #inten[[paste(as.character(k), "_", names(freq_df)[j],"_",name_samples[k],sep="")]] = f_cutoff$Mean_Intensity
      mz[[paste(names(freq_df)[j],"_",name_samples[k],"_scan_", name_scans[k], sep="")]] = f_cutoff$Mean_MZ
      inten[[paste(names(freq_df)[j],"_",name_samples[k],"_scan_", name_scans[k], sep="")]] = f_cutoff$Mean_Intensity
      sps_df = data.frame(fragments = f_cutoff$Mean_MZ, intensity = f_cutoff$Mean_Intensity)
      dir_path = paste(path_txt, "/", names(freq_df)[j],"/", sep="")
      if (!dir.exists(dir_path)) {
        dir.create(dir_path)
      }
      write.table(sps_df, paste(dir_path, paste(names(freq_df)[j],"_",name_samples[k],"_scan_", name_scans[k], sep=""), sep=""), sep="\t", col.names = T, row.names = F, quote = F)


    }
  }

  for(w in 1:length(unique(unlist(nm_s)))){
    samp = unique(unlist(nm_s))[w]
    if(ion_mode == "pos"){
      pol = 1L
    }else{
      pol = 0L
    }
    MZ = mz[grep(samp, names(mz))]
    INTEN = inten[grep(samp, names(mz))]
    spd <- DataFrame(
      msLevel = c(rep(2L, length(MZ))),
      polarity = c(rep(pol, length(MZ))),
      # id = c(names(freq_df)[j]),
      name = c(samp))

    spd$mz<- MZ
    spd$intensity <- INTEN
    sps <- Spectra::Spectra(spd)
    sps$spectrumId = names(mz)[grep(samp, names(mz))]
    fl = paste(path_mzML, "/", samp, sep="")
    export(sps, MsBackendMzR(), file = fl)

  }


  }

