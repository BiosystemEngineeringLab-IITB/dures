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
  if (!dir.exists(file.path(folder_path, "Denoised_spectra"))) {
    dir.create(file.path(folder_path, "Denoised_spectra/"))
  }
  path = file.path(folder_path, "Denoised_spectra/")
  sps_top_tic_2 <- get_sps_top_tic_2(); nm_s = list()
  freq_df = aggregate_list

  mz = list(); inten = list()

  for(j in 1:length(freq_df)){
    idx = which(names(sps_top_tic_2) %in% names(freq_df)[j])
    if(length(sps_top_tic_2[[idx]]) <= 25){
      threshold = (3/length(sps_top_tic_2[[idx]]))
    } else if(length(sps_top_tic_2[[idx]]) >= 416){
      threshold = (50/length(sps_top_tic_2[[idx]]))
    } else if(length(custom_threshold)!=0) {
      threshold = custom_threshold
    }else{
      threshold = 0.12
    }
    #for feature number j, I have m scans belonging to n samples
    name_samples = unlist(lapply(names(freq_df[[j]]), function(x) strsplit(x,"_scan")[[1]][1]))
    name_scans = unlist(lapply(names(freq_df[[j]]), function(x) strsplit(x,"_scan_")[[1]][2]))
    nm_s[[j]] = name_samples

    #sps_tic_20_all_freq = list();
    for(k in 1:length(freq_df[[j]])){
      f = freq_df[[j]][[k]]
      colnames(f)[2:4] = c("Mean_MZ", "Mean_Intensity", "Frequency")
      f_cutoff = subset(f, f$Frequency >= threshold)
      f_cutoff = f_cutoff[order(f_cutoff$Mean_MZ, decreasing = F),]
      #name_exp= paste(path, name, sep="")
      #sps_df = data.frame(fragments = f_cutoff$Mean_MZ, intensity = f_cutoff$Mean_Intensity)
      #print(f_cutoff)
      #mz[[paste(as.character(k), "_", names(freq_df)[j],"_",name_samples[k],sep="")]] = f_cutoff$Mean_MZ
      #inten[[paste(as.character(k), "_", names(freq_df)[j],"_",name_samples[k],sep="")]] = f_cutoff$Mean_Intensity
      mz[[paste(names(freq_df)[j],"_",name_samples[k],"_scan_", name_scans[k], sep="")]] = f_cutoff$Mean_MZ
      inten[[paste(names(freq_df)[j],"_",name_samples[k],"_scan_", name_scans[k], sep="")]] = f_cutoff$Mean_Intensity


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

    fl = paste(path, samp, sep="")
    export(sps, MsBackendMzR(), file = fl)

  }


  }

    # if(length(freq_df[[j]]) >=3){
    #   idx = which(names(sps_top_tic_2) %in% names(freq_df)[j])
    #   if(length(sps_top_tic_2[[idx]]) <= 25){
    #     threshold = (3/length(sps_top_tic_2[[idx]]))
    #   } else if(length(sps_top_tic_2[[idx]]) >= 416){
    #     threshold = (50/length(sps_top_tic_2[[idx]]))
      # } else else{
      #   threshold = 0.12
      # }

    #for every feature calculate the optimal threshold

    #print(length(sps_top_tic_2[[idx]]))

  # #try to create the one folder with just the MS2 spectra inside original samples in mzml format.
  # #for now no need to incorporate the MS1 spectra. Just the MS2 spectra and that too only for those samples which have MS2 in them
  # for(j in 1:length(freq_df)){
  #   if (!dir.exists(file.path(paste(folder_path,names(freq_df)[j],"/",sep="")))) {
  #     dir.create(file.path(paste(folder_path,"Denoised_spectra/",names(freq_df)[j],"/",sep="")))
  #     path = file.path(paste(folder_path,"Denoised_spectra/", names(freq_df)[j],"/",sep=""))
  #   }
  #   idx = which(names(sps_top_tic_2) %in% names(freq_df)[j])
  #   #print(length(sps_top_tic_2[[idx]]))
  #   if(length(sps_top_tic_2[[idx]]) <= 25){
  #     threshold = (3/length(sps_top_tic_2[[idx]]))
  #   } else if(length(sps_top_tic_2[[idx]]) >= 416){
  #     threshold = (50/length(sps_top_tic_2[[idx]]))
  #   } else if(length(custom_threshold)!=0) {
  #     threshold = custom_threshold
  #   }else{
  #     threshold = 0.12
  #   }
  #   #file_names = paste(names(freq_df)[[j]], names(freq_df[[j]]), sep="/")
  #   file_names = names(freq_df[[j]]);sps_tic_20_all_freq = list(); mz = list(); inten = list()
  #
  #
  #   for(k in 1:length(freq_df[[j]])){
  #     f = freq_df[[j]][[k]]
  #     sps_tic_20_all_freq = tune_parameters(f, threshold, file_names[k], path)
  #     mz[[k]] = sps_tic_20_all_freq[[1]]
  #     inten[[k]] = sps_tic_20_all_freq[[2]]
  #   }



