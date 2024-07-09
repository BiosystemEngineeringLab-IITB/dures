#' extract_features_from_scans_raw_data
#'
#' Helper function for extract_raw_spectra. Need not be run separately. For every feature, reads the MS2 spectra, groups fragments within a given tolerance, extracts the resulting sprectrum to a text file. Based on user input, it will retain only the top x% TIC spectra.
#'
#' @param sps_j: a list of concatenated MS2 spectra for a given feature
#' @param name: name of feature (ID)
#' @param f_path: path to the a folder with the same name as the feature ID is created inside the folder titled MS2_scans_before_denoising
#' @param cutoff: top x% TIC cutoff
#' @param tol: mass tolerance for grouping
#' @return A list containing the top x% TIC and number of scans before and after TIC filtering
#' @examples
#' # Example usage of the function
#' extract_features_from_scans_raw_data(sps_j, 1098, f_path, 0.8, 0.05)
#' @export

extract_features_from_scans_raw_data <- function(sps_j,name,f_path,cutoff,tol){
  mz=c(); inten=c(); name1=c(); scan =c(); tic=c(); nfg=c(); scan_info = c()
  for(i in 1:length(sps_j)){
    sps_i = sps_j[i]
    fls <- unique(Spectra::dataOrigin(sps_i))
    #n = paste(strsplit(strsplit(unique(dataOrigin(sps_i)), "/")[[1]][8], " ")[[1]], collapse = "_")
    sps_peaks = Spectra::peaksData(sps_i)[[1]]
    nf = dim(sps_peaks)[1]
    tic_1 = sum(sps_peaks[,2])
    colnames(sps_peaks)[1] = "fragments"
    mz = c(mz, sps_peaks[,1])
    inten = c(inten, sps_peaks[,2])
    #n1 = rep(n, dim(sps_peaks)[1])
    nf = rep(nf, dim(sps_peaks)[1])
    nfg = c(nfg, nf)
    sc = sps_i$spectrumId
    sc_1 = rep(sc, dim(sps_peaks)[1])
    tic_2 = rep(tic_1, dim(sps_peaks)[1])
    tic=c(tic, tic_1)
    scan = c(scan, sc_1)
    scan_info = c(scan_info, sc)
    #name1 = c(name1,  n1)
    file = paste(f_path,"/",sc,".txt",sep="")
    sps_peaks = as.data.frame(sps_peaks)
    sps_peaks_intrascan_grouped = intrascan_grouping(sps_peaks,tol)
    write.table(sps_peaks_intrascan_grouped, file, col.names = T, row.names = F, sep="\t", quote = F )
  }
  tic_20_df = data.frame(Scan = scan_info, TIC = tic)
  ranked = round(length(sps_j)*cutoff)
  tic_20_df = tic_20_df[order(tic_20_df$TIC, decreasing = T),]
  #print("Number of scans before TIC cutoff applied")
  n_scans_before = dim(tic_20_df)[1]
  tic_20_df = tic_20_df[c(1:ranked), ]
  #print("Number of scans after TIC cutoff applied")
  n_scans_after = dim(tic_20_df)[1]
  sps_new = sps_j[which(sps_j$spectrumId %in% tic_20_df$Scan)]
  #print("Number of scans left for aggregating")
  #print(length(sps_new))
  full_paths <- file.path(f_path, setdiff(list.files(f_path), paste(sps_new$spectrumId,".txt",sep="")))
  deleted = file.remove(full_paths)
  #print(length(list.files(f_path)))
  return(list(sps_new, n_scans_before, n_scans_after))
}
