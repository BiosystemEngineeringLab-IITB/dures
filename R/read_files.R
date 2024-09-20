#' read_files
#'
#' Helper function for preprocess(). Need not be run separately. Runs several checks to see if features file has required columns or if user-defined features file has been provided
#' Two options are available for inputting feature list from the user. First, the file titled feature_list.txt must be present in the folder_path with the columns "ID", "mz" and "RT". Second, if an additional file titled Sample-and-feature-wise-RT-tolerance.txt is present, it should contain the same number of rows and double the number of columns (RT_min_ and RT_max_ for every sample) as feature_list.txt.
#' @param folder_path: same folder path used in preprocess(). The extracted fragment-grouped spectra in .txt format will be stored here
#' @importFrom utils read.delim
#' @param tol_mz: mass tolerance (default: 0.05 Da)
#' @param tol_rt: Rt tolerance (default: 0.1667 minutes)
#' @return a list containing all MS2 spectra for every sample, an edited feature list file for subsequent steps and names of samples analyzed
#' @examples
#' # Example usage of the function
#' read_files(folder_path, 0.05, 0.1667)
read_files <- function(folder_path, tol_mz, tol_rt){
  #tolerance in minutes and mz_tolerance in ppm
  if(length(grep("feature_list.txt", list.files(folder_path))) > 0)
    {
    cat("\u2713 Feature list file found!\n")
    fil = read.delim(paste(folder_path, list.files(folder_path)[grep("feature_list.txt", list.files(folder_path))], sep=""), check.names = FALSE, comment.char = "#")
    # Assuming your dataframe is named df
    required_columns <- c("ID", "mz", "RT")

    # Check if the dataframe has exactly the required columns
    if (all(required_columns %in% colnames(fil)) && length(colnames(fil)) == length(required_columns)) {
      cat("\u2713 Feature list file has the correct columns\n")
    } else {
      stop("Error: Feature list file does not have the exact required columns.")
    }
  }else{
    stop("Feature list file not found!")
  }

  fil$ID = as.character(fil$ID)
  if(dir.exists(paste(folder_path, "mzml_files/",sep=""))){
    cat("\u2713 mzml folder exists\n")
    fls = list.files(paste(folder_path, "mzml_files/",sep=""))
  }else{
    stop("mzml folder not found. pls check if you have named your folder mzml_files/")
  }

  if(file.exists(paste(folder_path, "Sample-and-feature-wise-RT-tolerance.txt",sep=""))){
    cat("\u2713 RT tolerance file exists! Using custom tolerance\n")
    RT_tolerance = read.delim(paste(folder_path, "Sample-and-feature-wise-RT-tolerance.txt",sep=""),check.names = FALSE)
    if(dim(RT_tolerance)[2] == 2*length(fls)){
      cat("\u2713 Input file is in correct format\n")
    }
    else{
      stop("Input file should contain RT_min and RT_max for all files present in the mzml folder")
    }

  }else{
    cat(paste("RT tolerance file not found. Proceeding with the default mz tolerance of: ", tol_mz, " ppm and RT tolerance of: ", tol_rt, " minutes for all samples\n", sep = ""))
    RT_tolerance = matrix(0, nrow = dim(fil)[1], ncol = length(fls)*2)
    RT_tolerance[, c(1:length(fls))] = fil$RT - tol_rt
    RT_tolerance[, c((length(fls)+1) : dim(RT_tolerance)[2])] = fil$RT + tol_rt

    colnames(RT_tolerance) = c(paste("RT_min", unlist(lapply(fls, function(x) strsplit(x,".mzML")[[1]][1])), sep="_"), paste("RT_max", unlist(lapply(fls, function(x) strsplit(x,".mzML")[[1]][1])), sep="_"))
  }
  #converting RT start and end to seconds
  RT_tolerance = RT_tolerance * 60
  #fil$mz_tolerance = (fil$mz * ppm)/1000000
  fil$mz_up = fil$mz + (fil$mz * tol_mz)/1000000
  #this tolerance for m/z for MS2 is specific to the WTC dataset and the instrument used to sequence it
  fil$mz_down = fil$mz - (fil$mz * tol_mz)/1000000
  fil_rt = cbind(fil, RT_tolerance)
  # n_rows = dim(fil)[1]
  # n_cols = length(fls)*2 + 1
  # column_names = c("ID", paste("RT_min_", RT_tolerance$file_name,sep=""), paste("RT_max_", RT_tolerance$file_name,sep=""))
  # fil_rt <- data.frame(matrix(0, nrow = n_rows, ncol = n_cols))
  # colnames(fil_rt) = column_names
  # fil_rt$ID = fil$ID
  # for(j in 1:((n_cols-1)/2)){
  #   idx = grep(RT_tolerance$file_name[j], colnames(fil_rt))
  #   fil_rt[, idx[1]] = c(fil$RT - RT_tolerance$RT_tolerance[j])
  #   fil_rt[, idx[2]] = c(fil$RT + RT_tolerance$RT_tolerance[j])
  # }
  # fil_rt$ID = NULL
  # fil_combined = cbind(fil, fil_rt)

  cat("Reading mzml files\n")
  # Pre-define the path to the mzML folder
  path_to_mzml_folder <- paste(folder_path, "mzml_files/",sep="")

  # Initialize variables
  flag <- character(length(fls))  # Preallocate flag as a character vector
  spectral_files <- vector("list", length(fls))  # Preallocate the list
  flag_count <- 1  # Counter for flag entries

  # Initialize the progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(fls), style = 3)
  backend <- Spectra::MsBackendMzR()
  # Loop over the files
  for (i in seq_along(fls)) {
    #print(i)
    sps <- Spectra::Spectra(file.path(path_to_mzml_folder, fls[i]), source = backend)
    sps = Spectra::filterMsLevel(sps, msLevel = 2)
    if(length(sps) != 0){
      #sps$spectrumId = paste("scan_",seq(1,length(sps)), sep="")
      spectral_files[[i]] <- sps # Filtering only MS2 level spectra
    }


    if (length(spectral_files[[i]]) == 0) {
      flag[flag_count] <- fls[i]
      flag_count <- flag_count + 1
    }

    # Update the progress bar
    utils::setTxtProgressBar(pb, i)
  }

  # if(sapply(spectral_files, length) >3000){
  #   cat("yay")
  # }


  fls_f = unlist(lapply(fls, function(x) strsplit(x,".mzML")[[1]][1]))
  idx_ms2_absent = which(sapply(spectral_files, length) == 0)
  f=0
  if(length(idx_ms2_absent) > 0){
    f = 1
    cat("\n MS2 data are not present in the following files:\n")
    print(paste(fls_f[idx_ms2_absent], collapse = ","))
    pattern <- paste0("RT_(min|max)_(", paste(fls_f[idx_ms2_absent], collapse = "|"), ")")
    # Identify columns to remove
    cols_to_remove <- grepl(pattern, colnames(fil_rt))
    spectral_files_filtered = spectral_files[-c(idx_ms2_absent)]
    fil_rt_filtered <- fil_rt[, !cols_to_remove]
    fls_filtered = fls[-c(idx_ms2_absent)]
    names(spectral_files_filtered) = fls_filtered
  }else{
    names(spectral_files) = fls
  }

  cat("\n \u2713 Spectra has been extracted from the mzml files\n")
  if(dir.exists(paste(folder_path, "R_objects/",sep="")) == FALSE){
    dir.create(paste(folder_path, "R_objects/",sep=""))
  }

  if(f == 0){
    for(l in 1:length(spectral_files)){
      spectral_files[[l]]$spectrumId = paste(names(spectral_files)[l],"_scan_", 1:length(spectral_files[[l]]),sep="")
    }
    #spectral_files_f <<- spectral_files
    assign("spectral_files_f", spectral_files, envir = .dures_env)
    saveRDS(spectral_files, paste(folder_path, "R_objects/", "spectral_files.rds",sep=""))
    return(list(spectral_files = spectral_files, stats_file_analyzed = fil_rt, file_names = fls))
  }
  else{
    for(l in 1:length(spectral_files_filtered)){
      spectral_files_filtered[[l]]$spectrumId = paste(names(spectral_files_filtered)[l],"_scan_", 1:length(spectral_files_filtered[[l]]),sep="")
    }
    assign("spectral_files_f", spectral_files_filtered, envir = .dures_env)
    #spectral_files_f <<- spectral_files_filtered
    saveRDS(spectral_files_filtered, paste(folder_path, "R_objects/", "spectral_files.rds",sep=""))
    return(list(spectral_files = spectral_files_filtered, stats_file_analyzed = fil_rt_filtered, file_names = fls_filtered))
    }

}
