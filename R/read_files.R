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
read_files <- function(folder_path, tol_mz, tol_rt) {
  # 1. Read and validate feature list
  feature_list_path <- list.files(folder_path, pattern = "feature_list.txt", full.names = TRUE)
  if (length(feature_list_path) == 0) stop("Feature list file not found!")
  cat("\u2713 Feature list file found!\n")

  fil <- read.delim(feature_list_path, check.names = FALSE, comment.char = "#")
  required_columns <- c("ID", "mz", "RT")
  if (!all(required_columns %in% colnames(fil)) || length(colnames(fil)) != length(required_columns)) {
    stop("Error: Feature list file does not have the exact required columns.")
  }
  cat("\u2713 Feature list file has the correct columns\n")
  fil$ID <- as.character(fil$ID)

  # 2. Validate mzML/mzXML folder and list files
  mzml_path <- file.path(folder_path, "mzml_files")
  if (!dir.exists(mzml_path)) stop("mzml folder not found. Please name your folder 'mzml_files/'")
  cat("\u2713 mzml folder exists\n")

  fls <- list.files(mzml_path, pattern = "\\.(mzML|mzXML)$", full.names = FALSE)
  if (length(fls) == 0) stop("No .mzML or .mzXML files found in the mzml_files folder.")

  # 3. Load or compute RT tolerance
  rt_tol_path <- file.path(folder_path, "Sample-and-feature-wise-RT-tolerance.txt")
  if (file.exists(rt_tol_path)) {
    cat("\u2713 RT tolerance file exists! Using custom tolerance\n")
    RT_tolerance <- read.delim(rt_tol_path, check.names = FALSE)
    if (ncol(RT_tolerance) != 2 * length(fls)) {
      stop("RT tolerance file should contain RT_min and RT_max for all files.")
    }
    cat("\u2713 RT tolerance file format is valid\n")
  } else {
    cat(paste("RT tolerance file not found. Proceeding with default RT tolerance of ", tol_rt, " min and mz tolerance of ", tol_mz, " ppm.\n", sep = ""))
    RT_tolerance <- matrix(0, nrow = nrow(fil), ncol = 2 * length(fls))
    RT_tolerance[, 1:length(fls)] <- fil$RT - tol_rt
    RT_tolerance[, (length(fls)+1):(2*length(fls))] <- fil$RT + tol_rt
    clean_fls <- sub("\\.(mzML|mzXML)$", "", fls)
    colnames(RT_tolerance) <- c(paste("RT_min", clean_fls, sep = "_"), paste("RT_max", clean_fls, sep = "_"))
  }
  RT_tolerance <- RT_tolerance * 60  # Convert to seconds
  fil$mz_up <- fil$mz + (fil$mz * tol_mz) / 1e6
  fil$mz_down <- fil$mz - (fil$mz * tol_mz) / 1e6
  fil_rt <- cbind(fil, RT_tolerance)

  # 4. Read all spectra
  cat("Reading spectra files...\n")
  spectral_files <- vector("list", length(fls))
  flag <- character(length(fls))
  flag_count <- 1
  pb <- utils::txtProgressBar(min = 0, max = length(fls), style = 3)
  backend <- Spectra::MsBackendMzR()

  for (i in seq_along(fls)) {
    file_path <- file.path(mzml_path, fls[i])
    sps <- Spectra::Spectra(file_path, source = backend)
    sps <- Spectra::filterMsLevel(sps, msLevel = 2)
    spectral_files[[i]] <- if (length(sps) != 0) sps else NULL

    if (is.null(spectral_files[[i]]) || length(spectral_files[[i]]) == 0) {
      flag[flag_count] <- fls[i]
      flag_count <- flag_count + 1
    }
    utils::setTxtProgressBar(pb, i)
  }

  # 5. Post-processing
  close(pb)
  clean_fls <- sub("\\.(mzML|mzXML)$", "", fls)
  idx_ms2_absent <- which(sapply(spectral_files, function(x) is.null(x) || length(x) == 0))

  if (length(idx_ms2_absent) > 0) {
    cat("\nMS2 data not present in the following files:\n")
    print(clean_fls[idx_ms2_absent])

    pattern <- paste0("RT_(min|max)_(", paste(clean_fls[idx_ms2_absent], collapse = "|"), ")")
    cols_to_remove <- grepl(pattern, colnames(fil_rt))
    spectral_files_filtered <- spectral_files[-idx_ms2_absent]
    fil_rt_filtered <- fil_rt[, !cols_to_remove]
    fls_filtered <- fls[-idx_ms2_absent]
    names(spectral_files_filtered) <- fls_filtered
  } else {
    names(spectral_files) <- fls
  }

  cat("\n\u2713 Spectra have been extracted from the MS2 files\n")
  ro_path <- file.path(folder_path, "R_objects")
  if (!dir.exists(ro_path)) dir.create(ro_path)

  if (length(idx_ms2_absent) == 0) {
    for (l in seq_along(spectral_files)) {
      spectral_files[[l]]$spectrumId <- paste(names(spectral_files)[l], "_scan_", seq_along(spectral_files[[l]]), sep = "")
    }
    assign("spectral_files_f", spectral_files, envir = .dures_env)
    #assign("spectral_files_f", spectral_files, envir = .GlobalEnv)
    saveRDS(spectral_files, file.path(ro_path, "spectral_files.rds"))
    return(list(spectral_files = spectral_files, stats_file_analyzed = fil_rt, file_names = fls))
  } else {
    for (l in seq_along(spectral_files_filtered)) {
      spectral_files_filtered[[l]]$spectrumId <- paste(names(spectral_files_filtered)[l], "_scan_", seq_along(spectral_files_filtered[[l]]), sep = "")
    }
    assign("spectral_files_f", spectral_files_filtered, envir = .dures_env)
    #assign("spectral_files_f", spectral_files_filtered, envir = .dure)
    saveRDS(spectral_files_filtered, file.path(ro_path, "spectral_files.rds"))
    return(list(spectral_files = spectral_files_filtered, stats_file_analyzed = fil_rt_filtered, file_names = fls_filtered))
  }
}
