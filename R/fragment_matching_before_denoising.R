#' fragment_matching_before_denoising
#'
#' Matches experimental and reference spectra using a predefined MS/MS fragment tolerance
#' @param l4 Contain the final set of features which matched with the reference at the MS1 level
#' @param tolerance fragment ion tolerance at the MS/MS level, defaults to 0.05 Da
#' @param folder_path folder containing the input directory (mzml) feature list files where the results of the matching will be stored.
#' @param ionization_mode positive or negative modes
#' @return A dataframe with matching metrics, annotations and both experimental and reference spectrum identifiers. Features with matching score zero remain unannotated
#' @examples
#' # Example usage of the function
#' fragment_matching_before_denoising(folder_path, l4, tolerance = 0.05, "positive")
#' @export


fragment_matching_before_denoising <- function(folder_path, l4, tolerance, ionization_mode){

  output_dir <- file.path(folder_path, "Before_denoising_matches")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ref_dir <- file.path(folder_path, "Reference")
  scan_dir <- file.path(folder_path, "MS2_scans_before_denoising")

  an <- l4

  # Iterate through all metabolites
  for (k in seq_len(nrow(an))) {
    metabolite <- an$ID[k]
    mz <- an$mz[k]
    message("Processing Metabolite: ", metabolite, " Entry: ", k)

    scan_path <- file.path(scan_dir, metabolite)
    ref_path <- file.path(ref_dir, metabolite)

    scans <- lapply(list.files(scan_path, pattern = "\\.txt$", full.names = TRUE), data.table::fread)
    scan_names <- gsub("\\.txt$", "", list.files(scan_path, pattern = "\\.txt$"))

    refs <- lapply(list.files(ref_path, pattern = "\\.txt$", full.names = TRUE), data.table::fread)
    ref_names <- list.files(ref_path, pattern = "\\.txt$")

    dfs <- list()

    # Match scans and references
    for (i in seq_along(refs)) {
      for (j in seq_along(scans)) {
        file <- scans[[j]]
        sl <- paste(scan_names[j], ref_names[i], sep = "_")
        colnames(file) <- c("fragments", "intensity")

        sps_df <- calculate_SS_dataframe_no_thresholding(refs[[i]], file, tolerance = tolerance)
        sps_df <- sps_df %>%
          dplyr::mutate(
            Total_Library_fragments = nrow(refs[[i]]),
            Feature_ID = metabolite,
            Scan_Number = sub(".*(scan_.*)$", "\\1", scan_names[j]),
            Library_name = sub("^[0-9]+\\.[0-9]+_(.*?)_splash.*$", "\\1", ref_names[i]),
            Splash_key = sub(".*(splash.*)$", "\\1", ref_names[i]),
            ID = sl,
            mz = mz
          ) %>%
          dplyr::arrange(desc(Modified_Dot_Product))

        dfs[[length(dfs) + 1]] <- sps_df
      }
    }

    # Combine results and save
    final_df <- dplyr::bind_rows(dfs) %>%
      dplyr::arrange(desc(Modified_Dot_Product))

    final_df$Matched_Ref_Fragments <- NULL
    final_df$Matched_Exp_Fragments <- NULL

    data.table::fwrite(final_df, file.path(output_dir, paste0(metabolite, ".csv")))
  }

  ids <- an$ID

  # Ensure process_file is called correctly with IDs as argument
  top_matches_list <- lapply(seq_along(ids), function(i) process_file(i, output_dir, ids))

  # Remove NULL entries
  top_matches_list <- Filter(Negate(is.null), top_matches_list)

  if(length(top_matches_list) == 0) {
    stop("No valid matching results found; check input files or paths.")
  }

  # Combine using dplyr
  Top_raw_scan_level_matches_before_denoising <- dplyr::bind_rows(top_matches_list)

  # Verify object explicitly
  if (!is.data.frame(Top_raw_scan_level_matches_before_denoising)) {
    stop("Combined object is not a valid data.frame. Check your input files.")
  }

  # Ensure 'ID' column exists
  if (!"ID" %in% colnames(Top_raw_scan_level_matches_before_denoising)) {
    stop("Column 'ID' missing. Check process_file function.")
  }

  # Add Feature_ID and mz safely
  #Top_raw_scan_level_matches_before_denoising$Feature_ID <- ids[match(Top_raw_scan_level_matches_before_denoising$ID, paste0(ids, ".csv"))]
  #Top_raw_scan_level_matches_before_denoising$mz <- an$mz[match(Top_raw_scan_level_matches_before_denoising$Feature_ID, an$ID)]

  # Remove NA Matching_Score rows safely
  Top_raw_scan_level_matches_before_denoising <- Top_raw_scan_level_matches_before_denoising[!is.na(Top_raw_scan_level_matches_before_denoising$Matching_Score), ]

  # Select appropriate library based on ionization mode
  lib <- if (ionization_mode == "positive") .GlobalEnv$.positive_lib else .GlobalEnv$.negative_lib

  # Initialize Metabolite_Name with NA
  Top_raw_scan_level_matches_before_denoising$Metabolite_Name <- NA_character_

  # Identify rows with Matching_Score > 0
  positive_score_idx <- which(Top_raw_scan_level_matches_before_denoising$Matching_Score > 0)

  # Remove .txt extension from Splash_key
  Top_raw_scan_level_matches_before_denoising$Splash_key <- sub("\\.txt$", "", Top_raw_scan_level_matches_before_denoising$Splash_key)

  # Map metabolite names for positive Matching_Score rows
  matched_names <- lib$name[match(
    Top_raw_scan_level_matches_before_denoising$Splash_key[positive_score_idx],
    lib$splash_keys
  )]

  # Assign these matched names back to the appropriate rows
  Top_raw_scan_level_matches_before_denoising$Metabolite_Name[positive_score_idx] <- matched_names

  return(Top_raw_scan_level_matches_before_denoising)

}


# Function to process each file and return the top match in terms of matching score
process_file <- function(i, output_dir, ids) {
  file_path <- file.path(output_dir, paste0(ids[i], ".csv"))
  if (file.exists(file_path)) {
    f <- data.table::fread(file_path)
    data.table::setorder(f, -Matching_Score)  # Sort by Similarity_Score in descending order
    return(f[1, ])  # Return the first row
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
}

# Function to get the nearest index
getNearestIdx <- function(array, value) {
  idx <- findInterval(value, array)
  if (idx > 0 && (idx == length(array) || abs(value - array[idx]) < abs(value - array[idx + 1]))) {
    return(idx)
  } else {
    return(idx + 1)
  }
}

# Function to calculate forward and reverse dot products with matched fragments
calculateFDotnRDot <- function(ref_spectra, exp_spectra, tol) {
  mz_r <- ref_spectra$fragments
  mz_e <- exp_spectra$fragments
  i_r <- sqrt(ref_spectra$intensity / max(ref_spectra$intensity, na.rm = TRUE))
  i_e <- sqrt(exp_spectra$intensity / max(exp_spectra$intensity, na.rm = TRUE))

  w_r <- (sum(i_r) - 0.5) / (sum(i_r) - 0.5 + i_r)
  w_e <- (sum(i_e) - 0.5) / (sum(i_e) - 0.5 + i_e)

  mask <- (mz_r > (min(mz_e) - tol)) & (mz_r < (max(mz_e) + tol))
  mz_r <- mz_r[mask]
  i_r <- i_r[mask]
  w_r <- w_r[mask]

  wi_r <- w_r * i_r
  wi_e <- w_e * i_e

  num <- 0
  denr_e <- 0
  denf_r <- 0
  cnt <- 0

  # Initialize lists to store matched fragments
  matched_ref_fragments <- c()
  matched_exp_fragments <- c()

  for (i in seq_along(mz_e)) {
    idx <- getNearestIdx(mz_r, mz_e[i])

    # Check if idx or mz_r[idx] is NA
    if (!is.na(idx) && idx > 0 && !is.na(mz_r[idx]) && abs(mz_r[idx] - mz_e[i]) <= tol) {
      idx2 <- getNearestIdx(mz_e, mz_r[idx])

      if (!is.na(idx2) && idx2 == i) {
        num <- num + wi_r[idx] * wi_e[i]
        denr_e <- denr_e + wi_e[i]^2
        denf_r <- denf_r + wi_r[idx]^2
        cnt <- cnt + 1

        # Store the matched fragments
        matched_ref_fragments <- c(matched_ref_fragments, mz_r[idx])
        matched_exp_fragments <- c(matched_exp_fragments, mz_e[i])
      }
    }
  }

  if (denr_e == 0) {
    return(list(0, 0, 0, 0, nrow(exp_spectra), matched_ref_fragments, matched_exp_fragments))
  }

  denr_r <- sum(wi_r^2)
  denf_e <- sum(wi_e^2)
  num <- num^2
  forward_dot <- num / (denf_e * denf_r)
  reverse_dot <- num / (denr_e * denr_r)
  mdp <- (forward_dot + reverse_dot) * 0.5

  return(list(mdp, forward_dot, reverse_dot, cnt, nrow(exp_spectra), matched_ref_fragments, matched_exp_fragments))
}

# Function to calculate similarity score including matched fragments
calculateSimilarityScore <- function(ref_spectra, exp_spectra, tol) {
  res <- calculateFDotnRDot(ref_spectra, exp_spectra, tol)
  modified_dot_product_score <- res[[1]]
  forward_dot <- res[[2]]
  reverse_dot <- res[[3]]
  number_of_fragments_matched <- res[[4]]
  tot_exp_fragments <- res[[5]]

  # Add matched fragments to the output
  matched_ref_fragments <- res[[6]]
  matched_exp_fragments <- res[[7]]

  matching_score <- 5 * (modified_dot_product_score * 100 + 20.0 * log2(max(number_of_fragments_matched, 1)))
  return(list(matching_score, modified_dot_product_score, number_of_fragments_matched,
              forward_dot, reverse_dot, tot_exp_fragments,
              matched_ref_fragments, matched_exp_fragments))
}

# Function to calculate SS dataframe without thresholding and include matched fragments
calculate_SS_dataframe_no_thresholding <- function(ref_spectra, file, tolerance) {
  scores <- calculateSimilarityScore(ref_spectra, file, tolerance)

  sps_df <- data.frame(
    Matching_Score = scores[[1]],
    Modified_Dot_Product = scores[[2]],
    Number_of_Matching_Fragments = scores[[3]],
    Forward_Dot_Product = scores[[4]],
    Reverse_Dot_Product = scores[[5]],
    Total_Experimental_fragments = scores[[6]],
    Matched_Ref_Fragments = I(list(scores[[7]])), # Store as list to handle multiple values
    Matched_Exp_Fragments = I(list(scores[[8]]))
  )

  return(sps_df[1, , drop = FALSE])
}

