#' tuning module
#'
#' Matches experimental and reference spectra using a predefined MS/MS fragment tolerance
#' @param l5 Contain the final set of features which matched with the reference at the MS1 level, output of precursor_matching.R
#' @param tolerance fragment ion tolerance at the MS/MS level, defaults to 0.05 Da
#' @param folder_path folder containing the input directory (mzml) feature list files where the results of the matching will be stored.
#' @return A dataframe with matching metrics, annotations and both experimental and reference spectrum identifiers. Features with matching score zero remain unannotated
#' @examples
#' # Example usage of the function
#' fragment_matching_before_denoising(folder_path, l4, tolerance = 0.05, "positive")
#' @export


tuning_module <- function(folder_path, l4, l5, l6, tolerance){

  freq_df_1 = list()
  for(i in 1:dim(l6)[1]){
    metabolite = l6$Feature_ID[i]
    sp = l4[[which(names(l4) %in% metabolite)]]
    sp = sp[grep(l6$Scan_Number[i], names(sp))]
    freq_df_1 = c(freq_df_1, sp)
  }
  thresholds=as.numeric(as.character(seq(0,1,0.01)))
  if (!dir.exists(file.path(folder_path, "Denoised_spectra_for_tuning"))) {
    dir.create(file.path(folder_path, "Denoised_spectra_for_tuning"))
  }

  names(freq_df_1) = paste(l6$Feature_ID, names(freq_df_1), sep="_")

  cat("Creating subspectra for every threshold...\n")
  pb <- txtProgressBar(min = 0, max = length(freq_df_1), style = 3)

  for (j in seq_along(freq_df_1)) {
    path <- file.path(folder_path, "Denoised_spectra_for_tuning")
    dir.create(file.path(path, l6$Feature_ID[j]), showWarnings = FALSE)

    for (i in seq_along(thresholds)) {
      sps_tic_20_all_freq <- tune_parameters_for_denoising(
        freq_df_1[[j]], thresholds[i], l6$Feature_ID[j], path
      )
    }

    gc()
    setTxtProgressBar(pb, j)
  }

  close(pb)



  output_dir <- file.path(folder_path, "After_denoising_matches")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ref_dir <- file.path(folder_path, "Reference")
  scan_dir <- file.path(folder_path, "Denoised_spectra_for_tuning")

  an <- l5

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

    #best_scan_bd = l6[which(l6$Feature_ID %in% metabolite), "Scan_Number"]
    best_ref_bd = l6[which(l6$Feature_ID %in% metabolite), "Splash_key"]

    #scan_n = scans[[grep(best_scan_bd, scan_names)]]
    ref_n = refs[grep(best_ref_bd, ref_names)]

    #scan_names_n = scan_names[[grep(best_scan_bd, scan_names)]]
    ref_names_n = ref_names[[grep(best_ref_bd, ref_names)]]

    dfs <- list()

    # Match scans and references
    for (i in seq_along(ref_n)) {
      for (j in seq_along(scans)) {
        file <- scans[[j]]
        sl <- paste(scan_names[j], ref_names_n[i], sep = "_")
        colnames(file) <- c("fragments", "intensity")

        sps_df <- calculate_SS_dataframe_no_thresholding(ref_n[[i]], file, tolerance = tolerance)
        sps_df <- sps_df %>%
          dplyr::mutate(
            Total_Library_fragments = nrow(ref_n[[i]]),
            Feature_ID = metabolite,
            Scan_Number = sub(".*(scan_.*)$", "\\1", scan_names[j]),
            Library_name = sub("^[0-9]+\\.[0-9]+_(.*?)_splash.*$", "\\1", ref_names_n[i]),
            Splash_key = sub(".*(splash.*)$", "\\1", ref_names_n[i]),
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

  for(i in 1:length(l6$Feature_ID)){
    dat = data.frame()
    met = l6$Feature_ID[i]
    d = read.csv(paste(folder_path, "After_denoising_matches/", met,".csv",sep=""))
    d$X = NULL
    signal_red = (l6$Number_of_Matching_Fragments[i] - d$Number_of_Matching_Fragments)/l6$Number_of_Matching_Fragments[i]
    percent_signal_red = signal_red * 100
    noise_frag_bd = l6$Total_Experimental_fragments[i] - l6$Number_of_Matching_Fragments[i]
    noise_frag_ad = d$Total_Experimental_fragments - d$Number_of_Matching_Fragments
    noise_red = (noise_frag_bd - noise_frag_ad)/noise_frag_bd
    percent_noise_red = noise_red * 100
    d$Signal_reduction = percent_signal_red
    d$Noise_reduction =percent_noise_red
    #d$SNR_bd = Top_raw_scan_level_matches_before_denoising_scores_above700_latest$NMF[i]/(Top_raw_scan_level_matches_before_denoising_scores_above700_latest$`Tot. Exp fragments`[i] - Top_raw_scan_level_matches_before_denoising_scores_above700_latest$NMF[i])
    #d$SNR_ad = d$NMF/(d$Tot..Exp.fragments - d$NMF)
    d$freq = as.numeric(unlist(lapply(d$ID, function(x) strsplit(x, "_")[[1]][2])))
    #d$SignalbyNoisereduction = d$Signal_reduction/d$Noise_reduction
    #d$SNR_difference = d$SNR_ad - d$SNR_bd
    d$FMR_ad = d$Number_of_Matching_Fragments/d$Total_Experimental_fragments
    d$FMR_bd = l6$Number_of_Matching_Fragments[i]/l6$Total_Experimental_fragments[i]
    d$Feature_ID = met
    #zeros_df <- data.frame(matrix(0, nrow = (101-length(which(d$ID!=0))), ncol = ncol(d)))
    #colnames(zeros_df) = colnames(d)
    #zeros_df$freq = setdiff(thresholds, d$freq)
    #zeros_df$Feature_ID = met
    #dat = rbind(d, )
    write.csv(d, paste(folder_path, "After_denoising_matches/",met,".csv",sep=""), quote = F)

  }

  #pareto optimization
  cat("Starting pareto optimization..")

  if (!dir.exists(file.path(folder_path, "pareto_results"))) {
    dir.create(file.path(folder_path, "pareto_results"))
  }
  dir.create(file.path(folder_path, "pareto_results", "pdf"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(folder_path, "pareto_results", "csv"), recursive = TRUE, showWarnings = FALSE)

  # Initialize data frames to store results
  freq <- data.frame()
  id <- list()
  k <- 1

  # Optimization Function to Minimize Signal Loss
  optimize_weights <- function(weights, df) {
    weight_similarity <- weights[1]
    weight_noise <- weights[2]
    weight_signal <- weights[3]

    df$Optimized_Score <- (weight_similarity * df$Similarity_Score) +
      (weight_noise * df$Noise_reduction) -
      (weight_signal * df$Signal_reduction)

    selected_features <- df[which.max(df$Optimized_Score), ]
    return(mean(selected_features$Signal_reduction))
  }

  # Define bounds for optimization (Weights must sum to 1)
  weight_constraints <- function(weights) {
    return(sum(weights) - 1)
  }

  # Progress bar setup
  pb <- txtProgressBar(min = 0, max = length(l6$ID), style = 3)

  # Loop through all features
  for (i in seq_along(l6$ID)) {
    setTxtProgressBar(pb, i)
    met <- l6$Feature_ID[i]

    # Read the CSV file containing post-denoising solutions
    test <- read.csv(file.path(folder_path, "After_denoising_matches", paste0(met, ".csv")))
    test$X <- NULL
    test$Noise_reduction <- 100 - test$Noise_reduction

    if (any(is.na(test$Signal_reduction)) || any(is.na(test$Noise_reduction))) {
      id[[k]] <- met
      k <- k + 1
      next
    }

    sky1 <- rPref::psel(test, rPref::low(Signal_reduction) * rPref::low(Noise_reduction))
    sky1$metabolite <- met

    if (nrow(sky1) >= 3) {
      sky1_sorted <- sky1[!duplicated(sky1$Signal_reduction), ]
      sky1_sorted <- sky1_sorted[order(sky1_sorted$Signal_reduction), ]

      if (nrow(sky1_sorted) >= 3) {
        x <- sky1_sorted$Signal_reduction
        y <- sky1_sorted$Noise_reduction

        dx <- diff(x)
        dy <- diff(y)
        dx[dx == 0] <- 1e-6

        first_derivative <- dy / dx
        d2x <- diff(x[-length(x)])
        d2x[d2x == 0] <- 1e-6

        second_derivative <- diff(first_derivative) / d2x

        if (all(is.finite(second_derivative))) {
          knee_index <- which.max(abs(second_derivative)) + 1
          final_solution <- sky1_sorted[knee_index, ]
        } else {
          final_solution <- sky1_sorted[which.min(sky1_sorted$Signal_reduction), ]
        }
      } else {
        final_solution <- sky1_sorted[which.min(sky1_sorted$Signal_reduction), ]
      }
    } else {
      opt_result <- DEoptim::DEoptim(optimize_weights, lower = c(0, 0, 0), upper = c(1, 1, 1),
                                     control = list(itermax = 100), df = sky1)
      best_weights <- opt_result$optim$bestmem

      sky1$Weighted_Score <- (best_weights[1] * sky1$Similarity_Score) +
        (best_weights[2] * sky1$Noise_reduction) -
        (best_weights[3] * sky1$Signal_reduction)

      final_solution <- sky1[which.max(sky1$Weighted_Score), ]
    }

    freq <- rbind(freq, final_solution)

    # Save Pareto front
    write.csv(sky1, file.path(folder_path, "pareto_results", "csv", paste0(met, ".csv")), quote = FALSE)

    # Save plot
    g <- ggplot(test, aes(x = Signal_reduction, y = Noise_reduction, label = freq)) +
      labs(y = "100-Noise_reduction") +
      geom_point(shape = 21) +
      geom_point(data = sky1, size = 3, col = "red") +
      geom_text(hjust = 0, vjust = 0) +
      ggtitle(met)

    ggsave(file.path(folder_path, "pareto_results", "pdf", paste0(met, ".pdf")), plot = g)
  }

  close(pb)


  final_freq <- merge(l6, freq,  by="Feature_ID")

  #final_freq$percentage_increase_in_SS = 100 * (final_freq$Matching_Score.y - final_freq$Matching_Score.x)/final_freq$

  cat("summary statistics of optimal frequencies..\n")

  print(summary(final_freq$freq))


  return(final_freq)

}


# Function to process each file
process_file <- function(i, output_dir, ids) {
  file_path <- file.path(output_dir, paste0(ids[i], ".csv"))
  if (file.exists(file_path)) {
    f <- data.table::fread(file_path)
    setorder(f, -Matching_Score)  # Sort by Similarity_Score in descending order
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






tune_parameters_for_denoising <- function(agg_df, threshold, name, path){
  colnames(agg_df)[2:4] = c("Mean_MZ", "Mean_Intensity", "Frequency")
  agg_df_cutoff = subset(agg_df, agg_df$Frequency >= threshold)
  agg_df_cutoff = agg_df_cutoff[order(agg_df_cutoff$Mean_MZ, decreasing = F),]
  #creating the composite spectra
  spd <- DataFrame(
    msLevel = c(2L),
    polarity = c(1L),
    id = c(name),
    name = c(name))
  spd$mz <- list(
    c(agg_df_cutoff$Mean_MZ))
  spd$intensity <- list(
    c(agg_df_cutoff$Mean_Intensity))
  sps <- Spectra(spd)
  #plotSpectra(sps, main = paste(sps$name, " (Freq_cutoff_", as.character(threshold),")", sep=""))
  # folder_create = dir.create(paste("~/metabolomics/WTC+DIELdata/results_ppm=5/Denoised_spectra/",name,"/", sep=""))
  #setwd(paste("~/metabolomics/Sportomics/Database_metabolites/Aggregate_top_5_files_SS_score/Files_of_interest_Sportomic/After_filtering_stage_2_spectra/", name,sep=""))
  name_exp= paste(path,"/", name,"/", "denoised_",threshold, ".txt", sep="")
  sps_df=data.frame(fragments = peaksData(sps)[[1]][,1], intensity =  peaksData(sps)[[1]][,2])
  write.table(sps_df, name_exp, sep="\t", col.names = T, row.names = F, quote = F)
  return(sps)
}
