#' sensitivity analysis
#'
#'
#' @param l7 Contain the final set of features with matching metrics before and after denoising
#' @param folder_path where the user wants to store the sensitivity plots
#' @return A dataframe with significance values for pairwise comparisons in matching metrics across different thresholds
#' @examples
#' # Example usage of the function
#' sensitivity_analysis(folder_path, l7)
#' @export

sensitivity_analysis <- function(l7, folder_path){

  cat("Checking if significant improvement in matching score was observed after denoising for atleast one feature...\n")
  combined_freq_selected <- l7[which(l7$percentage_increase_in_SS > 0), ]
  if(dim(combined_freq_selected)[1] == 0){
    cat("After denoising no single feature showed improvement in matching score. Hence proceeding with the entire set of features..\n")
    combined_freq_selected = l7
  }else{
    print(paste0("Improvement noted in ", as.character(dim(combined_freq_selected)[1]), " number of features out of ", as.character(dim(l7)[1])))
  }

  cat("Starting Sensitivity analysis..\n")
  path = paste0(folder_path, "After_denoising_matches/")
  dat_extra_cols = data.frame()
  for(i in 1:length(combined_freq_selected$Feature_ID)){
    met = combined_freq_selected$Feature_ID[i]
    dat = read.csv(paste(path,"/", met, ".csv", sep=""))
    dat$before_denoising_score = combined_freq_selected$Similarity_Score.x[i]
    dat_extra_cols = rbind(dat_extra_cols, dat)
  }

  filtered_data = dat_extra_cols

  # Count occurrences of each frequency value
  freq_counts <- filtered_data %>%
    dplyr::group_by(freq) %>%
    dplyr::summarise(freq_count = dplyr::n())

  # Compute scaling factor: count per frequency divided by total unique Feature_IDs
  total_unique_features <- dplyr::n_distinct(filtered_data$Feature_ID)
  freq_counts <- freq_counts %>%
    dplyr::mutate(scaling_factor = freq_count / total_unique_features)

  # Merge scaling factor back into the original dataset
  filtered_data <- dplyr::left_join(filtered_data, freq_counts, by = "freq")
  #View(filtered_data)

  filtered_data$Matching_Score_new = filtered_data$Matching_Score * filtered_data$scaling_factor
  filtered_data$noise_fragments = filtered_data$Total_Experimental_fragments - filtered_data$Number_of_Matching_Fragments
  filtered_data$noise_fragments = filtered_data$noise_fragments * filtered_data$scaling_factor

  df_filtered = filtered_data

  # Step 1: Compute median values for Similarity Score and MNR across frequency thresholds
  median_similarity <- df_filtered %>% dplyr::group_by(freq) %>% dplyr::summarize(Median_matching_score = median(Matching_Score_new))
  median_nf <- df_filtered %>% dplyr::group_by(freq) %>% dplyr::summarize(Median_nf = median(noise_fragments))
  # Calculate Percent Reduction in Noise Fragments relative to freq = 0
  median_nf <- median_nf %>%
    mutate(Percent_Reduction_Noise = (Median_nf - Median_nf[freq == 0]) / Median_nf[freq == 0] * 100)
  # Step 2: Identify frequency where median is maximum
  max_median_similarity_freq <- median_similarity$freq[which.max(median_similarity$Median_matching_score)]
  max_median_nf_freq <- median_nf$freq[which.max(median_nf$Median_nf)]

  #max_median_nf_freq

  # Step 3: Compute feature count per frequency
  feature_counts <- df_filtered %>% dplyr::group_by(freq) %>% dplyr::summarize(Feature_Count = dplyr::n_distinct(Feature_ID))

  # Step 4: Define the feature count threshold (80% of max feature count)
  max_feature_count <- max(feature_counts$Feature_Count)
  feature_count_threshold <- 0.8 * max_feature_count

  p_values_similarity <- compute_p_values(df_filtered, "Matching_Score_new", max_median_similarity_freq)
  #p_values_mnr <- compute_p_values(df_filtered, "MNR_new", max_median_mnr_freq)
  #p_value_texp <- compute_p_values(df_filtered, "texp_new", max_median_texp_freq)
  p_value_nf <- compute_p_values(df_filtered, "noise_fragments", max_median_nf_freq)

  median_nf = median_nf %>%
    dplyr::left_join(p_value_nf, by="freq")

  # Step 6: Merge feature count and p-values
  # verification_table <- feature_counts %>%
  #   left_join(p_values_similarity, by = "freq") %>%
  #   rename(P_Value_Similarity = P_Value) %>%
  #   left_join(p_values_mnr, by = "freq") %>%
  #   rename(P_Value_MNR = P_Value) %>%
  #   left_join(p_value_nf, by = "freq") %>%
  #   rename(P_value_noise_frag = P_Value)

  # Step 7: Function to find stable range based on stopping criteria
  # Step 8: Find stable ranges for Similarity Score and MNR
  stable_range_similarity <- find_stable_range(median_similarity, max_median_similarity_freq, feature_counts, p_values_similarity, feature_count_threshold = feature_count_threshold)
  #stable_range_mnr <- find_stable_range(median_mnr, max_median_mnr_freq, feature_counts, p_values_mnr)
  stable_range_nf <- find_stable_range(median_nf, max_median_nf_freq, feature_counts, p_value_nf, feature_count_threshold)

  g1 = # Create the plot
    ggplot(median_similarity, aes(x = freq, y = Median_matching_score)) +
    geom_line(color = "blue") +
    geom_point() +
    geom_vline(xintercept = max_median_similarity_freq, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = seq(min(median_similarity$freq), max(median_similarity$freq), by = 0.1)) +
    annotate("rect", xmin = stable_range_similarity[1], xmax = stable_range_similarity[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightgreen") +
    annotate("text", x = mean(stable_range_similarity), y = max(median_similarity$Median_matching_score),
             label = paste0("Stable Range: [", round(stable_range_similarity[1], 2), ", ", round(stable_range_similarity[2], 2), "]"),
             size = 5, hjust = 0.5, vjust = -1, color = "black", fontface = "bold") +
    labs(
      title = "Sensitivity Analysis of Similarity Score",
      x = "Recurrence Frequency",
      y = "Median Similarity Score"
    ) +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14)
    )

  # First, create the dataframe explicitly
  significant_points <- median_nf %>%
    dplyr::filter(!is.na(P_Value), P_Value < 0.05)

  g2 <- ggplot(median_nf, aes(x = freq, y = Median_nf)) +
    geom_line(color = "blue") +
    geom_point() +
    geom_vline(xintercept = max_median_nf_freq, linetype = "dashed", color = "red") +
    annotate("rect", xmin = stable_range_nf[1], xmax = stable_range_nf[2],
             ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightgreen") +

    # Add significance stars only if significant points exist
    {if(nrow(significant_points) > 0) geom_text(data = significant_points,
                                                aes(x = freq, y = Median_nf, label = "*"),
                                                color = "red", size = 5, vjust = -0.5)} +

    scale_x_continuous(breaks = seq(min(median_nf$freq), max(median_nf$freq), by = 0.1)) +
    labs(
      title = "Sensitivity Analysis of Noise Fragments",
      x = "Recurrence Frequency",
      y = "Median Noise Fragments"
    ) +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14)
    )
  print(g1 / g2)

  ggsave(
    filename = paste0(folder_path, "/Sensitivity_Analysis_Matching_Score.pdf"),
    plot = g1,
    device = "pdf",
    width = 8, height = 6, units = "in"
  )

  ggsave(
    filename = paste0(folder_path, "/Sensitivity_Analysis_Noise_Fragments.pdf"),
    plot = g2,
    device = "pdf",
    width = 8, height = 6, units = "in"
  )



}



find_stable_range <- function(median_df, max_freq, feature_counts, p_values, feature_count_threshold) {
  stable_freqs <- median_df$freq
  lower_bound <- max_freq
  upper_bound <- max_freq

  for (freq in rev(stable_freqs)) {
    if (freq < max_freq) {
      p_val_current <- p_values$P_Value[p_values$freq == freq]
      feature_count_current <- feature_counts$Feature_Count[feature_counts$freq == freq]

      condition_p_value <- ifelse(is.na(p_val_current), FALSE, p_val_current < 0.05)
      condition_feature_count <- ifelse(is.na(feature_count_current), FALSE, feature_count_current < feature_count_threshold)

      if (condition_p_value || condition_feature_count) {
        break
      }
      lower_bound <- freq
    }
  }

  for (freq in stable_freqs) {
    if (freq > max_freq) {
      p_val_current <- p_values$P_Value[p_values$freq == freq]
      feature_count_current <- feature_counts$Feature_Count[feature_counts$freq == freq]

      condition_p_value <- ifelse(is.na(p_val_current), FALSE, p_val_current < 0.05)
      condition_feature_count <- ifelse(is.na(feature_count_current), FALSE, feature_count_current < feature_count_threshold)

      if (condition_p_value || condition_feature_count) {
        break
      }
      upper_bound <- freq
    }
  }

  return(c(lower_bound, upper_bound))
}


# Step 5: Compute Wilcoxon p-values
compute_p_values <- function(df, metric, max_freq) {
  p_values <- sapply(unique(df$freq), function(freq) {
    if (freq == max_freq) return(NA)

    group1 <- df[df$freq == freq, metric]
    group2 <- df[df$freq == max_freq, metric]

    if(length(unique(group1)) <= 1 && length(unique(group2)) <= 1){
      return(NA) # Wilcox test won't run properly if there isn't enough variability
    }

    test_result <- wilcox.test(group1, group2, exact = FALSE)
    return(test_result$p.value)
  })
  return(data.frame(freq = unique(df$freq), P_Value = p_values))
}
