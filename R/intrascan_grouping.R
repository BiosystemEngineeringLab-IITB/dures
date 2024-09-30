#' intrascan_grouping
#'
#' A function that accepts a spectrum as input and returns a spectrum with fragments grouped within a given mass tolerance
#'
#' @param data A data frame of a spectrum containing two columns (fragments and intensity).
#' @param tol A tolerance value in daltons
#' @return A data frame with the same columns but mass fragments grouped.
#' @importFrom magrittr %>%
#' @importFrom dplyr %>% mutate

#' @export
intrascan_grouping <- function(scan, tol) {
  colnames(scan) = c("fragments", "intensity")
  if (length(scan[,1]) == 0) {
    return(NULL)
  } else {
    scan <- scan[order(scan[,1]), ]
    scan$unique_id <- paste("frag", scan[,1], sep = "_")

    scan <- scan %>%
      dplyr::arrange(fragments) %>%
      dplyr::mutate(group = cumsum(c(1, diff(fragments) >= tol)))

    agg_df <- aggregate(cbind(fragments = scan$fragments, intensity = scan$intensity),
                        by = list(group = scan$group),
                        FUN = function(x) c(mean = mean(x), sum = sum(x)))

    new_scan <- data.frame(fragments = agg_df[, 2][,1], intensity = agg_df[, 3][,2]) #print the agg_df dataframe to look for the second and fifth columns. Second column contain mean mz and fifth contain sum intensities
    new_scan <- new_scan[new_scan$intensity != 0, ]  # remove fragments with zero intensities

    return(new_scan)
  }
}
