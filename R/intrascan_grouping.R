#' Example function using the pipe operator
#'
#' @param data A data frame.
#' @return A data frame with the first 10 rows.
#' @importFrom magrittr %>%
#' @importFrom dplyr %>% mutate

#' @export
intrascan_grouping <- function(scan, tol) {
  if (length(scan$fragments) == 0) {
    return(NULL)
  } else {
    scan <- scan[order(scan$fragments, decreasing = TRUE), ]
    scan$unique_id <- paste("frag", scan$fragments, sep = "_")

    scan <- scan %>%
      dplyr::arrange(fragments) %>%
      dplyr::mutate(group = cumsum(c(1, diff(fragments) > tol)))

    agg_df <- aggregate(cbind(fragments = scan$fragments, intensity = scan$intensity),
                        by = list(group = scan$group),
                        FUN = function(x) c(mean = mean(x), sum = sum(x)))

    new_scan <- data.frame(fragments = agg_df$fragments[, "mean"], intensity = agg_df$intensity[, "sum"])
    new_scan <- new_scan[new_scan$intensity != 0, ]  # remove fragments with zero intensities

    return(new_scan)
  }
}
