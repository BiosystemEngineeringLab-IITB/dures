#' prepare_reference
#'
#' Stores the reference spectra as a txt file for further matching with current experimental spectra and returns the reference spectra as a Spectra object. This neednot be run independently. This is a companion function to precursor_matching.R
#'
#' @param index index
#' @param name feature ID
#' @param  lib positive or negative library
#' @param path folder_path to store the spectra as txt files
#' @return reference spectra in Spectra format
#' @examples
#' # Example usage of the function
#' prepare_reference(ind, n, lib, p)
#' @export

prepare_reference <- function(index, name, lib, path) {
  # Extract spectra data
  ref <- data.frame(
    mz = unlist(lib$ms2f[index]),
    intensity = unlist(lib$ms2i[index])
  )
  ref <- ref[order(ref$mz), ]  # Sort by m/z

  # Create Spectra object
  spd_ref <- DataFrame(
    msLevel = 2L,
    polarity = 1L,
    id = index,
    name = name,
    mass = list(ref$mz),
    inten = list(ref$intensity)
  )
  sps_ref <- Spectra(spd_ref)

  # Save the spectra data
  sps_ref_df <- data.frame(
    fragments =ref$mz,
    intensity = ref$intensity
  )
  name_ref <- file.path(path, paste0(name, ".txt"))
  write.table(sps_ref_df, name_ref, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  return(sps_ref)
}
