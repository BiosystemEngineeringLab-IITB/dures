
# Prepare individual references
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
