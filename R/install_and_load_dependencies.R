#' Install and Load All Required Packages for 'dures'
#'
#' Installs missing packages from CRAN or Bioconductor and loads all required packages.
#' @export
install_and_load_dures_dependencies <- function() {
  # All required CRAN and Bioc packages
  cran_pkgs <- c("dplyr", "readr", "data.table", "pbapply", "magrittr",
                 "utils", "stats", "rPref", "ggplot2", "DEoptim", "patchwork")
  bioc_pkgs <- c("S4Vectors", "Spectra")
  suggested_pkgs <- c("BiocManager", "knitr", "markdown")

  all_pkgs <- c(cran_pkgs, bioc_pkgs, suggested_pkgs)

  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager from CRAN...")
    install.packages("BiocManager")
  }

  # Install and load each package
  for (pkg in all_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% bioc_pkgs) {
        message(sprintf("Installing Bioconductor package: %s", pkg))
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        message(sprintf("Installing CRAN package: %s", pkg))
        install.packages(pkg)
      }
    } else {
      message(sprintf("✓ %s is already installed.", pkg))
    }

    # Try loading the package
    success <- suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
    if (success) {
      message(sprintf("✓ %s loaded successfully.", pkg))
    } else {
      warning(sprintf("⚠️ Failed to load %s after installation.", pkg))
    }
  }

  message("✅ All dependencies checked and loaded.")
}
