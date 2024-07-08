.dures_env <- new.env(parent = emptyenv())
utils::globalVariables(c("spectral_files_f"))
# Define the .onLoad function to run any additional setup code
.onLoad <- function(libname, pkgname) {
  .dures_env <<- new.env(parent = emptyenv())
  # This function is called when the package is loaded
  # Initialize any other package-wide settings here if needed
}
