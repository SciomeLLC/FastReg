.onLoad <- function(libname, pkgname) {
  if (Sys.info()["sysname"] == "Windows") {
    # Windows may complain about mismatch HDF5 headers and library due to multiple installations of HDF5 library
    # assign("HDF5_DISABLE_VERSION_CHECK", 2, envir = .GlobalEnv)
    Sys.setenv(HDF5_DISABLE_VERSION_CHECK='2')
  }
}

.onUnload <- function(libname, pkgname) {
  if (Sys.info()["sysname"] == "Windows") {
    # rm("HDF5_DISABLE_VERSION_CHECK", envir = .GlobalEnv)
    Sys.unsetenv(HDF5_DISABLE_VERSION_CHECK)
  }
}