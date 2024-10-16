#include <blas_library_manager.h>

//////////////////////////////////////////////////
// @brief detects which BLAS library is linked by checking for specific function
// symbols
// @return  none
//////////////////////////////////////////////////
void BLASLibraryManager::detect_lib() {
  blas_handle = load_blas_library();
  if (!blas_handle) {
    Rcpp::Rcout << "Failed to load the BLAS library!" << std::endl;
    return;
  }

  // Check for OpenBLAS
  openblas_get_num_threads =
      (openblas_get_num_threads_func)resolve_symbol("openblas_get_num_threads");
  openblas_set_num_threads =
      (openblas_set_num_threads_func)resolve_symbol("openblas_set_num_threads");
  if (openblas_get_num_threads && openblas_set_num_threads) {
    blas_type = "OpenBLAS";
    Rcpp::Rcout << "Detected OpenBLAS/GotoBLAS library." << std::endl;
    original_num_threads = openblas_get_num_threads();
    return;
  }

  // Check for Intel MKL
  mkl_get_max_threads =
      (mkl_get_max_threads_func)resolve_symbol("MKL_Get_Max_Threads");
  mkl_set_num_threads =
      (mkl_set_num_threads_func)resolve_symbol("MKL_Set_Num_Threads");
  if (mkl_get_max_threads && mkl_set_num_threads) {
    blas_type = "Intel MKL";

    Rcpp::Rcout << "Detected Intel MKL library." << std::endl;
    original_num_threads = mkl_get_max_threads();
    return;
  }

  // Check for BLIS
  blis_get_num_threads =
      (blis_get_num_threads_func)resolve_symbol("bli_thread_get_num_threads");
  blis_set_num_threads =
      (blis_set_num_threads_func)resolve_symbol("bli_thread_set_num_threads");
  if (blis_get_num_threads && blis_set_num_threads) {
    blas_type = "BLIS";
    Rcpp::Rcout << "Detected BLIS library." << std::endl;
    original_num_threads = blis_get_num_threads();
    return;
  }

#ifdef __APPLE__
  // If on macOS, assume Apple's Accelerate
  blas_type = "Apple Accelerate";
  Rcpp::Rcout << "Detected Apple Accelerate library.\n" << std::endl;
  // Apple Accelerate does not provide thread management API, assume 1 thread
  original_num_threads = 1;
#endif
  original_num_threads = 1;
  Rcpp::Rcout << "BLAS library is unknown or does not support thread detection."
              << std::endl;
}

//////////////////////////////////////////////////
// @brief gets the number of BLAS threads by calling the underlying BLAS
// libraries function
// @return integer count of BLAS threads
//////////////////////////////////////////////////
int BLASLibraryManager::get_num_threads() const {
  if (blas_type == "OpenBLAS" && openblas_get_num_threads) {
    return openblas_get_num_threads();
  }
  if (blas_type == "Intel MKL" && mkl_get_max_threads) {
    return mkl_get_max_threads();
  }

  if (blas_type == "BLIS" && blis_get_num_threads) {
    return blis_get_num_threads();
  }
  return original_num_threads; // Default to 1
}
//////////////////////////////////////////////////
// @brief sets the number of BLAS threads by calling the underlying BLAS
// libraries function
// @return none
//////////////////////////////////////////////////
void BLASLibraryManager::set_num_threads(int num_threads) {
  if (num_threads <= 0) {
    return;
  }
  if (blas_type == "OpenBLAS" && openblas_set_num_threads) {
    openblas_set_num_threads(num_threads);
    // Rcpp::Rcout << "OpenBLAS threads set to " << num_threads << "."
    //             << std::endl;
  } else if (blas_type == "Intel MKL" && mkl_set_num_threads) {
    mkl_set_num_threads(num_threads);
    // Rcpp::Rcout << "MKL threads set to " << num_threads << "."
    //             << std::endl;
  } else if (blas_type == "BLIS" && blis_set_num_threads) {
    blis_set_num_threads(num_threads);
    // Rcpp::Rcout << "BLIS threads set to" << num_threads << "."
    //             << std::endl;
  } else {
    Rcpp::Rcout << "Cannot set threads: BLAS library is either unknown or does "
                   "not support this operation."
                << std::endl;
  }
}

//////////////////////////////////////////////////
// @brief tries to load the BLAS dll / shared lib
// @return pointer to BLAS library
//////////////////////////////////////////////////
void *BLASLibraryManager::load_blas_library() {
#ifdef WIN32
  return LoadLibrary("Rblas.dll");
#else
  return dlopen(NULL, RTLD_NOW | RTLD_GLOBAL); // Load current process libraries
#endif
}

void *BLASLibraryManager::resolve_symbol(const char *symbol) {
#ifdef WIN32
  return (void *)GetProcAddress((HMODULE)blas_handle, symbol);
#else
  return dlsym(blas_handle, symbol);
#endif
}