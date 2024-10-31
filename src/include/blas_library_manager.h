#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <RcppArmadillo.h>
#ifdef WIN32
#include <windows.h>
#else
#include <dlfcn.h> // For POSIX systems (Linux/macOS)
#endif

// BLAS function pointers for MKL/OpenBLAS
typedef void (*openblas_set_num_threads_func)(int);
typedef int (*openblas_get_num_threads_func)(void);
typedef void (*mkl_set_num_threads_func)(int);
typedef int (*mkl_get_max_threads_func)(void);
typedef void (*blis_set_num_threads_func)(int);
typedef int (*blis_get_num_threads_func)(void);

/**
 * @brief Manages the BLAS library and its threading functionality.
 *
 * This class detects the linked BLAS library (OpenBLAS, MKL, BLIS, etc.),
 * provides functions to get and set the number of BLAS threads, and handles
 * the unloading of the BLAS library when no longer needed.
 */
class BLASLibraryManager
{
public:
  /**
   * @brief Constructs the BLASLibraryManager with default settings.
   *
   * Initializes the `blas_handle` to `nullptr`, `blas_type` to "Unknown",
   * and `original_num_threads` to 1.
   */
  BLASLibraryManager()
      : blas_handle(nullptr), blas_type("Unknown"), original_num_threads(1) {}
  /**
   * @brief Detects which BLAS library is linked by checking for specific function symbols.
   *
   * This function loads the BLAS library and determines whether it's OpenBLAS, MKL, or BLIS.
   */
  void detect_lib();
  /**
   * @brief Gets the number of BLAS threads.
   *
   * Calls the appropriate function for the detected BLAS library to return the number of threads.
   *
   * @return The number of BLAS threads currently in use.
   */
  int get_num_threads() const;
  /**
   * @brief Sets the number of BLAS threads.
   *
   * Calls the appropriate function for the detected BLAS library to set the number of threads.
   *
   * @param num_threads The number of threads to set for BLAS operations.
   */
  void set_num_threads(int num_threads);
  /**
   * @brief Destructor for BLASLibraryManager.
   *
   * Restores the original number of BLAS threads and unloads the BLAS library when the manager is destroyed.
   */
  ~BLASLibraryManager()
  {
    // Restore to previous number of blas threads
    set_num_threads(original_num_threads);
    if (blas_handle)
    {
#ifdef WIN32
      FreeLibrary((HMODULE)blas_handle); // On Windows, unload the library
#else
      dlclose(blas_handle); // On POSIX systems, unload the library
#endif
    }
  };

private:
  void *blas_handle;
  std::string blas_type;
  int original_num_threads;

  openblas_get_num_threads_func openblas_get_num_threads = nullptr;
  openblas_set_num_threads_func openblas_set_num_threads = nullptr;
  mkl_get_max_threads_func mkl_get_max_threads = nullptr;
  mkl_set_num_threads_func mkl_set_num_threads = nullptr;
  blis_set_num_threads_func blis_set_num_threads = nullptr;
  blis_get_num_threads_func blis_get_num_threads = nullptr;

  /**
   * @brief Tries to load the BLAS shared library (DLL or SO file).
   *
   * This function attempts to load the appropriate BLAS library from the system.
   *
   * @return A handle to the loaded BLAS library, or `nullptr` if loading fails.
   */
  void *load_blas_library();

  /**
   * @brief Resolves a function symbol from the loaded BLAS library.
   *
   * This function finds the address of a specific function in the BLAS library (e.g., `set_num_threads`).
   *
   * @param symbol The name of the symbol to resolve.
   * @return A pointer to the resolved function, or `nullptr` if the symbol could not be found.
   */
  void *resolve_symbol(const char *symbol);
};