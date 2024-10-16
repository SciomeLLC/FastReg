#include <cstdio>
#include <cstdlib>
#include <cstring>

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

class BLASLibraryManager {
public:
  BLASLibraryManager()
      : blas_handle(nullptr), blas_type("Unknown"), original_num_threads(1) {}
  //////////////////////////////////////////////////
  // @brief detects which BLAS library is linked by checking for specific
  // function symbols
  // @return  none
  //////////////////////////////////////////////////
  void detect_lib();
  //////////////////////////////////////////////////
  // @brief gets the number of BLAS threads by calling the underlying BLAS
  // libraries function
  // @return integer count of BLAS threads
  //////////////////////////////////////////////////
  int get_num_threads() const;
  //////////////////////////////////////////////////
  // @brief sets the number of BLAS threads by calling the underlying BLAS
  // libraries function
  // @return none
  //////////////////////////////////////////////////
  void set_num_threads(int num_threads);

  ~BLASLibraryManager() {
    // Restore to previous number of blas threads
    set_num_threads(original_num_threads);
    if (blas_handle) {
#ifdef WIN32
      FreeLibrary((HMODULE)blas_handle); // On Windows, unload the library
#else
      dlclose(blas_handle); // On POSIX systems, unload the library
#endif
    }
  };

private:
  void *blas_handle;
  const char *blas_type;
  int original_num_threads;

  openblas_get_num_threads_func openblas_get_num_threads = nullptr;
  openblas_set_num_threads_func openblas_set_num_threads = nullptr;
  mkl_get_max_threads_func mkl_get_max_threads = nullptr;
  mkl_set_num_threads_func mkl_set_num_threads = nullptr;
  blis_set_num_threads_func blis_set_num_threads = nullptr;
  blis_get_num_threads_func blis_get_num_threads = nullptr;

  //////////////////////////////////////////////////
  // @brief tries to load the BLAS dll / shared lib
  // @return pointer to BLAS library
  //////////////////////////////////////////////////
  void *load_blas_library();

  //////////////////////////////////////////////////
  // @brief Resolve the symbol for a function from the loaded BLAS library
  // @symbol name of the symbol to resolve
  // @return pointer to BLAS library
  //////////////////////////////////////////////////
  void *resolve_symbol(const char *symbol);
};