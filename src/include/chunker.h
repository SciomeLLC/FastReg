

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif


#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <thread>

unsigned long long getTotalSystemMemory() {
#ifdef _WIN32
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullAvailPhys;
#elif defined(__APPLE__)
    uint64_t memsize;
    size_t len = sizeof(memsize);
    if (sysctlbyname("hw.memsize", &memsize, &len, NULL, 0) == 0) {
        return (unsigned long long) memsize;
    } else {
        return 0;
    }
#else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#endif
}


std::vector<int> estimate_poi_block_size(int num_poi, int num_ind, std::string poi_type, int max_threads, int poi_block_size) {
    std::vector<int> res;
    Rcpp::Rcout << "Estimating block size" << std::endl;
    int num_threads = std::thread::hardware_concurrency();
    if (max_threads != 0 && num_threads > max_threads) {
        num_threads = max_threads;
    } else {
        num_threads--;
        Rcpp::Rcout << "Keeping 1 thread idle, num_threads: " << num_threads << std::endl;
    }

    std::string os = "";
    #ifdef _WIN32
    os = "Windows";
    #elif defined(__APPLE__)
    os = "Darwin";
    #else
    os = "Linux";
    #endif
    if (os == "Darwin") {
        Rcpp::warning("**********\n " \
                    "Mac detected - using 1 thread.\n" \
                    "**********");
        num_threads = 1;
    }
    unsigned long long memfree = getTotalSystemMemory();
    
    Rcpp::Rcout << "Free memory: " << memfree/(1024*1024*1024) << "GB" << std::endl;

    double matrix_size = std::exp(std::log(num_poi) + std::log(num_ind));
    int max_num_matrix = 6;
    int float_size = 8; // 8 bytes per number assuming 64-bit numbers
    double data_size = std::exp(std::log(matrix_size) + std::log(float_size) + std::log(max_num_matrix));
    unsigned long long master_thread_memory = 524288000ULL; // 500mb
    double chunks = (data_size + master_thread_memory) / static_cast<double>(memfree);
    int chunked_dim1 = std::floor(num_poi / chunks);

    if (poi_block_size > 0 && chunked_dim1 > poi_block_size) {
        chunked_dim1 = poi_block_size;
    }

    if (chunked_dim1 > num_poi) {
        chunked_dim1 = num_poi;
    }
    res.push_back(chunked_dim1);
    res.push_back(num_threads);
    return res;
}
