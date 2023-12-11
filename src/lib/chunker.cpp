
#include <chunker.h>
#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif


#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif


unsigned long long Chunker::getTotalSystemMemory() {
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

void Chunker::estimate_chunks_threads() {
    // TEMPORARY HARD CODE
    num_procs = 10;
    Rcpp::Rcout << "Estimating block size" << std::endl;
    num_threads = std::thread::hardware_concurrency();
    Rcpp::Rcout << "Totals threads available: " << num_threads << std::endl;
    if (_max_threads != 0 && num_threads > _max_threads) {
        num_threads = _max_threads;
    } else {
        num_threads--;
        Rcpp::Rcout << "Keeping 1 thread idle, num_threads: " << num_threads << std::endl;
    }

    
    if (_max_threads > num_procs || num_threads > num_procs) {
        num_threads = num_threads - num_procs;
    }

    num_threads = std::floor(num_threads/num_procs);

    std::string os = "";
    #ifdef _WIN32
    os = "Windows";
    #elif defined(__APPLE__)
    os = "Darwin";
    #else
    os = "Linux";
    #endif
    // Hard code
    num_threads = 2;
    if (os == "Darwin") {
        Rcpp::warning("**********\n " \
                    "Mac detected - using 1 thread.\n" \
                    "**********");
        num_threads = 1;
    }
    double matrix_size = std::exp(std::log(_num_poi) + std::log(_num_ind));
    int max_num_matrix = 3;
    int float_size = 8; // 8 bytes per number assuming 64-bit numbers
    double data_size = std::exp(std::log(matrix_size) + std::log(float_size) + std::log(max_num_matrix)) + 5;
    double chunks = (data_size) / static_cast<double>(memfree);
    int chunked_dim1 = std::floor(_num_poi / chunks);

    if (_poi_block_size > 0 && chunked_dim1 > _poi_block_size) {
        chunked_dim1 = _poi_block_size;
    }

    if (chunked_dim1 > _num_poi) {
        chunked_dim1 = _num_poi;
    }
    chunk_size = chunked_dim1;
}

int Chunker::get_chunk_size() {
    return chunk_size;
}

int Chunker::get_parallel_chunk_size() {
    return (int) std::floor((double)chunk_size/(double)num_procs);
}

int Chunker::get_threads() {
    return num_threads;
};

int Chunker::get_procs() {
    return num_procs;
};