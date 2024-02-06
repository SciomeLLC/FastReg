
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


void Chunker::get_num_cpus() {
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    num_procs = sysinfo.dwNumberOfProcessors;
    #elif defined(__APPLE__)
    int mib[2];
    size_t len = sizeof(mib[0]);
    mib[0] = CTL_HW;
    mib[1] = HW_AVAILCPU;
    int cpus = 0;
    sysctl(mib, 2, &cpus, &len, NULL, 0);
    if (cpus < 1) 
    {
        mib[1] = HW_NCPU;
        sysctl(mib, 2, &cpus, &len, NULL, 0);
        if (cpus < 1)
            cpus = 1;
    }
    num_procs = cpus;
    #else
    num_procs = sysconf(_SC_NPROCESSORS_ONLN);
    #endif
}

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
    // Rcpp::Rcout << "Estimating block size" << std::endl;
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
    num_threads = 1;
    if (os == "Darwin") {
        Rcpp::warning("**********\n " \
                    "Mac detected - using 1 thread.\n" \
                    "**********");
        num_threads = 1;
    }
    double matrix_size = std::exp(std::log(_num_poi) + std::log(_num_ind));
    
    int num_workers = 0;
    if (num_procs > _num_files) {
        num_workers = _num_files;
    } else {
        num_workers = num_procs;
    }

    int max_num_matrix = 4 * num_workers;
    int float_size = 8; // 8 bytes per number assuming 64-bit numbers
    double data_size = std::exp(std::log(matrix_size) + std::log(float_size) + std::log(max_num_matrix));
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
    int parallel_chunk_size = 0;
    if(num_procs > _num_files) {
        parallel_chunk_size = (int) std::floor((double)chunk_size/(double)_num_files);
    }
    else {
        parallel_chunk_size = (int) std::floor((double)chunk_size/(double)num_procs);
    }
    return parallel_chunk_size;
}

ChunkConfig Chunker::estimate_num_files(int num_poi, int num_ind) {
    int num_procs = 0;
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    num_procs = sysinfo.dwNumberOfProcessors;
    #elif defined(__APPLE__)
    int mib[2];
    size_t len = sizeof(mib[0]);
    mib[0] = CTL_HW;
    mib[1] = HW_AVAILCPU;
    int cpus = 0;
    sysctl(mib, 2, &cpus, &len, NULL, 0);
    if (cpus < 1) 
    {
        mib[1] = HW_NCPU;
        sysctl(mib, 2, &cpus, &len, NULL, 0);
        if (cpus < 1)
            cpus = 1;
    }
    num_procs = cpus;
    #else
    num_procs = sysconf(_SC_NPROCESSORS_ONLN);
    #endif
    unsigned long long memfree = 0.0;
    #ifdef _WIN32
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    memfree = status.ullAvailPhys;
    #elif defined(__APPLE__)
    uint64_t memsize;
    size_t len = sizeof(memsize);
    if (sysctlbyname("hw.memsize", &memsize, &len, NULL, 0) == 0) {
        memfree = (unsigned long long) memsize;
    }
    #else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    memfree = pages * page_size;
    #endif
    
    int threads = std::thread::hardware_concurrency();

    int max_num_matrix = 4 * num_procs;
    double matrix_size = std::exp(std::log(num_poi) + std::log(num_ind));
    int float_size = 8; // 8 bytes per number assuming 64-bit numbers
    double data_size = std::exp(std::log(matrix_size) + std::log(float_size) + std::log(max_num_matrix));
    double chunks = (data_size) / static_cast<double>(memfree);
    int chunked_dim1 = std::floor(num_poi / chunks);

    ChunkConfig chunk_config;
    chunk_config.num_files = std::ciel(chunks);
    chunk_config.num_poi = chunked_dim1;
    return chunk_config;
}   

int Chunker::get_threads() {
    return num_threads;
};

int Chunker::get_procs() {
    return num_procs;
};