
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

void Chunker::get_num_threads() {
    num_available_threads = std::thread::hardware_concurrency();
}

int Chunker::get_num_workers() {
    return num_workers;
}

int Chunker::get_openmp_threads() {
    return num_threads;
}

int Chunker::get_chunk_size() {
    return chunk_size;
}

int Chunker::get_total_workers() {
    return _num_files;
}

void Chunker::get_num_cpus() {
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    num_available_workers = sysinfo.dwNumberOfProcessors;
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
    num_available_workers = cpus;
    #else
    num_available_workers = sysconf(_SC_NPROCESSORS_ONLN);
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

double Chunker::get_data_size() {
    double matrix_size = std::exp(std::log(_num_poi) + std::log(_num_ind));
    int max_num_matrix = 3;
    int float_size = 8; // 8 bytes per number assuming 64-bit numbers
    double data_size = std::exp(std::log(matrix_size) + std::log(float_size) + std::log(max_num_matrix)) + master_thread_memory;
    return data_size;
}
void Chunker::estimate_chunks() {
    num_threads = std::floor(num_available_threads/num_available_workers);
    this->num_workers = std::min(std::min(_num_files, num_available_workers), _max_workers);

    if (num_workers == num_available_workers) {
        this->num_workers--; // keep 1 process idle 
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

    if (os == "Windows") {
        this->num_workers = 1;
    }

    double data_size = get_data_size();
    int chunks = std::ceil(data_size / static_cast<double>(memfree));
    int chunked_dim1 = std::ceil((double)_num_poi / chunks);
    int chunked_parallel_dim1 = std::floor((double)chunked_dim1/num_workers);

    if (_poi_block_size > 0 && chunked_parallel_dim1 > _poi_block_size) {
        chunked_parallel_dim1 = _poi_block_size;
    }

    if (chunked_parallel_dim1 > _num_poi) {
        chunked_parallel_dim1 = _num_poi;
    }
    chunk_size = chunked_parallel_dim1;
}

ChunkConfig Chunker::estimate_num_files(int num_poi, int num_ind, int usr_threads, float usr_mem) {
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
    size_t len_mem = sizeof(memsize);
    if (sysctlbyname("hw.memsize", &memsize, &len_mem, NULL, 0) == 0) {
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
	
	int max_processes;
	unsigned long long max_mem;
	unsigned long long single_poi_size;
	int max_chunk;
	int blas_threads = 2;
	
	//total available threads minus 1 to run OS and 1 to run fastR main process divided by count of BLAS threads per process
	
	if(usr_threads != -1) {
		max_processes = (usr_threads - 2) / blas_threads;
	} else {
		max_processes = (num_procs - 2) / blas_threads;
	}
	
	if(usr_mem != -1.0) {
		max_mem = static_cast<unsigned long long>((usr_mem * 1024 * 1024 *1024) - 524288000ULL);
	} else {
		max_mem = memfree - 524288000ULL;
	}
	
	single_poi_size = 4 * num_ind * float_size;
	max_chunk = max_mem / (single_poi_size * max_processes);
	
	if(max_mem > (single_poi_size * num_poi)) {
		max_chunk = num_poi / max_processes;
	}

    ChunkConfig chunk_config;
    chunk_config.num_files = (num_poi / max_chunk) + ((num_poi % max_chunk) > 0 ? 1 : 0);
    chunk_config.num_poi = max_chunk;
    return chunk_config;
}
