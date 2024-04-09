// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <thread>


struct ChunkConfig {
    int num_files;
    int num_poi;
};

class Chunker{
    public:
    int num_threads, chunk_size, _num_poi, _num_ind, _poi_block_size, _max_threads, num_procs, num_chunks, _num_files, _max_workers, num_workers, total_workers_required;
    unsigned long long memfree, master_thread_memory;
    Chunker(int num_poi, int num_ind, int max_openmp_threads, int poi_block_size, int num_files, int max_workers) {
        _num_poi = num_poi;
        _num_ind = num_ind;
        _max_threads = max_openmp_threads;
        _max_workers = max_workers;
        _poi_block_size = poi_block_size;
        master_thread_memory = 524288000ULL; // 500mb
        memfree = getTotalSystemMemory();
        _num_files = num_files;
        chunk_size = 100;
        get_num_cpus();
        
        if (_max_workers == 0) {
            _max_workers = num_available_workers;
        }
        get_num_threads();
        Rcpp::Rcout << num_available_workers << " processors detected." << std::endl;
        Rcpp::Rcout << num_available_threads << " threads detected." << std::endl;
        Rcpp::Rcout << "Free memory: " << memfree/(1024*1024*1024) << "GB" << std::endl;
        estimate_chunks();
    }
    unsigned long long getTotalSystemMemory();
    int get_chunk_size();
    int get_num_workers();
    int get_openmp_threads();
    int get_total_workers();
    static ChunkConfig estimate_num_files(int num_poi, int num_ind, int usr_threads, float usr_mem);
    private:
    int num_available_workers, num_available_threads;
    double get_data_size();
    void estimate_chunks();
    void get_num_threads();
    void get_num_cpus();
};
