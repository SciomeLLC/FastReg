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
    int num_threads, chunk_size, _num_poi, _num_ind, _poi_block_size, _max_threads, num_procs, num_chunks, _num_files;
    unsigned long long memfree, master_thread_memory;
    Chunker(int num_poi, int num_ind, int max_threads, int poi_block_size, int num_files) {
        _num_poi = num_poi;
        _num_ind = num_ind;
        _max_threads = max_threads;
        _poi_block_size = poi_block_size;
        master_thread_memory = 524288000ULL; // 500mb
        memfree = getTotalSystemMemory();
        _num_files = num_files;
        chunk_size = 100;
        get_num_cpus();
        Rcpp::Rcout << num_procs << " processors detected." << std::endl;
        Rcpp::Rcout << "Free memory: " << memfree/(1024*1024*1024) << "GB" << std::endl;
        estimate_chunks_threads();
    }
    unsigned long long getTotalSystemMemory();
    int get_chunk_size();
    int get_parallel_chunk_size();
    int get_threads();
    int get_procs();
    void get_num_cpus();
    private:
    void estimate_chunks_threads();
    static ChunkConfig estimate_num_files(int num_poi, int num_ind);
};