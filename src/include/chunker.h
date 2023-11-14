// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <thread>

class Chunker{
    public:
    int num_threads, chunk_size, _num_poi, _num_ind, _poi_block_size, _max_threads, num_procs, num_chunks;
    unsigned long long memfree, master_thread_memory;
    Chunker(int num_poi, int num_ind, int max_threads, int poi_block_size) {
        _num_poi = num_poi;
        _num_ind = num_ind;
        _max_threads = max_threads;
        _poi_block_size = poi_block_size;
        master_thread_memory = 524288000ULL; // 500mb
        memfree = getTotalSystemMemory();
        
        Rcpp::Rcout << "Free memory: " << memfree/(1024*1024*1024) << "GB" << std::endl;
        estimate_chunks_threads();
    }
    unsigned long long getTotalSystemMemory();
    int get_chunk_size();
    int get_parallel_chunk_size();
    int get_threads();
    int get_procs();
    private:
    void estimate_chunks_threads();
};