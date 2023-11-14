#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <h5file.h>
#include <fr_matrix.h>
#include <regression.h>
#include <config.h>
#if !defined(__APPLE__) && !defined(__MACH__)
  #include <omp.h>
#endif
class Job {
    public:
    int id;
    FRMatrix covar_matrix;
    FRMatrix pheno_matrix;
    FRMatrix interactions_matrix;
    Config config;
    double nonconvergence_status;
    double total_filtered_pois;
    int num_threads;
    int num_poi;
    int chunk_size;
    int stratum;
    int poi_reading_time = 0;
    int memory_allocation_time = 0;
    int regression_time = 0;
    int file_writing_time = 0;

    std::vector<std::string> strat_individuals;
    std::vector<std::string> poi_names;
    int num_poi_blocks;
    Job(
        int id, 
        const FRMatrix& covar_matrix, 
        const FRMatrix& pheno_matrix, 
        const FRMatrix& interactions_matrix,
        const Config& config, 
        double nonconvergence_status, 
        double total_filtered_pois, 
        int num_threads, 
        int num_poi, 
        int chunk_size,
        int stratum,
        const std::vector<std::string> strat_individuals,
        const std::vector<std::string> poi_names,
        int num_poi_blocks
    )
    : id(id), 
      covar_matrix(covar_matrix), 
      pheno_matrix(pheno_matrix), 
      interactions_matrix(interactions_matrix),
      config(config), 
      nonconvergence_status(nonconvergence_status), 
      total_filtered_pois(total_filtered_pois), 
      num_threads(num_threads), 
      num_poi(num_poi), 
      chunk_size(chunk_size), 
      stratum(stratum), 
      strat_individuals(strat_individuals),
      poi_names(poi_names),
      num_poi_blocks(num_poi_blocks) {};
    Job(){};
};