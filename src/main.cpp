#define ARMA_WARN_LEVEL 0
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <ctime>
#include <job.h>
#include <chunker.h>
#include <config.h>
#include <covariate.h>
#include <fr_matrix.h>
#include <h5file.h>
#include <regression.h>
#include <strata.h>

#ifndef UTILS_H
#include <utils.h>
#endif

#include <worker.h>
#if !defined(__APPLE__) && !defined(__MACH__)
  #include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#  if __cplusplus >= 201703L && __has_include(<filesystem>)
#    include <filesystem>
namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#  endif
#endif


// [[Rcpp::export]]
void FastRegCpp(
    const std::string phenotype,
    const std::string regression_type,
    const std::string pvalue_dist,
    bool output_exclude_covar,
    double maf_threshold,
    double hwe_threshold,
    bool no_intercept,
    double colinearity_rsq,
    int poi_block_size,
    int max_iter,
    double rel_conv_tolerance,
    double abs_conv_tolderance, 
    int max_threads,
    const std::string pheno_file,
    const std::string pheno_rowname_cols,
    const std::string pheno_file_delim,
    const std::string covar_file,
    const std::string covar_rowname_cols,
    const std::string covar_file_delim,
    const std::string poi_file, 
    const std::string poi_file_delim,
    const std::string poi_file_format,
    const std::string poi_type,
    const std::string poi_effect_type,
    const Rcpp::StringVector covariates,
    const Rcpp::StringVector covariate_type,
    const Rcpp::LogicalVector covariate_standardize,
    const Rcpp::StringVector covariate_levels,
    const Rcpp::StringVector covariate_ref_level,
    const Rcpp::StringVector POI_covar_interactions_str,
    const Rcpp::StringVector split_by_str,
    const std::string output_dir,
    bool compress_results
    ) {
    Config config(
        phenotype,
        regression_type,
        pvalue_dist,
        output_exclude_covar,
        maf_threshold,
        hwe_threshold,
        no_intercept,
        colinearity_rsq,
        poi_block_size,
        max_iter,
        rel_conv_tolerance,
        abs_conv_tolderance, 
        max_threads,
        pheno_file,
        pheno_rowname_cols,
        pheno_file_delim,
        covar_file,
        covar_rowname_cols,
        covar_file_delim,
        poi_file, 
        poi_file_delim,
        poi_file_format,
        poi_type,
        poi_effect_type,
        covariates,
        covariate_type,
        covariate_standardize,
        covariate_levels,
        covariate_ref_level,
        POI_covar_interactions_str,
        split_by_str,
        output_dir,
        compress_results
    );
    FRMatrix pheno_df(config.pheno_file, config.pheno_file_delim, config.pheno_rowname_cols);
    FRMatrix covar_df(config.covar_file, config.covar_file_delim, config.covar_rowname_cols);
    H5File poi(config.POI_file);
    poi.open_file();
    poi.get_values_dataset_id();
    poi.get_POI_names();
    poi.get_POI_individuals();
    
    // Clean up previous run 
    if(dir_exists(config.output_dir)) {
        delete_dir(config.output_dir);
    }

    // Find common individuals
    std::vector<std::string> poi_names = poi.names;
    std::vector<std::string> common_ind = intersect_row_names(pheno_df.sort_map(true), covar_df.sort_map(true));
    std::vector<std::string> intersected_ind = intersect_row_names(common_ind, poi.individuals);

    Rcout << intersected_ind.size() << " unique subjects were found to be common in pheno.file, covar.file, and POI.file" << std::endl;
    if (intersected_ind.empty()) {
        stop("No overlapping individuals found in POI, pheno, and covar files");
    }

    // if (!config.POI_subset_file.empty()) {
        
    //     FRMatrix poi_subset(config.POI_subset_file, config.POI_subset_file_delim, config.POI_subset_rowname_col);
    //     std::vector<std::string> poi_names = intersect_row_names(poi_subset.str_data[0], poi.names);
    // }

    int num_poi = poi_names.size();
    if (num_poi == 0) {
        stop("No overlapping individuals found in POI, pheno, covar files");
    }

    // Stratify data
    Strata stratums;
    stratums.stratify(config.split_by, covar_df, intersected_ind);
    Rcpp::Rcout << "Successfully Stratified" << std::endl;


    for (int stratum = 0; stratum < stratums.nstrata; ++stratum) {
        std::string outfile_suffix = stratums.ids[stratum];
        if (!config.split_by[0].empty()) {
            Rcpp::Rcout << "Processing stratum: " << outfile_suffix.substr(1) << std::endl;
        }
        std::vector<std::string> ind_set = stratums.index_list[outfile_suffix];
        int ct = 0;
        std::vector<int> ind_set_idx(ind_set.size());
        for (std::string ind: ind_set) {
            ind_set_idx[ct] = pheno_df.get_row_idx(ind);
            ct++;
        }

        // initialize matrices
        FRMatrix pheno_matrix = pheno_df.get_submat_by_cols(ind_set_idx, {config.phenotype}); // check col names

        Rcout << "Creating Design Matrix for stratum: " << stratum + 1 << std::endl;
        FRMatrix covar_matrix = create_design_matrix(
            covar_df,
            config.covs,
            config.no_intercept,
            config.colinearity_rsq
        ); // n individuals x # covariates

        FRMatrix covar_poi_interaction_matrix;
        covar_poi_interaction_matrix.data = arma::mat(covar_matrix.data.n_rows, 1, arma::fill::ones);
        covar_poi_interaction_matrix.row_names = covar_matrix.row_names;
        covar_poi_interaction_matrix.col_names = {{"poi", 0}};
        create_Z_matrix(covar_matrix, config.POI_covar_interactions, covar_poi_interaction_matrix); // n individuals x 1 or 1 + num interacting poi covars

        // int num_threads, chunk_size;
        Chunker chunker = Chunker(num_poi, ind_set.size(), config.max_threads, config.poi_block_size);
        std::vector<int> nan_idx; 
        std::vector<std::string> ind_set_filtered; 
        // identify missing values for covar, pheno matrix
        for (size_t i = 0; i < covar_matrix.data.n_rows; i++) {

            arma::uvec covar_nan_idx = arma::find_nonfinite(covar_matrix.data.row(i));
            arma::uvec pheno_nan_idx = arma::find_nonfinite(pheno_matrix.data.row(i));
            if (covar_nan_idx.size() > 0 || pheno_nan_idx.size() > 0) {
                nan_idx.push_back(i);
            }
            else {
                ind_set_filtered.push_back(ind_set[i]);
            }
        }

        // remove from covar, pheno
        covar_matrix.data.shed_rows(arma::conv_to<arma::uvec>::from(nan_idx));
        pheno_matrix.data.shed_rows(arma::conv_to<arma::uvec>::from(nan_idx));
        covar_poi_interaction_matrix.data.shed_rows(arma::conv_to<arma::uvec>::from(nan_idx));
        covar_matrix.row_names = std::unordered_map<std::string, int>(ind_set_filtered.size());
        pheno_matrix.row_names = std::unordered_map<std::string, int>(ind_set_filtered.size());

        for (size_t j = 0; j < covar_matrix.data.n_rows; j++) {
            covar_matrix.row_names[ind_set_filtered[j]] = j;
            pheno_matrix.row_names[ind_set_filtered[j]] = j;
            covar_poi_interaction_matrix.row_names[ind_set_filtered[j]] = j;
        }

        std::vector<std::string> strat_individuals(ind_set_filtered.size());
        std::transform(
            ind_set_filtered.begin(),
            ind_set_filtered.end(),
            strat_individuals.begin(),
            [&intersected_ind](const std::string& elem) {
                return intersected_ind[
                    std::distance(
                        intersected_ind.begin(),
                        std::find(intersected_ind.begin(), intersected_ind.end(), elem)
                    )
                ];
            }
        );

        int num_threads = chunker.get_threads();
        int parallel_chunk_size = chunker.get_parallel_chunk_size();
        int num_processes = chunker.get_procs();
        // Rcout << "POI block size: " << parallel_chunk_size << std::endl;
        // int num_poi_blocks = (int) std::ceil((double)num_poi/(double)chunk_size);
        // Rcout << "POIs will be processed in " << num_poi_blocks << " blocks each of size " << chunk_size << std::endl;

        Rcout << "processes: " << num_processes << std::endl;
        double nonconvergence_status = 0.0;
        double total_filtered_pois = 0.0;
        
        int num_parallel_poi_blocks = (int) std::ceil((double)num_poi/(double)parallel_chunk_size);
        Rcout << "POIs will be processed in " << num_parallel_poi_blocks << " blocks each of size " << parallel_chunk_size << std::endl;
        Rcpp::Rcout << "Closing main thread h5 file" << std::endl;
        poi.close_all();
        // std::atomic<double> nonconvergence_status_accumulator(0.0);
        // std::atomic<double> total_filtered_pois_accumulator(0.0);
        std::mutex mtx;
        std::condition_variable cv;
        std::queue<Job> job_queue;
        for (int i = 0; i < num_parallel_poi_blocks; i++) {
            job_queue.emplace(Job(
                i, 
                covar_matrix, 
                pheno_matrix, 
                covar_poi_interaction_matrix,
                config, 
                nonconvergence_status, 
                total_filtered_pois, 
                num_threads, 
                num_poi, 
                parallel_chunk_size,
                stratum,
                strat_individuals,
                poi_names, 
                num_parallel_poi_blocks)
            );
        }

        if (num_parallel_poi_blocks < num_processes) {
            num_processes = num_parallel_poi_blocks;
        }
        std::vector<std::thread> threads;
        auto start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        
        std::atomic<int> memory_allocation_time(0);
        std::atomic<int> file_writing_time(0);
        std::atomic<int> poi_reading_time(0);
        std::atomic<int> regression_time(0);
        Rcpp::Rcout << "Total processes spawned: " << num_processes << " at: " << std::ctime(&start) << std::endl;
        for (int i = 0; i < num_processes; i++) {
            threads.emplace_back([&]() {
                try {
                Worker worker(i, job_queue, mtx, cv, memory_allocation_time, file_writing_time, poi_reading_time, regression_time);
                worker();
                }  catch (const std::exception& ex) {
                    std::cerr << "Worker thread caught an exception: " << ex.what() << std::endl;
                }
            });
        }

        // Notify all threads that there are no more jobs
        {
            std::lock_guard<std::mutex> lock(mtx);
            cv.notify_all();
        }

        // Wait for all threads to finish
        for (std::thread& thread : threads) {
            thread.join();
        }
        

        Rcpp::Rcout << "Timing Summary: " << std::endl;
        Rcpp::Rcout << "Reading HDF5: " << poi_reading_time / 1000.0 << "s" << std::endl;
        Rcpp::Rcout << "Writing results: " << file_writing_time / 1000.0 << "s" << std::endl;
        Rcpp::Rcout << "Memory allocation: " << memory_allocation_time / 1000.0 << "s" << std::endl;
        Rcpp::Rcout << "Regression: " << regression_time / 1000.0 << "s" << std::endl;
        auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        // double noncovergence_percent = (nonconvergence_status_accumulator.load() / total_filtered_pois_accumulator.load()) * 100;
        // if (noncovergence_percent > 0.0) {
        //     Rcpp::Rcout << nonconvergence_status << " out of " << total_filtered_pois << " (" << std::setprecision(2) << noncovergence_percent << "%) POIs did not meet relative and absolute convergence threshold." << std::endl;
        //     Rcpp::Rcout << "See convergence_" << stratum << ".tsv for additional details." << std::endl;
        // }
    }
    if (config.compress_results) {
        FRMatrix::zip_results(config.output_dir);
    }
}
