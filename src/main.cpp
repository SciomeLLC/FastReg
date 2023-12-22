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
#include <unistd.h>

#ifndef UTILS_H
#include <utils.h>
#endif

#if !defined(__APPLE__) && !defined(__MACH__)
  #include <omp.h>
#endif

#ifdef _WIN32
#include <windows.h>
#else 
#include <sys/wait.h>
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

struct ProcResult {
    int timing_results[4] = {0, 0, 0, 0};
    double process_nonconvergence_status = 0.0;
    double process_total_filtered_pois = 0.0;
};

// exception handler - per thread. print the stack trace (check out Boost)

// parameter by reference - Config* config, int* process_id
void process_regression(
    Config& config,
    int process_id,
    int chunk_size,
    int num_poi,
    std::vector<std::string>& poi_names,
    std::vector<std::string>& strat_individuals,
    int num_threads,
    int stratum,
    double& total_filtered_pois,
    int num_poi_blocks,
    FRMatrix& interactions_matrix,
    FRMatrix& covar_matrix,
    FRMatrix& pheno_matrix,
    double& nonconvergence_status,
    int timing_results[4]
    ) {
    FRMatrix poi_matrix;
    H5File poi_child(config.POI_file);
    std::cout << "Processing POIs block: " << process_id + 1 << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    poi_child.open_file(true);
    poi_child.get_values_dataset_id();
    
    std::cout << "started Reading POI individual ids for: " << process_id << std::endl;
    poi_child.get_POI_individuals();
    std::cout << "completed Reading individual ids names for: " << process_id << std::endl;
    std::cout << "started Reading POI names for: " << process_id << std::endl;
    poi_child.get_POI_names();
    std::cout << "completed Reading POI names for: " << process_id << std::endl;
    poi_child.set_memspace(poi_child.individuals.size(), chunk_size);
    int start_chunk = process_id*chunk_size;
    int end_chunk = start_chunk + chunk_size;
    if (end_chunk > num_poi) {
        end_chunk = num_poi;
        poi_child.set_memspace(poi_child.individuals.size(), end_chunk - start_chunk);
    }
    
    std::cout << "started Reading POIs for: " << process_id << std::endl;
    std::vector<std::string> poi_names_chunk(poi_names.begin() + start_chunk, poi_names.begin() + end_chunk);
    poi_child.get_POI_matrix(poi_matrix, poi_child.individuals, poi_names_chunk, chunk_size);
    std::cout << "Completed Reading POIs for: " << process_id << std::endl;
    poi_child.close_all();
    std::cout << "Closed POI file for: " << process_id << std::endl;
    #if !defined(__APPLE__) && !defined(__MACH__)
    omp_set_num_threads(num_threads);
    #endif
    std::vector<std::string> srt_cols_2 = poi_matrix.sort_map(false);
    std::vector<std::string> drop_rows = set_diff(poi_child.individuals, strat_individuals);

    int num_dropped = poi_child.individuals.size() - strat_individuals.size();
    arma::uvec drop_row_idx(drop_rows.size());
    for (size_t i = 0; i < drop_rows.size(); i++) {
        drop_row_idx[i] = poi_matrix.row_names[drop_rows[i]];
    }
    
    std::unordered_map<std::string, int> new_row_names(strat_individuals.size());
    for(auto& ind : strat_individuals) {
        new_row_names[ind] = poi_matrix.row_names[ind] - num_dropped;
    }
    poi_matrix.data.shed_rows(drop_row_idx);
    poi_matrix.row_names = new_row_names;
    
    srt_cols_2 = poi_matrix.sort_map(false);
    auto end_time = std::chrono::high_resolution_clock::now();
    timing_results[2] = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count(); // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
    //std::cout << "Reading POI timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";

    if (config.POI_type == "genotype") {
        std::cout << "Filtering MAF and HWE" << std::endl;

        FRMatrix filtered = filter_poi(poi_matrix, config.maf_threshold, config.hwe_threshold);
        arma::uvec filtered_col = arma::find(filtered.data.row(5) == 0);

        if (filtered.data.n_cols == 0 || filtered_col.n_elem == poi_matrix.data.n_cols) {
            std::cout << "no POI passed filtering" << std::endl;
            return;
        }
        
        std::vector<std::string> poi_col_names = filtered.sort_map(false);
        int cols_erased = 0;

        for(unsigned int i = 0; i < poi_col_names.size(); i++) {
            if ((unsigned) cols_erased < filtered_col.n_elem && filtered_col[cols_erased] == i) {
                poi_matrix.col_names.erase(poi_col_names[i]);
                cols_erased++;
            }
            else {
                poi_matrix.col_names[poi_col_names[i]] = poi_matrix.col_names[poi_col_names[i]] - cols_erased;
            }
        }

        poi_matrix.data.shed_cols(filtered_col);

        //std::cout << "transforming POI and writing statistics" << std::endl;
        srt_cols_2 = poi_matrix.sort_map(false);
        transform_poi(poi_matrix, config.POI_effect_type);
        start_time = std::chrono::high_resolution_clock::now();
        std::string summary_name = "POI_Summary_" + std::to_string(process_id);
        filtered.write_summary(config.output_dir, summary_name, stratum);
        end_time = std::chrono::high_resolution_clock::now();
        // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
        timing_results[1] = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
        //std::cout << "Wrting POI Summary timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";  
        //std::cout << "Effective block size after filtering: " << poi_matrix.data.n_cols << std::endl;
    }

    total_filtered_pois = poi_matrix.data.n_cols;
    
    start_time = std::chrono::high_resolution_clock::now();
    // 1 phenotype 
    
    std::cout << "Creating matrices for regression for: " << process_id << std::endl;
    FRMatrix beta_est;
    FRMatrix se_beta;
    int num_parms = interactions_matrix.data.n_cols + covar_matrix.data.n_cols;
    beta_est.data = arma::mat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
    arma::colvec beta_rel_errs = arma::colvec(poi_matrix.data.n_cols, arma::fill::zeros);
    arma::colvec beta_abs_errs = arma::colvec(poi_matrix.data.n_cols, arma::fill::zeros);
    se_beta.data = arma::mat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
    
    FRMatrix neglog10_pvl;
    neglog10_pvl.data = arma::mat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);


    for (auto& col_name : covar_matrix.col_names) {
        beta_est.row_names[col_name.first] = col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] = beta_est.row_names[col_name.first];
    }
    for (auto& col_name : interactions_matrix.col_names) {
        beta_est.row_names[col_name.first] = covar_matrix.col_names.size() + col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] = beta_est.row_names[col_name.first];
    }

    beta_est.col_names = covar_matrix.row_names;
    se_beta.col_names = beta_est.col_names;
    neglog10_pvl.col_names = beta_est.col_names;
    
    std::cout << "Completed allocating matrices for regression for: " << process_id << std::endl;
    std::vector<std::string> srt_cols = poi_matrix.sort_map(false);

    
    std::cout << "Creating weigth matrices for regression for: " << process_id << std::endl;
    // create weight mask matrix    
    arma::umat W2 = arma::umat(poi_matrix.data.n_rows, poi_matrix.data.n_cols, arma::fill::ones);
    for (arma::uword v = 0; v < poi_matrix.data.n_cols; v++) {
        arma::uvec G_na = arma::find_nonfinite(poi_matrix.data.col(v));
        for (arma::uword i = 0; i < G_na.n_elem; i++) {
            W2(G_na(i), v) = 0;
            poi_matrix.data(G_na(i), v) = 0;
        }
    }
    
    std::cout << "Completed allocating weigth matrices for regression for: " << process_id << std::endl;
    end_time = std::chrono::high_resolution_clock::now();
    // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
    timing_results[0] = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    //std::cout << "Memory allocation timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";

    start_time = std::chrono::high_resolution_clock::now();
    std::unique_ptr<RegressionBase> regression;


    if (config.regression_type == "logistic") {
        //std::cout << "Started logistic regression" << std::endl;
        regression.reset(new LogisticRegression());
    }
    else {
        //std::cout << "Started linear regression" << std::endl;
        regression.reset(new LinearRegression ());
    }
    std::cout << "Started regression for " << process_id << std::endl;
    regression->run(
        covar_matrix, 
        pheno_matrix, 
        poi_matrix, 
        interactions_matrix, 
        W2, 
        beta_est, 
        se_beta, 
        neglog10_pvl, 
        beta_rel_errs,
        beta_abs_errs,
        config.max_iter, 
        config.p_value_type == "t.dist"
    );
    
    std::cout << "Completed regression for " << process_id << std::endl;
    end_time = std::chrono::high_resolution_clock::now();
    // Rcpp::Rcout << "regression timing for block " << process_id + 1 <<  "/" << num_poi_blocks << ": " <<std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
    // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
    timing_results[3] = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    start_time = std::chrono::high_resolution_clock::now();
    std::string results_file_prefix = "Results_" + std::to_string(process_id);
    FRMatrix::write_results(beta_est, se_beta, neglog10_pvl, W2, srt_cols, config.output_dir, results_file_prefix, stratum, config.output_exclude_covar);
    end_time = std::chrono::high_resolution_clock::now();
    // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
    timing_results[1] += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    // std::cout << "Writing results for block " << process_id + 1 <<  "/" << num_poi_blocks << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
    poi_matrix.col_names.clear();

    if (config.regression_type == "logistic") {
        start_time = std::chrono::high_resolution_clock::now();
        std::string convergence_file_prefix = "convergence_" + std::to_string(process_id);
        FRMatrix::write_convergence_results(beta_est, srt_cols, config.output_dir, convergence_file_prefix, beta_rel_errs, beta_abs_errs, stratum);
        end_time = std::chrono::high_resolution_clock::now();
        // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
        timing_results[1] += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    }

    arma::colvec convergence = conv_to<colvec>::from((beta_rel_errs > config.rel_conv_tolerance) && (beta_abs_errs > config.abs_conv_tolerance));
    nonconvergence_status = arma::sum(convergence);
    auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Completed process " << process_id << " at: " << std::ctime(&end) << std::endl;
}

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
    poi.open_file(true);
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

        Rcpp::Rcout << "processes: " << num_processes << std::endl;
        double nonconvergence_status = 0.0;
        double total_filtered_pois = 0.0;
        
        int num_parallel_poi_blocks = (int) std::ceil((double)num_poi/(double)parallel_chunk_size);
        Rcpp::Rcout << "POIs will be processed in " << num_parallel_poi_blocks << " blocks each of size " << parallel_chunk_size << std::endl;
        Rcpp::Rcout << "Closing main thread h5 file" << std::endl;
        
        poi.close_all();
        // auto start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        int total_nonconvergence_status = 0.0;
        double sum_total_filtered_pois = 0.0;
        int memory_allocation_time = 0;
        int file_writing_time = 0;
        int poi_reading_time = 0;
        int regression_time = 0;
        // Rcpp::Rcout << "Total processes spawned: " << num_processes << " at: " << std::ctime(&start) << std::endl;
        

        int num_processes_started = 0;
        int num_processes_completed = 0;
        #ifdef _WIN32
        std::vector<PROCESS_INFORMATION> process_info(num_processes);
        std::vector<HANDLE> h_read_pipe(num_processes);
        SECURITY_ATTRIBUTES security_attributes;
        security_attributes.nLength = sizeof(SECURITY_ATTRIBUTES);
        security_attributes.bInheritHandle = TRUE;
        security_attributes.lpSecurityDescriptor = NULL;

        for (int i = 0; i < num_parallel_poi_blocks; i++) {
            // HANDLE h_write_pipe;
            // if(!CreatePipe(&h_read_pipe[i], &h_write_pipe, &security_attributes, 0)) {
            //     Rcpp::Rcout << "CreatePipe failed (" << GetLastError() << ")." << std::endl;
            //     return -1;
            // }

            // STARTUPINFO siStartInfo;
            // ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));
            // ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
            // siStartInfo.cb = sizeof(STARTUPINFO);
            // siStartInfo.hStdError = hWritePipe;
            // siStartInfo.hStdOutput = hWritePipe;
            // siStartInfo.dwFlags |= STARTF_USESTDHANDLES;
            // std::string cmd = "regression_process " + std::to_string(i);

            // int rval = CreateProcess(
            //     NULL,
            //     &cmd[0],

            // )
        }
        #else
        // create 2 pipes per process - one each for read and write
        std::vector<int> pipe_file_descriptors(num_processes * 2); 
        std::vector<pid_t> process_ids(num_processes);
        std::vector<bool> has_completed(num_parallel_poi_blocks, false);
        // std::unique_lock<std::mutex> lock(_mtx);
        while(num_processes_completed < num_parallel_poi_blocks) {
            while(num_processes_started < num_parallel_poi_blocks && num_processes_started - num_processes_completed < num_processes) {
                int i = num_processes_started;
                if (pipe(&pipe_file_descriptors[i*2]) == -1) {
                    perror("pipe");
                    return;
                }
                process_ids[i] = fork();
                if(process_ids[i] == -1) {
                    perror("fork");
                    return;
                }

                if(process_ids[i] == 0) { // child process
                    close(pipe_file_descriptors[i*2]); // close read pipe
                    int timing_results[4] = {memory_allocation_time, file_writing_time, poi_reading_time, regression_time};
                    
                    Rcpp::Rcout << "Started process_regression for " << i << std::endl;
                    process_regression(
                        config, 
                        i, 
                        parallel_chunk_size, 
                        num_poi, 
                        poi_names, 
                        strat_individuals, 
                        num_threads, 
                        stratum, 
                        total_filtered_pois, 
                        num_parallel_poi_blocks, 
                        covar_poi_interaction_matrix, 
                        covar_matrix, 
                        pheno_matrix, 
                        nonconvergence_status, 
                        timing_results
                    );
                    Rcpp::Rcout << "Completed process_regression for " << i << std::endl;

                    ProcResult proc_result;
                    size_t num_elements = sizeof(timing_results) / sizeof(timing_results[0]);
                    std::copy(timing_results, timing_results + num_elements, proc_result.timing_results);
                    proc_result.process_nonconvergence_status = nonconvergence_status;
                    proc_result.process_total_filtered_pois = total_filtered_pois;
                    ssize_t bytes_written = write(pipe_file_descriptors[i*2 + 1], &proc_result, sizeof(proc_result));
                    if (bytes_written == -1) {
                        perror("Unable write to pipe fd.");
                    }
                    close(pipe_file_descriptors[i*2 + 1]);
                    _exit(EXIT_SUCCESS);
                    return;
                } else {
                    close(pipe_file_descriptors[i*2 + 1]);
                    num_processes_started++;
                }
            }
            // Check for finished processes
            for(int i = 0; i < num_processes_started; i++) {
                if(process_ids[i] != 0) { // parent process
                    if(!has_completed[i] && waitpid(process_ids[i], NULL, WNOHANG) > 0) {
                        // {memory_allocation_time, file_writing_time, poi_reading_time, regression_time}
                        ProcResult proc_result;
                        ssize_t bytes_read = read(
                            pipe_file_descriptors[i*2], 
                            &proc_result, 
                            sizeof(proc_result)
                        );
                        if (bytes_read == -1) {
                            stop("Unable to read pipe fd.");
                        }
                        close(pipe_file_descriptors[i*2]);

                        memory_allocation_time += proc_result.timing_results[0];
                        file_writing_time += proc_result.timing_results[1];
                        poi_reading_time += proc_result.timing_results[2];
                        regression_time += proc_result.timing_results[3];

                        Rcpp::Rcout << "Memory allocation time for process " << i + 1 << ":" << proc_result.timing_results[0]  << " ms" << std::endl;
                        Rcpp::Rcout << "File writing time for process " << i + 1 << ":" << proc_result.timing_results[1]  << " ms" << std::endl;
                        Rcpp::Rcout << "POI reading time for process " << i + 1 << ":" << proc_result.timing_results[2]  << " ms" << std::endl;
                        Rcpp::Rcout << "Regression time for process " << i + 1 << ":" << proc_result.timing_results[3]  << " ms" << std::endl;

                        total_nonconvergence_status += proc_result.process_nonconvergence_status;
                        sum_total_filtered_pois += proc_result.process_total_filtered_pois;
                        has_completed[i] = true;
                        num_processes_completed++;
                        // delete proc_result;
                    }
                }
            }
        }
        #endif

        Rcpp::Rcout << "Timing Summary: " << std::endl;
        Rcpp::Rcout << "Reading HDF5: " << poi_reading_time / 1000.0 << "s" << std::endl;
        Rcpp::Rcout << "Writing results: " << file_writing_time / 1000.0 << "s" << std::endl;
        Rcpp::Rcout << "Memory allocation: " << memory_allocation_time / 1000.0 << "s" << std::endl;
        Rcpp::Rcout << "Regression: " << regression_time / 1000.0 << "s" << std::endl;
        // auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        double noncovergence_percent = (total_nonconvergence_status / sum_total_filtered_pois) * 100;
        if (noncovergence_percent > 0.0) {
            Rcpp::Rcout << total_nonconvergence_status << " out of " << sum_total_filtered_pois << " (" << std::setprecision(2) << noncovergence_percent << "%) POIs did not meet relative and absolute convergence threshold." << std::endl;
            Rcpp::Rcout << "See convergence_" << stratum << ".tsv for additional details." << std::endl;
        }
    }
    if (config.compress_results) {
        FRMatrix::zip_results(config.output_dir);
    }
}
