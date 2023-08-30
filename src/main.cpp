
#define ARMA_WARN_LEVEL 0
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <utils.h>
#include <regression.h>
#include <chunker.h>
#include <h5file.h>
#include <fr_matrix.h>
#include <covariate.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <config.h>
#include <strata.h>
#include <chrono>
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

bool dir_exists(const std::string& path) {
    fs::path directory(path);
    return fs::is_directory(directory);
}

void delete_dir(const std::string& path) {
    fs::path directory(path);
    
    if (fs::exists(directory) && fs::is_directory(directory)) {
        fs::remove_all(directory);
        std::cout << "Directory deleted: " << path << std::endl;
    }
    else {
        std::cout << "Directory does not exist: " << path << std::endl;
    }
}

template<typename T, typename U>
bool isin(const T& value, const U& container) {
    return std::find(std::begin(container), std::end(container), value) != std::end(container);
}

std::vector<std::string> intersect_row_names(const std::vector<std::string>& a, const std::vector<std::string>& b) {
    std::vector<std::string> result;
    result.reserve(std::min(a.size(), b.size()));
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));
    return result;
}

std::vector<std::string> set_diff(const std::vector<std::string>& a, const std::vector<std::string>& b) {
    std::vector<std::string> result;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));
    return result;
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
    double file_writing_time = 0.0;
    double poi_reading_time = 0.0;
    double regression_time = 0.0;
    double memory_allocation_time = 0.0;
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

        int num_threads, chunk_size;
        std::vector<int> chunk_size_num_threads = estimate_poi_block_size(num_poi, ind_set.size(), config.POI_type, config.max_threads, config.poi_block_size);
        chunk_size = chunk_size_num_threads[0];
        num_threads = chunk_size_num_threads[1];
        
        Rcout << "POI block size: " << chunk_size << std::endl;
        Rcout << "threads: " << num_threads << std::endl;
        int num_poi_blocks = (int) std::ceil((double)num_poi/(double)chunk_size);
        Rcout << "POIs will be processed in " << num_poi_blocks << " blocks each of size " << chunk_size << std::endl;

        
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

        FRMatrix poi_matrix;
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

        poi.set_memspace(poi.individuals.size(), chunk_size);
        double nonconvergence_status = 0.0;
        double total_filtered_pois = 0.0;
        for (int block = 0; block < num_poi_blocks; block ++) {
            Rcpp::Rcout << "Processing POIs block: " << block + 1 << std::endl;

            int start_chunk = block*chunk_size;
            int end_chunk = start_chunk + chunk_size;
            if (end_chunk > (int)poi_names.size()) {
                end_chunk = (int)poi_names.size();
                poi.set_memspace(poi.individuals.size(), end_chunk - start_chunk);
            }
            std::vector<std::string> poi_names_chunk(poi_names.begin() + start_chunk, poi_names.begin() + end_chunk);
            #if !defined(__APPLE__) && !defined(__MACH__)
            omp_set_num_threads(num_threads);
            #endif
            auto start_time = std::chrono::high_resolution_clock::now();
            poi.get_POI_matrix(poi_matrix, poi.individuals, poi_names_chunk, chunk_size);

            std::vector<std::string> srt_cols_2 = poi_matrix.sort_map(false);
            std::vector<std::string> drop_rows = set_diff(poi.individuals, strat_individuals);

            int num_dropped = poi.individuals.size() - strat_individuals.size();
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
            poi_reading_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
            Rcpp::Rcout << "Reading POI timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";

            if (config.POI_type == "genotype") {
                Rcout << "Filtering MAF and HWE" << std::endl;

                FRMatrix filtered = filter_poi(poi_matrix, config.maf_threshold, config.hwe_threshold);
                arma::uvec filtered_col = arma::find(filtered.data.row(5) == 0);

                if (filtered.data.n_cols == 0 || filtered_col.n_elem == poi_matrix.data.n_cols) {
                    Rcout << "no POI passed filtering" << std::endl;
                    continue;
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

                Rcout << "transforming POI and writing statistics" << std::endl;
                srt_cols_2 = poi_matrix.sort_map(false);
                transform_poi(poi_matrix, config.POI_effect_type);
                start_time = std::chrono::high_resolution_clock::now();
                filtered.write_summary(config.output_dir, "POI_Summary", stratum);
                end_time = std::chrono::high_resolution_clock::now();
                file_writing_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
                Rcpp::Rcout << "Wrting POI Summary timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
                Rcpp::Rcout << "Effective block size after filtering: " << poi_matrix.data.n_cols << std::endl;
            }

            total_filtered_pois += poi_matrix.data.n_cols;
            
            start_time = std::chrono::high_resolution_clock::now();
            // 1 phenotype 
            FRMatrix beta_est;
            FRMatrix se_beta;
            int num_parms = covar_poi_interaction_matrix.data.n_cols + covar_matrix.data.n_cols;
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
            for (auto& col_name : covar_poi_interaction_matrix.col_names) {
                beta_est.row_names[col_name.first] = covar_matrix.col_names.size() + col_name.second;
                se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
                neglog10_pvl.row_names[col_name.first] = beta_est.row_names[col_name.first];
            }

            beta_est.col_names = covar_matrix.row_names;
            se_beta.col_names = beta_est.col_names;
            neglog10_pvl.col_names = beta_est.col_names;
            std::vector<std::string> srt_cols = poi_matrix.sort_map(false);


            // create weight mask matrix    
            arma::umat W2 = arma::umat(poi_matrix.data.n_rows, poi_matrix.data.n_cols, arma::fill::ones);
            for (arma::uword v = 0; v < poi_matrix.data.n_cols; v++) {
                arma::uvec G_na = arma::find_nonfinite(poi_matrix.data.col(v));
                for (arma::uword i = 0; i < G_na.n_elem; i++) {
                    W2(G_na(i), v) = 0;
                    poi_matrix.data(G_na(i), v) = 0;
                }
            }
            
            end_time = std::chrono::high_resolution_clock::now();
            memory_allocation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
            Rcpp::Rcout << "Memory allocation timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";

            start_time = std::chrono::high_resolution_clock::now();
            std::unique_ptr<RegressionBase> regression;


            if (config.regression_type == "logistic") {
                Rcout << "Started logistic regression" << std::endl;
                regression.reset(new LogisticRegression());
            }
            else {
                Rcout << "Started linear regression" << std::endl;
                regression.reset(new LinearRegression ());
            }
            regression->run(
                covar_matrix, 
                pheno_matrix, 
                poi_matrix, 
                covar_poi_interaction_matrix, 
                W2, 
                beta_est, 
                se_beta, 
                neglog10_pvl, 
                beta_rel_errs,
                beta_abs_errs,
                config.max_iter, 
                config.p_value_type == "t.dist"
            );
            end_time = std::chrono::high_resolution_clock::now();
            Rcpp::Rcout << "regression timing for block " << block + 1 <<  "/" << num_poi_blocks << ": " <<std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
            regression_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
            start_time = std::chrono::high_resolution_clock::now();
            FRMatrix::write_results(beta_est, se_beta, neglog10_pvl, W2, srt_cols, config.output_dir, "Results", stratum, config.output_exclude_covar);
            end_time = std::chrono::high_resolution_clock::now();
            file_writing_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
            Rcpp::Rcout << "Writing results for block " << block + 1 <<  "/" << num_poi_blocks << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
            poi_matrix.col_names.clear();

            if (config.regression_type == "logistic") {
                start_time = std::chrono::high_resolution_clock::now();
                FRMatrix::write_convergence_results(beta_est, srt_cols, config.output_dir, "convergence", beta_rel_errs, beta_abs_errs, stratum);
                end_time = std::chrono::high_resolution_clock::now();
                file_writing_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
            }

            arma::colvec convergence = conv_to<colvec>::from((beta_rel_errs > config.rel_conv_tolerance) && (beta_abs_errs > config.abs_conv_tolerance));
            nonconvergence_status += arma::sum(convergence);
        }
        double noncovergence_percent = (nonconvergence_status/total_filtered_pois) * 100;
        if (noncovergence_percent > 0.0) {
            Rcpp::Rcout << nonconvergence_status << " out of " << total_filtered_pois << " (" << std::setprecision(2) << noncovergence_percent << "%) POIs did not meet relative and absolute convergence threshold." << std::endl;
            Rcpp::Rcout << "See convergence_" << stratum << ".tsv for additional details." << std::endl;
        }
    }

    Rcpp::Rcout << "Timing Summary: " << std::endl;
    Rcpp::Rcout << "Reading HDF5: " << poi_reading_time / 1000.0 << "s" << std::endl;
    Rcpp::Rcout << "Writing results: " << file_writing_time / 1000.0 << "s" << std::endl;
    Rcpp::Rcout << "Memory allocation: " << memory_allocation_time / 1000.0 << "s" << std::endl;
    Rcpp::Rcout << "Regression: " << regression_time / 1000.0 << "s" << std::endl;
    if (config.compress_results) {
        FRMatrix::zip_results(config.output_dir);
    }
}
