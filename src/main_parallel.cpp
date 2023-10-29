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
}