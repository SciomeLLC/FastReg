#ifndef CONFIG_H
#define CONFIG_H
#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <list>
#include <algorithm>
#include <covariate.h>

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

class Config {
    public:
    
    std::string pheno_file, POI_file_format, POI_file_delim, POI_effect_type, regression_type, pheno_file_delim,
        covar_file_delim, output_file_format, p_value_type, POI_file_dir, covar_file, pheno_rowname_cols,
        covar_rowname_cols, phenotype, output_dir, POI_type, covariate_terms;

    bool no_intercept, verbose, compress_results, output_exclude_covar;
    double hwe_threshold, maf_threshold, colinearity_rsq, rel_conv_tolerance, abs_conv_tolerance;
    int max_iter, max_openmp_threads, poi_block_size, max_workers;

    // optional
    std::string subject_subset_rowname_cols,
        subject_subset_file, subject_subset_delim,
        POI_subset_file, POI_subset_file_delim, POI_subset_rowname_col;
    Rcpp::StringVector POI_covar_interactions_str, split_by_str;
    std::vector<Covariate> covs;
    std::vector<std::string> poi_files;
    std::vector<std::string> covariates, covariate_levels, covariate_ref_level, covariate_type, split_by,
        POI_covar_interactions;
    std::vector<bool> covariate_standardize;
    Config(
        std::string phenotype_str,
        std::string regression_type_str,
        std::string pvalue_dist_str,
        bool output_exclude_covar_bool,
        double maf_threshold_dbl,
        double hwe_threshold_dbl,
        bool no_intercept_bool,
        double colinearity_rsq_dbl,
        int poi_block_size_int,
        int max_iter_int,
        double rel_conv_tolerance_dbl,
        double abs_conv_tolderance_dbl, 
        int max_openmp_threads_int,
        std::string pheno_file_str,
        std::string pheno_rowname_cols_str,
        std::string pheno_file_delim_str,
        std::string covar_file_str,
        std::string covar_rowname_cols_str,
        std::string covar_file_delim_str,
        std::string poi_file_dir_str, 
        std::string poi_file_delim_str,
        std::string poi_file_format_str,
        std::string poi_type_str,
        std::string poi_effect_type_str,
        Rcpp::StringVector covariates_str,
        Rcpp::StringVector covariate_type_str,
        Rcpp::LogicalVector covariate_standardize_str,
        Rcpp::StringVector covariate_levels_str,
        Rcpp::StringVector covariate_ref_level_str,
        Rcpp::StringVector POI_covar_interactions_str,
        Rcpp::StringVector split_by_str,
        std::string output_dir_str,
        bool compress_results_bool,
        int max_workers_int
    ) {
        phenotype = phenotype_str;
        regression_type = regression_type_str;
        p_value_type = pvalue_dist_str;
        output_exclude_covar = output_exclude_covar_bool;
        maf_threshold = maf_threshold_dbl;
        hwe_threshold = hwe_threshold_dbl;
        no_intercept = no_intercept_bool;
        colinearity_rsq = colinearity_rsq_dbl;
        poi_block_size = poi_block_size_int;
        max_iter = max_iter_int;
        rel_conv_tolerance = rel_conv_tolerance_dbl;
        abs_conv_tolerance = abs_conv_tolderance_dbl;
        max_openmp_threads = max_openmp_threads_int;
        pheno_file = pheno_file_str;
        pheno_rowname_cols = pheno_rowname_cols_str;
        pheno_file_delim = pheno_file_delim_str;
        covar_file = covar_file_str;
        covar_file_delim = covar_file_delim_str;
        covar_rowname_cols = covar_rowname_cols_str;
        POI_file_dir = poi_file_dir_str;
        POI_effect_type = poi_effect_type_str;
        POI_file_delim = poi_file_delim_str;
        POI_file_format = poi_file_format_str;
        POI_type = poi_type_str;
        this->POI_covar_interactions_str = POI_covar_interactions_str;
        this->split_by_str = split_by_str;
        split_by = convert_stringV_to_string_arr(split_by_str);
        output_dir = output_dir_str;
        compress_results = compress_results_bool;
        max_workers = max_workers_int;
        validate_covariate_config(covariates_str, covariate_type_str, covariate_standardize_str, covariate_levels_str, covariate_ref_level_str, POI_covar_interactions_str);

        validate_required_files();
        validate_args();
        get_poi_files();
    }
    Config(){};
    void print();
    private:
    std::unordered_map<std::string, std::string> values;
    template <typename T>
    T get(const std::string &key);
    bool has_key(const std::string &key);
    std::string get_value(std::string key, std::string def);
    std::vector<std::string> convert_stringV_to_string_arr(Rcpp::StringVector stringV);
    void trim(std::string &s);
    void parse(std::string file_path, const char delim = '\t', const char comment = '#');
    void validate_keys();
    void validate_required_files();
    void set_default_values();
    void validate_args();
    void validate_covariate_config(
        Rcpp::StringVector covariates_str, 
        Rcpp::StringVector covariate_type_str, 
        Rcpp::LogicalVector covariate_standardize_str, 
        Rcpp::StringVector covariate_levels_str, 
        Rcpp::StringVector covariate_ref_level_str,
        Rcpp::StringVector POI_covar_interactions_str    
    );
    void get_poi_files();
    std::vector<std::string> split(std::string val, std::string delim, std::string default_str_val, unsigned int size);
};
#endif