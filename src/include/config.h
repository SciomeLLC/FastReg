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

#pragma once
class Config {
    public:
    Config(const std::string &file_path) {
        parse(file_path);
        validate_required_files();
        validate_keys();
        set_default_values();
        validate_args();
    }
    std::string pheno_file, POI_file_format, POI_file_delim, POI_effect_type, regression_type, pheno_file_delim,
        covar_file_delim, output_file_format, p_value_type, POI_file, covar_file, pheno_rowname_cols,
        covar_rowname_cols, phenotype, output_dir, POI_type, covariate_terms;

    bool no_intercept, verbose, compress_results, output_exclude_covar;
    double hwe_threshold, maf_threshold, colinearity_rsq;
    int max_iter, max_cores, poi_block_size;

    // optional
    std::string subject_subset_rowname_cols,
        subject_subset_file, subject_subset_delim,
        POI_subset_file, POI_subset_file_delim, POI_subset_rowname_col;

    std::vector<Covariate> covs;
    std::vector<std::string> covariates, covariate_levels, covariate_ref_level, covariate_type, split_by,
        covariate_standardize, POI_covar_interactions;
    private:
    std::unordered_map<std::string, std::string> values;
    template <typename T>
    T get(const std::string &key);
    bool has_key(const std::string &key);
    std::string get_value(std::string key, std::string def);
    void trim(std::string &s);
    void parse(std::string file_path, const char delim = '\t', const char comment = '#');
    void validate_keys();
    void validate_required_files();
    void set_default_values();
    void validate_args();
    std::vector<std::string> split(std::string val, std::string delim, std::string default_str_val, unsigned int size);

};
