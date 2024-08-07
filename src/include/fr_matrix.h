// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <names_map.h>
#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>
#pragma once

class FRMatrix {
public:
    FRMatrix() {};

    arma::fmat data;
    // NamesMap row_names;
    // NamesMap col_names;
    // NamesMap col_names_str;
    std::vector<std::vector<std::string>> str_data;
    std::unordered_map<std::string, int> col_names_str;
    
    std::unordered_map<std::string, int> row_names;
    std::unordered_map<std::string, int> col_names;
    FRMatrix(std::string filename, std::string delim, std::string &names) {
        load_from_csv(filename, delim, names);
        // if (!validate_cols(names)) {
        //     Rcpp::stop("Invalid rowname.cols for file: %s.", filename);
        // }
        Rcpp::Rcout << "Discovered " << row_names.size() << " unique subjects in " << filename << " file." << std::endl;
    };

    FRMatrix(std::string filename, std::string delim, std::string id, std::string phenotype) {
        load_from_csv(filename, delim, id, phenotype);
        // if (!validate_cols(names)) {
        //     Rcpp::stop("Invalid rowname.cols for file: %s.", filename);
        // }
        Rcpp::Rcout << "Discovered " << row_names.size() << " unique subjects in " << filename << " file." << std::endl;
    };

    FRMatrix(std::string filename, std::string delim, std::string id, std::vector<std::string> covariates, std::vector<std::string> cov_type) {
        load_from_csv(filename, delim, id, covariates, cov_type);
        // if (!validate_cols(names)) {
        //     Rcpp::stop("Invalid rowname.cols for file: %s.", filename);
        // }
        Rcpp::Rcout << "Discovered " << row_names.size() << " unique subjects in " << filename << " file." << std::endl;
    };

    bool validate_cols(std::string &names);
    void load_from_csv(std::string& filename, std::string& delim, std::string& id);
    void load_from_csv(std::string& filename, std::string& delim, std::string& id, std::string& phenotype);
    void load_from_csv(std::string& filename, std::string& delim, std::string& id, std::vector<std::string> covariates, std::vector<std::string> cov_type);
    FRMatrix get_submat_by_cols(const std::vector<int>& row_idx, const std::vector<std::string>& names);
    FRMatrix get_submat_by_cols(const std::vector<int>& row_idx, const std::unordered_map<std::string, int>& names);
    std::vector<std::string> split(const std::string& str_tokens, char delim);
    int get_col_idx(const std::string &col_name);
    int get_row_idx(const std::string &row_name);
    std::vector<std::string> get_col_str(const std::string &col_name);
    void write_summary(std::string dir, std::string name, int stratum, int process_id);
    bool file_exists(const std::string& name);
    std::vector<std::string> sort_map(bool rows);
    void print();
    static void write_results(
        FRMatrix& beta, 
        FRMatrix& se_beta, 
        FRMatrix& neglog10, 
        arma::umat& W2, 
        std::vector<std::string> poi_names, 
        std::string dir, 
        std::string file_name, 
        int stratum,
        bool exclude_covars,
        int process_id
    );
    static void write_convergence_results(
        FRMatrix& beta, 
        std::vector<std::string> poi_names, 
        std::string dir, 
        std::string file_name, 
        arma::fcolvec& rel_err,
        arma::fcolvec& abs_err,
        int stratum,
        int process_id
    ); 
    static void concatenate_results(
        std::string output_dir,
        std::string file_name_prefix,
        std::string file_concatenation_prefix
    );
    static void zip_results(std::string output_dir);
};
