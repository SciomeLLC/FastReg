#ifndef COVARIATE_H
#define COVARIATE_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <fr_matrix.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace arma;
class Covariate
{
public:
    std::vector<std::string> levels;
    bool standardize;
    std::string name;
    std::string cov_type;
    std::string ref_level;
    FRMatrix frmat;
    std::unordered_map<std::string, int> row_names_map;
    int col_idx;
    std::vector<std::string> col_names_arr;

    Covariate(std::string cov, std::string covar_type, std::string cov_ref_level, std::string cov_levels, bool cov_standardize)
    {
        name = cov;
        Rcpp::Rcout << cov << " has type " << covar_type << " with reference level " << cov_ref_level << " and levels " << cov_levels << std::endl;
        cov_type = covar_type;
        ref_level = cov_ref_level;
        standardize = cov_standardize;
        levels = split(cov_levels, ',');

        if (levels.empty() && cov_type != "numeric")
        {
            Rcpp::warning("Warning: Covariate %s is categorical but levels aren't specified. FastReg will generate the levels and select the first row as reference.", name);
        }

        if (standardize && cov_type != "numeric")
        {
            Rcpp::stop("Standardization of non-numeric variable not permitted");
        }
    }
    std::vector<std::string> split(std::string &val, char delim);
    void add_to_matrix(FRMatrix &df, FRMatrix &X_mat, double colinearity_rsq);
    void set_col_idx(int idx);
    void create_matrix(std::vector<std::vector<std::string>> tokenized, std::vector<std::string> row_names);
    FRMatrix filter_colinear(FRMatrix &design_mat, double colinearity_rsq);

private:
    FRMatrix create_categorical_matrix(std::vector<std::string> col_vals);
    FRMatrix create_numeric_matrix(std::vector<std::string> col_vals);
    void standardize_matrix(FRMatrix &cov_mat);
    void generate_levels(std::vector<std::string> col_vals);
};
#endif