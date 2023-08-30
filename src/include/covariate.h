
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <fr_matrix.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace arma;
#pragma once
class Covariate {
    public:
    std::vector<std::string> levels;
    bool standardize;
    std::string name;
    std::string type;
    std::string ref_level;

    Covariate(std::string cov, std::string cov_type, std::string cov_ref_level, std::string cov_levels, bool cov_standardize) {
        name = cov;
        type = cov_type;
        ref_level = cov_ref_level;
        standardize = cov_standardize;
        levels = split(cov_levels, ",");

        if(standardize && type != "numeric") {
            Rcpp::stop("Standardization of non-numeric variable not permitted");
        }
    }
    std::vector<std::string> split(std::string& val, std::string delim);
    void add_to_matrix(FRMatrix& df, FRMatrix& X_mat, double colinearity_rsq);
};
