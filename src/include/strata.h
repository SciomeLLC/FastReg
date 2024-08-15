
#ifndef STRATA_H
#define STRATA_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <fr_matrix.h>

using namespace arma;

class Strata {
public:
    Strata() {};
    int nstrata;
    std::vector<std::string> ids;
    std::unordered_map<std::string, std::vector<std::string>> index_list;
    void stratify(const std::vector<std::string>& split_by, FRMatrix &covar_df, const std::vector<std::string>& common_individuals);
};
#endif