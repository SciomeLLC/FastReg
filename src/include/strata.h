
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
/**
 * @brief The Strata class is used to group individuals into strata based on specified covariates.
 *
 * This class allows for the stratification of a dataset by grouping individuals according
 * to unique combinations of values in specified covariate columns.
 */
class Strata
{
public:
    Strata() {};
    int nstrata;
    std::vector<std::string> ids;
    std::unordered_map<std::string, std::vector<std::string>> index_list;
    /**
     * @brief Stratifies individuals based on specified covariate columns.
     *
     * This method groups individuals into strata based on the unique combinations of values in the specified covariate columns.
     * If no columns are provided, all individuals are placed in a single stratum.
     *
     * @param split_by A vector of covariate column names by which to stratify the individuals.
     * @param covar_df The FRMatrix object containing the covariate data.
     * @param common_individuals A vector of individual identifiers present in the dataset.
     */
    void stratify(const std::vector<std::string> &split_by, FRMatrix &covar_df, const std::vector<std::string> &common_individuals);
};
#endif