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
/**
 * @brief The Covariate class manages individual covariates in a regression model.
 *
 * It handles both numeric and categorical covariates, including creating matrices,
 * generating levels for categorical variables, and standardizing numeric variables.
 */
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
    /**
     * @brief Constructs a Covariate object.
     *
     * @param cov Name of the covariate.
     * @param covar_type Type of the covariate (e.g., "numeric", "categorical").
     * @param cov_ref_level Reference level for categorical covariates.
     * @param cov_levels Comma-separated string of levels for the covariate.
     * @param cov_standardize Boolean flag to indicate if standardization is required.
     */
    Covariate(std::string cov, std::string covar_type, std::string cov_ref_level, std::string cov_levels, bool cov_standardize)
    {
        name = cov;
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
    /**
     * @brief Splits a string by a given delimiter.
     *
     * @param val The string to split.
     * @param delim The delimiter character.
     * @return A vector of strings after splitting.
     */
    std::vector<std::string> split(std::string &val, char delim);
    /**
     * @brief Adds the covariate to the design matrix.
     *
     * @param df The FRMatrix containing the design matrix.
     * @param X_mat The resulting FRMatrix to which the covariate will be added.
     * @param colinearity_rsq The threshold for collinearity filtering.
     */
    void add_to_matrix(FRMatrix &df, FRMatrix &X_mat, double colinearity_rsq);
    /**
     * @brief Sets the column index of the covariate.
     *
     * @param idx The index of the covariate column.
     */
    void set_col_idx(int idx);
    /**
     * @brief Prints the covariate details to the console.
     */
    void print();
    /**
     * @brief Creates a matrix representation for the covariate.
     *
     * @param values The values for the covariate from the dataset.
     * @param row_names The row names for the covariate data.
     */
    void create_matrix(std::vector<std::vector<std::string>> values, std::vector<std::string> row_names);
    /**
     * @brief Filters collinear covariate columns based on the specified R-squared threshold.
     *
     * @param design_mat The design matrix to check for collinearity.
     * @param colinearity_rsq The R-squared threshold for collinearity filtering.
     * @return An FRMatrix containing non-collinear covariate columns.
     */
    FRMatrix filter_colinear(FRMatrix &design_mat, double colinearity_rsq);

private:
    /**
     * @brief Creates a matrix for a categorical covariate.
     *
     * @param col_vals The column values of the categorical covariate.
     * @return A matrix representation of the categorical covariate.
     */
    FRMatrix create_categorical_matrix(std::vector<std::string> col_vals);
    /**
     * @brief Creates a matrix for a numeric covariate.
     *
     * @param col_vals The column values of the numeric covariate.
     * @return A matrix representation of the numeric covariate.
     */
    FRMatrix create_numeric_matrix(std::vector<std::string> col_vals);
    /**
     * @brief Standardizes the covariate matrix (mean 0, stddev 1).
     *
     * @param cov_mat The matrix to be standardized.
     */
    void standardize_matrix(FRMatrix &cov_mat);
    /**
     * @brief Generates unique levels for a categorical covariate if none are provided.
     *
     * @param col_vals The column values for the categorical covariate.
     */
    void generate_levels(std::vector<std::string> col_vals);
};
#endif