// utils.h
#ifndef UTILS_H
#define UTILS_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <fr_matrix.h>
#include <covariate.h>
#if !defined(__APPLE__) && !defined(__MACH__)
#include <omp.h>
#endif
#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#if __cplusplus >= 201703L && __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#endif
using namespace arma;

/**
 * @brief Function to check if an R interrupt has occurred.
 *
 * This function is called internally to check for user interrupts in R.
 */
static void chkIntFn(void *dummy);
/**
 * @brief Checks if the user has triggered an interrupt.
 *
 * This function checks for user interrupts and throws an exception if one is detected, stopping the current process.
 */
void checkInterrupt();
/**
 * @brief Checks if a directory exists.
 *
 * @param path The path to the directory.
 * @return true if the directory exists, false otherwise.
 */
bool dir_exists(const std::string &path);
/**
 * @brief Deletes a directory and its contents.
 *
 * @param path The path to the directory to be deleted.
 */
void delete_dir(const std::string &path);

/**
 * @brief Checks if a value is present in a container.
 *
 * @tparam T The type of the value.
 * @tparam U The type of the container.
 * @param value The value to search for.
 * @param container The container in which to search.
 * @return true if the value is found in the container, false otherwise.
 */
template <typename T, typename U>
bool isin(const T &value, const U &container);

/**
 * @brief Finds the intersection of two sets of row names.
 *
 * @param a The first vector of row names.
 * @param b The second vector of row names.
 * @return A vector containing the intersection of the two sets of row names.
 */
std::vector<std::string> intersect_row_names(const std::vector<std::string> &a, const std::vector<std::string> &b);

/**
 * @brief Computes the set difference between two sets of row names.
 *
 * @param a The first vector of row names.
 * @param b The second vector of row names.
 * @return A vector containing the elements of `a` that are not in `b`.
 */
std::vector<std::string> set_diff(const std::vector<std::string> &a, const std::vector<std::string> &b);

/**
 * @brief Transforms the Predictors of Interest (POI) matrix based on the specified effect type.
 *
 * This function modifies the POI matrix by applying transformations for dominant and recessive effect types.
 *
 * @param G The POI matrix to transform.
 * @param effect_type The type of effect to apply ("additive", "dominant", or "recessive").
 */
void transform_poi(FRMatrix &G, std::string effect_type = "additive");

/**
 * @brief Filters the POI matrix based on Minor Allele Frequency (MAF) and Hardy-Weinberg Equilibrium (HWE) thresholds.
 *
 * This function returns a matrix with filtered POI data and associated statistics such as MAF and HWE p-values.
 *
 * @param G The POI matrix to filter.
 * @param maf_threshold The MAF threshold.
 * @param hwe_threshold The HWE threshold.
 * @return A filtered FRMatrix containing the POI data and statistics.
 */
FRMatrix filter_poi(FRMatrix &G, double maf_threshold = 0.01, double hwe_threshold = 0.05);
/**
 * @brief Creates a matrix of interactions between POIs and covariates.
 *
 * @param df The FRMatrix containing the data.
 * @param poi_covar_interactions A vector of strings specifying the POI-covariate interactions.
 * @param Z The matrix where the interactions will be stored.
 */
void create_Z_matrix(FRMatrix &df, const std::vector<std::string> &poi_covar_interactions, FRMatrix &Z);

/**
 * @brief Returns the unique elements of a floating point vector.
 *
 * @tparam T The type of the elements.
 * @param vec The vector to process.
 * @return A vector containing the unique elements of the input vector.
 */
template <typename T>
std::vector<T> fr_unique(const frowvec &vec);

/**
 * @brief Creates a design matrix from covariates and filters based on collinearity.
 *
 * This function builds the design matrix from covariates and optionally adds an intercept term.
 *
 * @param df The FRMatrix containing the data.
 * @param covariates The vector of covariates to include in the design matrix.
 * @param no_intercept Whether to exclude the intercept term.
 * @param colinearity_rsq The collinearity threshold for filtering covariates.
 * @return An FRMatrix representing the design matrix.
 */
FRMatrix create_design_matrix(
    FRMatrix &df,
    std::vector<Covariate> covariates,
    bool no_intercept = false,
    double colinearity_rsq = 1.0);

/**
 * @brief Checks if a string contains only whitespace characters.
 *
 * @param str The string to check.
 * @return true if the string contains only whitespace, false otherwise.
 */
bool isWhitespace(const std::string &str);
#endif // UTILS_H
