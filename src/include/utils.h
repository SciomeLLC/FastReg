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
#  if __cplusplus >= 201703L && __has_include(<filesystem>)
#    include <filesystem>
namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#  endif
#endif
using namespace arma;

static void chkIntFn(void *dummy);
void checkInterrupt();
bool dir_exists(const std::string& path);

void delete_dir(const std::string& path);

template<typename T, typename U>
bool isin(const T& value, const U& container);

std::vector<std::string> intersect_row_names(const std::vector<std::string>& a, const std::vector<std::string>& b);

std::vector<std::string> set_diff(const std::vector<std::string>& a, const std::vector<std::string>& b);

void transform_poi(FRMatrix &G, std::string effect_type = "additive");

FRMatrix filter_poi(FRMatrix &G, double maf_threshold = 0.01, double hwe_threshold = 0.05);
void create_Z_matrix(FRMatrix& df, const std::vector<std::string>& poi_covar_interactions, FRMatrix& Z);


template <typename T>
std::vector<T> fr_unique(const frowvec &vec);


FRMatrix create_design_matrix(
    FRMatrix &df,
    std::vector<Covariate> covariates,
    bool no_intercept=false,
    double colinearity_rsq=1.0
);

/// @brief given a string it determines if it is only whitespace
/// @param str the string in question
/// @return true if only whitespace, false if true
bool isWhitespace(const std::string &str);
#endif // UTILS_H
