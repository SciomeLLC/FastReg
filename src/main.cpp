#define ARMA_WARN_LEVEL 0

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <regression.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "BEDMatrix.h"
#include <vla_result.h>

#if !defined(__APPLE__) && !defined(__MACH__)
#include <omp.h>
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/wait.h>
#endif

using namespace Rcpp;
using namespace arma;

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

#include <RcppEigen.h>

double recode_genotype2(int genotype) {
  double coding = arma::datum::nan; // missing
  if (genotype == 0) {
    coding = 2.0; // two copies of A1
  } else if (genotype == 3) {
    coding = 0.0; // zero copies of A1
  } else if (genotype == 2) {
    coding = 1.0; // one copy of A1
  }
  return coding;
}

//' Function to subset a BEDMatrix pointer and return an Integer Matrix
//' @param xptr denotes the slot(BEDMatrix, "xptr") from R
//' @param i denotes the row indices
//' @param j denotes the column indices
//' @return arma::fmat of variant calls
//' @export
//' @useDynLib FastVLA
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat scanBEDMatrix(SEXP xptr, arma::ivec i, arma::ivec j) {
  struct BEDMatrix *state =
      static_cast<struct BEDMatrix *>(R_ExternalPtrAddr(xptr));
  if (state == NULL) {
    Rf_error("BEDMatrix instance has been unmapped.");
  }
  int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
  int ni = i.n_rows;
  int nj = j.n_rows;
  arma::mat out(ni, nj);

  for (int cj = 0; cj < nj; cj++) {
    int jj = j[cj];
    for (int ci = 0; ci < ni; ci++) {
      int ii = i[ci];
      out(ci, cj) = recode_genotype2(extract_genotype_cartesian(
          state->data, ii - 1, jj - 1, num_bytes_per_variant));
    }
  }
  return out;
}
// [[Rcpp::export]]
int FastVLA_chunked_sota(const arma::mat &Y, SEXP Gptr,
                         const arma::ivec &v_index, const arma::ivec &i_index,
                         const arma::mat &X, const int &chunk_size,
                         const std::string &dir,
                         const std::vector<std::string> cnames,
                         const std::vector<std::string> vnames,
                         const std::string suffix, const double epss = 1e-6,
                         const double &mafthresh = 0.005) {
  arma::mat G = scanBEDMatrix(Gptr, v_index, i_index);
  arma::mat pheno = Y;
  arma::mat cov = X;
  arma::mat no_interactions = arma::mat(X.n_rows, 1, arma::fill::ones);
  arma::mat interactions = arma::join_rows(no_interactions, X);
  arma::mat interactions_sqrd = arma::join_rows(no_interactions, X);
  VLAResult res(cov, G, no_interactions, interactions, interactions_sqrd);
#if !defined(__APPLE__) && !defined(__MACH__)
  omp_set_num_threads(1);
#endif
  LogisticRegression regression;
  regression.run_vla(cov, pheno, G, res, 6, true);
  res.write_to_file(dir, suffix);
}
