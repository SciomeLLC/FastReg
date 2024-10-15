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
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "BEDMatrix.h"
#include <regression.h>
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

/// @brief Standardize (center/scale) X and remove features with low explainability; it is essentially PCA dimension reduction
/// @param x Matrix with zero'ed out columns for X/Y missing
/// @param dev_to_acc Amount of explained variance (summed Eigenvalues) to keep
/// @param N Number of non-zero rows (number of obs)
/// @return A clean X matrix for further calculations
arma::mat standardize_and_derank(arma::mat x,double dev_to_acc = 0.95, double N = 1.0){
  //scale (subtract mean and divide by SD) in place
  x.each_col([&N](arma::vec &y){
    y -= (arma::accu(y)/N);
    y /= sqrt(arma::accu(arma::square(y))/(N-1));
    });
  //get principal components
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::princomp(coeff, score, latent, x);
  //coeff = rotation in R, latent is the eigenvalues (sdev^2)
  // arma::uvec indices = sort_index(-1.0 * latent);
  // arma::vec sorted = sort(latent);
  arma::uvec ans = arma::find(cumsum(latent)/accu(latent) > dev_to_acc);
  return x * coeff.cols(0,ans[0]);
}


// [[Rcpp::export]]
void FastVLA_logistic(const arma::mat &Y, SEXP Gptr,
                          const arma::ivec &v_index, const arma::ivec &i_index,
                          const arma::mat &X, const int &chunk_size,
                          const std::string &dir,
                          const std::vector<std::string> cnames,
                          const std::vector<std::string> vnames,
                          const std::string suffix, const double epss = 1e-6,
                          const double &mafthresh = 0.005,
                          const double &pca_var_explained = 0.95) {

  int total_variants = v_index.n_elem;
  int num_chunks = total_variants / chunk_size;
  if (total_variants % chunk_size != 0) {
    num_chunks += 1;
  }

  // get boolean inclusion vectors
  arma::mat xy_bools(Y.n_rows, Y.n_cols);
  uvec badX = find_nan(arma::sum(X, 1.0));
  // put in each column booleans
  for (uword i = 0; i < Y.n_cols; i++) {
    // get X NANs for removal
    arma::vec x_bools(Y.n_rows, fill::ones);
    x_bools.elem(badX).fill(0.0);
    // add Y NANs for removal
    arma::uvec badY = find_nan(Y.col(i));
    x_bools.elem(badY).fill(0.0);
    xy_bools.col(i) = x_bools;
  }

#if !defined(__APPLE__) && !defined(__MACH__)
  omp_set_num_threads(1);
#endif
  // cnames = prefix per pheno - make into folder for each pheno
  // vnames = 1 vname per poi - aka rownames

  for (int i = 0; i < Y.n_cols; i++) {
    arma::mat cov = X;
    arma::ivec _index = i_index;
    arma::mat pheno = Y.col(i);
    arma::uvec nan_idx = arma::find(xy_bools.col(i) == 0.0);
    _index.shed_rows(nan_idx);
    pheno.shed_rows(nan_idx);
    cov.shed_rows(nan_idx);

    cov = standardize_and_derank(cov, pca_var_explained, double(cov.n_rows));

    arma::mat no_interactions = arma::mat(cov.n_rows, 1, arma::fill::ones);
    arma::mat interactions = arma::join_rows(no_interactions, cov);
    arma::mat interactions_sqrd = arma::join_rows(no_interactions, cov);
    for (int chunk = 0; chunk < num_chunks; ++chunk) {
      int start_idx = chunk * chunk_size;
      int end_idx = std::min(start_idx + chunk_size - 1, total_variants - 1);

      arma::ivec current_variant_indices = v_index.subvec(start_idx, end_idx);

      // Extract corresponding variant names
      std::vector<std::string> variant_names_chunk(
          vnames.begin() + start_idx, vnames.begin() + end_idx + 1);
      arma::mat G = scanBEDMatrix(Gptr, _index, current_variant_indices);
      VLAResult res(cov, G, no_interactions, interactions, interactions_sqrd);

      LogisticRegression regression;
      regression.run_vla(cov, pheno, G, res, 6, true, mafthresh);
      res.write_to_file(dir, suffix, cnames[i], variant_names_chunk);
    }
  }
}
