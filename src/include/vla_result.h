#ifndef VLARESULT_H
#define VLARESULT_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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

class VLAResult {
public:
  arma::mat beta_est;
  arma::mat beta_est2;
  arma::mat beta_est2_sqrd;
  arma::mat se_beta;
  arma::mat se_beta2;
  arma::mat se_beta2_sqrd;
  arma::mat neglog10_pvl;
  arma::mat neglog10_pvl2;
  arma::mat neglog10_pvl2_sqrd;
  arma::mat *interactions = NULL;
  arma::mat *no_interactions = NULL;
  arma::mat *interactions_sqrd = NULL;

  arma::colvec beta_rel_errs;
  arma::colvec beta_abs_errs;
  arma::colvec beta_rel_errs2;
  arma::colvec beta_abs_errs2;

  arma::mat iters;
  arma::mat lls;
  arma::mat poi_sqrd_mat;
  arma::mat W2;

  std::vector<std::string> srt_cols;
  std::vector<std::string> srt_rows;
  std::vector<std::string> srt_rows2;

  std::vector<std::string> cov_int_names;
  std::vector<std::string> cov_int_names_sqrd;
  std::vector<std::string> cov_no_int_names;

  std::unordered_map<std::string, int> row_names;
  std::unordered_map<std::string, int> row_names2;
  std::unordered_map<std::string, int> col_names;
  std::unordered_map<std::string, int> col_names2;

  int num_parms, num_parms2, num_parms2_sqrd;

  VLAResult(){};
  VLAResult(arma::mat &covar_matrix, arma::mat &poi_matrix,
            arma::mat &no_interactions, arma::mat &interactions,
            arma::mat &interactions_sqrd);
  // VLAResult
  void set_lls(double ll1, double ll2, double lrs, double lrs_pval, int num_g,
               int idx, int rank);
  void set_betas_fit1(arma::colvec &beta, arma::colvec &se, arma::colvec &pval,
                      int idx);
  void set_betas_fit2(arma::colvec &beta, arma::colvec &se, arma::colvec &pval,
                      int idx);
  void set_betas_fit2_sqrd(arma::colvec &beta, arma::colvec &se,
                           arma::colvec &pval, int idx);
  void write_to_file(std::string dir, std::string file_name);
  static void concatenate(std::string output_dir, std::string file_name_prefix,
                          std::string file_concatenation_prefix);
};
#endif // VLARESULT_H