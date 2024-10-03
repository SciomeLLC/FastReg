#ifndef FRRESULT_H
#define FRRESULT_H
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

class FRMatrix;
class FRResult
{
public:
  arma::fmat beta_est;
  arma::fmat beta_est2;
  arma::fmat beta_est2_sqrd;
  arma::fmat se_beta;
  arma::fmat se_beta2;
  arma::fmat se_beta2_sqrd;
  arma::fmat neglog10_pvl;
  arma::fmat neglog10_pvl2;
  arma::fmat neglog10_pvl2_sqrd;
  FRMatrix *interactions = NULL;
  FRMatrix *no_interactions = NULL;
  FRMatrix *interactions_sqrd = NULL;

  arma::fcolvec beta_rel_errs;
  arma::fcolvec beta_abs_errs;
  arma::fcolvec beta_rel_errs2;
  arma::fcolvec beta_abs_errs2;

  arma::fmat iters;
  arma::fmat lls;
  arma::fmat poi_sqrd_mat;
  arma::fmat W2;

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

  FRResult() {};
  FRResult(FRMatrix &covar_matrix, FRMatrix &poi_matrix,
           FRMatrix &no_interactions, FRMatrix &interactions, FRMatrix &interactions_sqrd);

  void set_lls(double ll1, double ll2, double lrs, double lrs_pval, int num_g,
               int idx);
  void set_betas_fit1(arma::fcolvec &beta, arma::fcolvec &se,
                      arma::fcolvec &pval, int idx);
  void set_betas_fit2(arma::fcolvec &beta, arma::fcolvec &se,
                      arma::fcolvec &pval, int idx);
  void write_to_file(std::string dir, std::string file_name, int stratum, int process_id);
  static void concatenate(std::string output_dir,
                          std::string file_name_prefix,
                          std::string file_concatenation_prefix);
};
#endif // FRRESULT_H