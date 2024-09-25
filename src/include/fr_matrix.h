#ifndef FRMATRIX_H
#define FRMATRIX_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <fstream>
#include <mutex>
#include <names_map.h>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>

class FRMatrix {
public:
  FRMatrix(){};

  arma::fmat data;
  std::vector<std::string> col_name_str_arr;
  std::vector<std::vector<std::string>> str_data;
  std::unordered_map<std::string, int> col_names_str;
  std::vector<std::string> col_names_arr;
  std::vector<std::string> row_names_arr;
  std::unordered_map<std::string, int> row_names;
  std::unordered_map<std::string, int> col_names;

  bool validate_cols(std::string &names);
  FRMatrix get_submat_by_cols(const std::vector<int> &row_idx,
                              const std::vector<std::string> &names);
  FRMatrix
  get_submat_by_cols(const std::vector<int> &row_idx,
                     const std::unordered_map<std::string, int> &names);
  std::vector<std::string> split(const std::string &str_tokens, char delim);
  int get_col_idx(const std::string &col_name);
  int get_row_idx(const std::string &row_name);
  std::vector<std::string> get_col_str(const std::string &col_name);
  void write_summary(std::string dir, std::string name, int stratum,
                     int process_id);
  bool file_exists(const std::string &name);
  std::vector<std::string> sort_map(bool rows);
  void join(FRMatrix &frmat);
  void print();
  static void write_results(FRMatrix &beta, FRMatrix &se_beta,
                            FRMatrix &neglog10, arma::umat &W2,
                            arma::fcolvec &rel_err, arma::fcolvec &abs_err,
                            arma::fcolvec &iters,
                            std::vector<std::string> poi_names, std::string dir,
                            std::string file_name, int stratum,
                            bool exclude_covars, int process_id);
  static void write_vla_results(
      FRMatrix &beta, FRMatrix &se_beta, FRMatrix &neglog10, arma::umat &W2,
      arma::fcolvec &rel_err, arma::fcolvec &abs_err, FRMatrix &beta2,
      FRMatrix &se_beta2, FRMatrix &neglog102, arma::fcolvec &rel_err2,
      arma::fcolvec &abs_err2, arma::fmat &lls, arma::fcolvec &iters,
      std::vector<std::string> poi_names, std::string dir,
      std::string file_name, int stratum, bool exclude_covars, int process_id);
  static void write_convergence_results(FRMatrix &beta,
                                        std::vector<std::string> poi_names,
                                        std::string dir, std::string file_name,
                                        arma::fcolvec &rel_err,
                                        arma::fcolvec &abs_err, int stratum,
                                        int process_id);
  static void concatenate_results(std::string output_dir,
                                  std::string file_name_prefix,
                                  std::string file_concatenation_prefix);
  static void zip_results(std::string output_dir);
};
#endif