#ifndef REGRESSION_H
#define REGRESSION_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <iterator>
#include <string>
#include <vector>
#define R_NO_REMAP
#include <RcppEigen.h>
#include <Rmath.h>
#include <utils.h>
#include <vla_result.h>
#if !defined(__APPLE__) && !defined(__MACH__)
#include <omp.h>
#endif

using namespace arma;
class LogisticRegression {

public:
  LogisticRegression() {}
  ~LogisticRegression() {}
  void run_vla(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
               VLAResult &result, const int max_iter, bool is_t_dist,
               double maf_thresh, double epss);
  void run_vla(arma::fmat &cov, arma::fmat &pheno, arma::fmat &poi_data,
               VLAResultf &result, const int max_iter, bool is_t_dist,
               double maf_thresh, float epss);
  void run_vla_2(arma::fmat &cov, arma::fmat &pheno, arma::fmat &poi_data,
                 VLAResultf &result, const int max_iter, bool is_t_dist,
                 std::vector<int> &poi_2_idx, Eigen::MatrixXf &W2f,
                 Eigen::MatrixXf &tphenoD, float epss);
  void run_vla_3(arma::fmat &cov, arma::fmat &pheno, arma::fmat &poi_data,
                 VLAResultf &result, const int max_iter, bool is_t_dist,
                 std::vector<int> &poi_3_idx, Eigen::MatrixXf &W2f,
                 Eigen::MatrixXf &tphenoD, float epss);
  void run_vla_2(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
                 VLAResult &result, const int max_iter, bool is_t_dist,
                 std::vector<int> &poi_2_idx, Eigen::MatrixXd &W2f,
                 Eigen::MatrixXd &tphenoD, double epss);
  void run_vla_3(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
                 VLAResult &result, const int max_iter, bool is_t_dist,
                 std::vector<int> &poi_3_idx, Eigen::MatrixXd &W2f,
                 Eigen::MatrixXd &tphenoD, double epss);
  void checkInterrupt() {
    if (R_ToplevelExec(chkIntFn, NULL) == FALSE) {
      Rcpp::stop("Received user interrupt. Stopping FastReg...");
    }
  }
};
#endif // REGRESSION_H