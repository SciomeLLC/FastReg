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
#include <vla_result.h>

using namespace arma;

class RegressionBase {
public:
  virtual ~RegressionBase() {}

  virtual void run_vla(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
                       VLAResult &result, int max_iter, bool is_t_dist) = 0;
};

class LogisticRegression : public RegressionBase {

public:
  LogisticRegression() {}
  ~LogisticRegression() {}
  void run_vla(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
               VLAResult &result, int max_iter, bool is_t_dist);

};

class LinearRegression : public RegressionBase {
public:
  LinearRegression() {}

  ~LinearRegression() {}
  void run_vla(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
               VLAResult &result, int max_iter, bool is_t_dist) {
    return;
  }
};
#endif // REGRESSION_H