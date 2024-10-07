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
#include <covariate.h>
#include <fr_matrix.h>

using namespace arma;

class RegressionBase {
public:
  virtual ~RegressionBase() {}
  virtual void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                   FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
                   FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                   arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
                   arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
                   arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
                   bool is_t_dist) = 0;
};

class LogisticRegression : public RegressionBase {

public:
  LogisticRegression() {}
  ~LogisticRegression() {}
  void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
           FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
           FRMatrix &se_beta, FRMatrix &neglog10_pvl,
           arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
           arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
           arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
           bool is_t_dist);

private:
  void run_BLAS(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
                FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
                arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
                arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
                bool is_t_dist);
};

class LinearRegression : public RegressionBase {
public:
  LinearRegression() {}

  ~LinearRegression() {}
  void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
           FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
           FRMatrix &se_beta, FRMatrix &neglog10_pvl,
           arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
           arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
           arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
           bool is_t_dist);
};
#endif // REGRESSION_H