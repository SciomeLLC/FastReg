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
#include <stats.hpp>
#pragma once

using namespace arma;

class RegressionBase {
public:
  virtual ~RegressionBase() {}
  virtual void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                   FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
                   FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                   arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
                   arma::fcolvec &iters, int max_iter, bool is_t_dist,
                   bool use_blas) = 0;
  virtual void run_vla(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                       arma::fmat &poi_sqrd, FRMatrix &interactions,
                       arma::umat &W2, FRMatrix &beta_est, FRMatrix &se_beta,
                       FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
                       arma::fcolvec &beta_abs_errs, FRMatrix &beta_est2,
                       FRMatrix &se_beta2, FRMatrix &neglog10_pvl2,
                       arma::fcolvec &beta_rel_errs2,
                       arma::fcolvec &beta_abs_errs2, arma::fcolvec &iters,
                       arma::fmat &lls, int max_iter, bool is_t_dist);
};

class LogisticRegression : public RegressionBase {

public:
  LogisticRegression() {}

  void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
           FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
           FRMatrix &se_beta, FRMatrix &neglog10_pvl,
           arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
           arma::fcolvec &iters, int max_iter, bool is_t_dist, bool use_blas);

  void run_vla(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
               arma::fmat &poi_sqrd, FRMatrix &interactions, arma::umat &W2,
               FRMatrix &beta_est, FRMatrix &se_beta, FRMatrix &neglog10_pvl,
               arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
               FRMatrix &beta_est2, FRMatrix &se_beta2, FRMatrix &neglog10_pvl2,
               arma::fcolvec &beta_rel_errs2, arma::fcolvec &beta_abs_errs2,
               arma::fcolvec &iters, arma::fmat &lls, int max_iter,
               bool is_t_dist);

private:
  void run_BLAS(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
                FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
                arma::fcolvec &iters, int max_iter, bool is_t_dist);
  void run_EIGEN(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                 FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
                 FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                 arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
                 arma::fcolvec &iters, int max_iter, bool is_t_dist);
};

class LinearRegression : public RegressionBase {
public:
  void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
           FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
           FRMatrix &se_beta, FRMatrix &neglog10_pvl,
           arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
           arma::fcolvec &iters, int max_iter, bool is_t_dist, bool use_blas);
  void run_vla(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                       arma::fmat &poi_sqrd, FRMatrix &interactions,
                       arma::umat &W2, FRMatrix &beta_est, FRMatrix &se_beta,
                       FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
                       arma::fcolvec &beta_abs_errs, FRMatrix &beta_est2,
                       FRMatrix &se_beta2, FRMatrix &neglog10_pvl2,
                       arma::fcolvec &beta_rel_errs2,
                       arma::fcolvec &beta_abs_errs2, arma::fcolvec &iters,
                       arma::fmat &lls, int max_iter, bool is_t_dist) {
                        return;
                       }
};
