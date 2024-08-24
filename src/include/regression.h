// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <chrono>
#define R_NO_REMAP
#include <Rmath.h>
#include <fr_matrix.h>
#include <covariate.h>
#include <stats.hpp>
#include <RcppEigen.h>
#pragma once

using namespace arma;

class RegressionBase
{
public:
    virtual ~RegressionBase() {}
    virtual void run(
        FRMatrix &cov,
        FRMatrix &pheno,
        FRMatrix &poi_data,
        FRMatrix &interactions,
        arma::umat &W2,
        FRMatrix &beta_est,
        FRMatrix &se_beta,
        FRMatrix &neglog10_pvl,
        arma::fcolvec &beta_rel_errs,
        arma::fcolvec &beta_abs_errs,
        arma::fcolvec &iters,
        int max_iter,
        bool is_t_dist, 
        bool use_blas) = 0;
};

class LogisticRegression : public RegressionBase
{

public:
    LogisticRegression(){}

    void run(
        FRMatrix &cov,
        FRMatrix &pheno,
        FRMatrix &poi_data,
        FRMatrix &interactions,
        arma::umat &W2,
        FRMatrix &beta_est,
        FRMatrix &se_beta,
        FRMatrix &neglog10_pvl,
        arma::fcolvec &beta_rel_errs,
        arma::fcolvec &beta_abs_errs,
        arma::fcolvec &iters,
        int max_iter,
        bool is_t_dist,
        bool use_blas);

private:
    void run_BLAS(
        FRMatrix &cov,
        FRMatrix &pheno,
        FRMatrix &poi_data,
        FRMatrix &interactions,
        arma::umat &W2,
        FRMatrix &beta_est,
        FRMatrix &se_beta,
        FRMatrix &neglog10_pvl,
        arma::fcolvec &beta_rel_errs,
        arma::fcolvec &beta_abs_errs,
        arma::fcolvec &iters,
        int max_iter,
        bool is_t_dist);
    void run_EIGEN(FRMatrix &cov,
                   FRMatrix &pheno,
                   FRMatrix &poi_data,
                   FRMatrix &interactions,
                   arma::umat &W2,
                   FRMatrix &beta_est,
                   FRMatrix &se_beta,
                   FRMatrix &neglog10_pvl,
                   arma::fcolvec &beta_rel_errs,
                   arma::fcolvec &beta_abs_errs,
                   arma::fcolvec &iters,
                   int max_iter,
                   bool is_t_dist);
};

class LinearRegression : public RegressionBase
{
public:
    void run(
        FRMatrix &cov,
        FRMatrix &pheno,
        FRMatrix &poi_data,
        FRMatrix &interactions,
        arma::umat &W2,
        FRMatrix &beta_est,
        FRMatrix &se_beta,
        FRMatrix &neglog10_pvl,
        arma::fcolvec &beta_rel_errs,
        arma::fcolvec &beta_abs_errs,
        arma::fcolvec &iters,
        int max_iter,
        bool is_t_dist,
        bool use_blas);
};
