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
#pragma once

using namespace arma;

class RegressionBase {
    public:
        virtual ~RegressionBase() {}
        virtual void run(
            FRMatrix& cov, 
            FRMatrix& pheno, 
            FRMatrix& poi_data,
            FRMatrix& interactions,
            arma::umat& W2,
            FRMatrix& beta_est,
            FRMatrix& se_beta,
            FRMatrix& neglog10_pvl,
            arma::colvec& beta_rel_errs,
            arma::colvec& beta_abs_errs,
            int max_iter, 
            bool is_t_dist
        ) = 0;
    protected:
        static arma::colvec t_dist(arma::colvec abs_z, int df);
        static arma::colvec norm_dist(arma::colvec abs_z, int df);
};

class LogisticRegression : public RegressionBase {
    public:
        void run(
            FRMatrix& cov, 
            FRMatrix& pheno, 
            FRMatrix& poi_data,
            FRMatrix& interactions,
            arma::umat& W2,
            FRMatrix& beta_est,
            FRMatrix& se_beta,
            FRMatrix& neglog10_pvl,
            arma::colvec& beta_rel_errs,
            arma::colvec& beta_abs_errs,
            int max_iter, 
            bool is_t_dist
        );
};

class LinearRegression : public RegressionBase {
    public:
        void run(
            FRMatrix& cov, 
            FRMatrix& pheno, 
            FRMatrix& poi_data,
            FRMatrix& interactions,
            arma::umat& W2,
            FRMatrix& beta_est,
            FRMatrix& se_beta,
            FRMatrix& neglog10_pvl,
            arma::colvec& beta_rel_errs,
            arma::colvec& beta_abs_errs,
            int max_iter, 
            bool is_t_dist
        );
};
