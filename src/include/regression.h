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
        int max_iter,
        bool is_t_dist) = 0;

protected:
    static arma::fcolvec t_dist(arma::fcolvec abs_z, int df);
    static arma::fcolvec norm_dist(arma::fcolvec abs_z, int df);
};

class LogisticRegression : public RegressionBase
{

public:
    LogisticRegression()
    {
        arma::fmat tempM(250, 500);
        tempM.randu(); // = arma::fmat::randu(500, 10000);
        Eigen::MatrixXf mddata = Eigen::Map<Eigen::MatrixXf>(tempM.memptr(),
                                                             tempM.n_rows,
                                                             tempM.n_cols);
        ///////////////
        auto start = std::chrono::system_clock::now();
        tempM = tempM.t() * tempM;
        auto end = std::chrono::system_clock::now();
        ///////////////

        ///////////////
        auto start_a = std::chrono::system_clock::now();
        mddata = mddata.transpose() * mddata;
        auto end_a = std::chrono::system_clock::now();
        ///////////////
        const auto ms_int_a = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        const auto ms_int_b = std::chrono::duration_cast<std::chrono::nanoseconds>(end_a - start_a);
        double BLAS = std::chrono::duration<double>(ms_int_a).count();
        double EIGEN = std::chrono::duration<double>(ms_int_b).count();

        if (BLAS / EIGEN < 1)
        {
            USE_BLAS = true;
            Rcpp::Rcout << "Using loaded BLAS ilbrary." << std::endl;
        }
        else
        {
            USE_BLAS = false;
            Rcpp::Rcout << "Using RcppEigen as it appears faster than loaded BLAS library. " << std::endl;
        }
    }

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
        int max_iter,
        bool is_t_dist);

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
                   int max_iter,
                   bool is_t_dist);
    bool USE_BLAS;
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
        int max_iter,
        bool is_t_dist);
};
