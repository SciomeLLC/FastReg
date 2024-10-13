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
/**
 * @brief The base class for all regression types.
 *
 * This class provides a virtual interface for running regressions with various data matrices.
 * It is meant to be extended by specific regression types such as `LogisticRegression` and `LinearRegression`.
 */
class RegressionBase
{
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

/**
 * @brief Logistic regression implementation.
 *
 * This class inherits from `RegressionBase` and implements logistic regression
 * using various matrices including covariates, phenotype data, and predictors of interest.
 */
class LogisticRegression : public RegressionBase
{

public:
  LogisticRegression() {}
  ~LogisticRegression() {}
  /**
   * @brief Runs logistic regression with the given data.
   *
   * This method performs logistic regression using covariates, phenotype data,
   * predictors of interest, and interactions, storing the resulting beta estimates, errors,
   * and p-values.
   *
   * @param cov Covariate matrix (independent variables).
   * @param pheno Phenotype matrix (dependent variable).
   * @param poi_data Predictors of interest data.
   * @param interactions Interaction matrix for the covariates and POI.
   * @param W2 Matrix indicating which rows to include.
   * @param beta_est Matrix to store beta estimates.
   * @param se_beta Matrix to store standard errors of beta estimates.
   * @param neglog10_pvl Matrix to store negative log10 p-values.
   * @param beta_rel_errs Vector to store relative errors of beta estimates.
   * @param beta_abs_errs Vector to store absolute errors of beta estimates.
   * @param iters Vector to store the number of iterations performed.
   * @param max_iter Maximum number of iterations for the regression.
   * @param x_mean Matrix of means for covariates.
   * @param x_sd Matrix of standard deviations for covariates.
   * @param xi_mean Matrix of means for interaction terms.
   * @param xi_sd Matrix of standard deviations for interaction terms.
   * @param is_t_dist Boolean indicating whether to use the t-distribution.
   */
  void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
           FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
           FRMatrix &se_beta, FRMatrix &neglog10_pvl,
           arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
           arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
           arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
           bool is_t_dist);

private:
  /**
   * @brief Runs logistic regression using BLAS operations.
   *
   * Performs the core logistic regression computation using Basic Linear Algebra Subprograms (BLAS).
   *
   * @param Same parameters as in the `run` method.
   */
  void run_BLAS(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
                FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
                arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
                arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
                bool is_t_dist);
};

/**
 * @brief Linear regression implementation.
 *
 * This class inherits from `RegressionBase` and implements linear regression
 * using various matrices including covariates, phenotype data, and predictors of interest.
 */
class LinearRegression : public RegressionBase
{
public:
  LinearRegression() {}

  ~LinearRegression() {}
  /**
   * @brief Runs linear regression with the given data.
   *
   * This method performs linear regression using covariates, phenotype data,
   * predictors of interest, and interactions, storing the resulting beta estimates, errors,
   * and p-values.
   *
   * @param cov Covariate matrix (independent variables).
   * @param pheno Phenotype matrix (dependent variable).
   * @param poi_data Predictors of interest data.
   * @param interactions Interaction matrix for the covariates and POI.
   * @param W2 Matrix indicating which rows to include.
   * @param beta_est Matrix to store beta estimates.
   * @param se_beta Matrix to store standard errors of beta estimates.
   * @param neglog10_pvl Matrix to store negative log10 p-values.
   * @param beta_rel_errs Vector to store relative errors of beta estimates.
   * @param beta_abs_errs Vector to store absolute errors of beta estimates.
   * @param iters Vector to store the number of iterations performed.
   * @param max_iter Maximum number of iterations for the regression.
   * @param x_mean Matrix of means for covariates.
   * @param x_sd Matrix of standard deviations for covariates.
   * @param xi_mean Matrix of means for interaction terms.
   * @param xi_sd Matrix of standard deviations for interaction terms.
   * @param is_t_dist Boolean indicating whether to use the t-distribution.
   */
  void run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
           FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
           FRMatrix &se_beta, FRMatrix &neglog10_pvl,
           arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
           arma::fcolvec &iters, int max_iter, arma::fmat &x_mean,
           arma::fmat &x_sd, arma::fmat &xi_mean, arma::fmat &xi_sd,
           bool is_t_dist);
};
#endif // REGRESSION_H