// [[Rcpp::depends(RcppArmadillo)]]
#include <regression.h>
#include <utils.h>
#include <iostream>
#include <chrono>

#include <RcppEigen.h>

template <typename T>
Rcpp::NumericVector arma2vec(const T &x)
{
  return Rcpp::NumericVector(x.begin(), x.end());
}

Rcpp::NumericVector t_dist_r(arma::fcolvec abs_z, int df)
{
  return Rcpp::dt(arma2vec(abs_z), df, true);
}

Rcpp::NumericVector norm_dist_r(arma::fcolvec abs_z, int df)
{
  return Rcpp::dnorm(arma2vec(abs_z), 1.0, 1.0, true);
}

typedef Eigen::Matrix<arma::uword, Eigen::Dynamic, Eigen::Dynamic> MatrixUd;
arma::fcolvec RegressionBase::t_dist(arma::fcolvec abs_z, int df)
{

  arma::fcolvec pvalues = conv_to<fcolvec>::from(
      -1 * (stats2::pt(abs_z, df, true) + log(2)) / log(10));
  return pvalues;
}

arma::fcolvec RegressionBase::norm_dist(arma::fcolvec abs_z, int df)
{
  return conv_to<fcolvec>::from(
      -1 * (stats2::pnorm(abs_z, 1.0, 1.0, true) + log(2)) / log(10));
}

arma::fcolvec pmean(arma::fcolvec a, arma::fcolvec b) { return (a + b) / 2; }

void LogisticRegression::run_BLAS(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                                  FRMatrix &interactions, arma::umat &W2,
                                  FRMatrix &beta_est, FRMatrix &se_beta,
                                  FRMatrix &neglog10_pvl,
                                  arma::fcolvec &beta_rel_errs,
                                  arma::fcolvec &beta_abs_errs, int max_iter,
                                  bool is_t_dist)
{

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  // create a pointer to the specified distribution function
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist : norm_dist;
#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++)
  {
    checkInterrupt();
    arma::fmat A(n_parms, n_parms, arma::fill::zeros);
    // Initialize beta
    arma::fcolvec beta(n_parms, arma::fill::zeros);
    arma::fcolvec beta_old(n_parms, arma::fill::ones);
    arma::fmat cov_w_mat = cov.data;
    arma::fmat int_w_mat = interactions.data;
    arma::ucolvec w2_col = W2.col(poi_col);
    int_w_mat = interactions.data % arma::repmat(poi_data.data.col(poi_col), 1,
                                                 interactions.data.n_cols);

    arma::uword first_chunk = cov_w_mat.n_cols - 1;
    arma::uword second_chunk = n_parms - 1;
    arma::span first = arma::span(0, first_chunk);
    arma::span second = arma::span(first_chunk + 1, second_chunk);
    arma::fcolvec beta_diff = arma::abs(beta - beta_old);
    // arma::span col_1 = arma::span(0,0);
    beta_rel_errs.at(poi_col) = 1;

    for (int iter = 0; iter < max_iter && beta_rel_errs.at(poi_col) > 1e-4; iter++)
    {

      arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
      eta += cov_w_mat * beta.subvec(first);
      eta += int_w_mat * beta.subvec(second);

      arma::fmat p = 1 / (1 + arma::exp(-eta));
      arma::fmat W1 = p % (1 - p) % w2_col;
      arma::fmat temp1 = cov_w_mat % arma::repmat(W1, 1, cov_w_mat.n_cols);
      arma::fmat temp2 = int_w_mat % arma::repmat(W1, 1, int_w_mat.n_cols);

      A.submat(first, first) = temp1.t() * cov_w_mat;        // C'C
      A.submat(second, second) = temp2.t() * int_w_mat;      // I'I
      A.submat(first, second) = temp1.t() * int_w_mat;       // C'I
      A.submat(second, first) = A.submat(first, second).t(); // I'C

      arma::fmat z = w2_col % (pheno.data - p);
      arma::fmat B = arma::fmat(n_parms, 1, arma::fill::zeros);
      B.submat(first, arma::span(0, 0)) = cov_w_mat.t() * z;
      B.submat(second, arma::span(0, 0)) = int_w_mat.t() * z;

      beta = beta + arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = arma::abs(beta - beta_old);
      beta_abs_errs.at(poi_col) = beta_diff.max();
      beta_rel_errs.at(poi_col) = (beta_diff / arma::abs(beta_old)).max();
      beta_old = beta;
    }

    if (beta_abs_errs.at(poi_col) > 1e-4)
    { // update the covariance matrix if there
      // was a maximum relative change in the
      // betas greater than 1e-4
      arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
      eta += cov_w_mat * beta.subvec(first);
      eta += int_w_mat * beta.subvec(second);

      arma::fmat p = 1 / (1 + arma::exp(-eta));
      arma::fmat W1 = p % (1 - p) % w2_col;
      arma::fmat temp1 = cov_w_mat % arma::repmat(W1, 1, cov_w_mat.n_cols);
      arma::fmat temp2 = int_w_mat % arma::repmat(W1, 1, int_w_mat.n_cols);

      A.submat(first, first) = temp1.t() * cov_w_mat;        // C'C
      A.submat(second, second) = temp2.t() * int_w_mat;      // I'I
      A.submat(first, second) = temp1.t() * int_w_mat;       // C'I
      A.submat(second, first) = A.submat(first, second).t(); // I'C
    }

    beta_abs_errs.at(poi_col) = beta_diff.max();
    beta_rel_errs.at(poi_col) = (beta_diff / arma::abs(beta)).max();
    int df = arma::as_scalar(arma::sum(w2_col, 0)) - n_parms;

    // arma::fcolvec temp_se = arma::sqrt(arma::diagvec(arma::inv_sympd(A)));
    arma::fcolvec temp_se = arma::sqrt(arma::diagvec(arma::pinv(A)));
    beta_est.data.col(poi_col) = beta;
    se_beta.data.col(poi_col) = temp_se;
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}

void LogisticRegression::run_EIGEN(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                                   FRMatrix &interactions, arma::umat &W2,
                                   FRMatrix &beta_est, FRMatrix &se_beta,
                                   FRMatrix &neglog10_pvl,
                                   arma::fcolvec &beta_rel_errs,
                                   arma::fcolvec &beta_abs_errs, int max_iter,
                                   bool is_t_dist)
{

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  // create a pointer to the specified distribution function
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist : norm_dist;
  MatrixUd W2_EIG = Eigen::Map<MatrixUd>(W2.memptr(),
                                         W2.n_rows,
                                         W2.n_cols);
  Eigen::MatrixXf W2f = W2_EIG.cast<float>();
  Eigen::MatrixXf tphenoD = Eigen::Map<Eigen::MatrixXf>(pheno.data.memptr(),
                                                        pheno.data.n_rows,
                                                        pheno.data.n_cols);

#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++)
  {
    checkInterrupt();
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(n_parms, n_parms);
    A(0, 0) = 1;
    A(n_parms - 1, n_parms - 1) = 10;
    // Initialize beta

    Eigen::VectorXf beta = Eigen::VectorXf::Zero(n_parms, 1);
    Eigen::VectorXf beta_old = beta;

    Eigen::MatrixXf cov_w_mat = Eigen::Map<Eigen::MatrixXf>(cov.data.memptr(),
                                                            cov.data.n_rows,
                                                            cov.data.n_cols);
    Eigen::MatrixXf int_w_mat = Eigen::Map<Eigen::MatrixXf>(interactions.data.memptr(),
                                                            interactions.data.n_rows,
                                                            interactions.data.n_cols);

    Eigen::VectorXf w2_col = W2f.col(poi_col);
    arma::fmat POI_ARMA = poi_data.data.col(poi_col);
    Eigen::MatrixXf POI = Eigen::Map<Eigen::MatrixXf>(POI_ARMA.memptr(),
                                                      POI_ARMA.n_rows,
                                                      1);
    Eigen::MatrixXf temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXf X = Eigen::MatrixXf::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXf beta_diff = (beta - beta_old).array().abs();
    beta_rel_errs.at(poi_col) = 1e8;

    for (int iter = 0; iter < max_iter && beta_rel_errs.at(poi_col) > 1e-4; iter++)
    {
      Eigen::MatrixXf eta = Eigen::MatrixXf::Zero(cov.data.n_rows, 1);
      // arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
      eta = X * beta;
      Eigen::MatrixXf p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();

      Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;

      Eigen::MatrixXf z = w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXf B = Eigen::MatrixXf::Zero(n_parms, 1);

      B = X.transpose() * z;

      beta = beta + A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      beta_abs_errs.at(poi_col) = beta_diff.array().maxCoeff();
      Eigen::MatrixXf mTemp = (beta_diff.array() / beta_old.array());
      beta_rel_errs.at(poi_col) = mTemp.maxCoeff();
      beta_old = beta;
    }

    if (beta_abs_errs.at(poi_col) > 1e-4)
    { // update the covariance matrix if there
      // was a maximum relative change in the
      // betas greater than 1e-4
      Eigen::MatrixXf p = -X * beta;
      p = p.array().exp() + 1;
      p = 1 / p.array();

      Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;
    }

    beta_abs_errs.at(poi_col) = beta_diff.maxCoeff();
    beta_rel_errs.at(poi_col) = (beta_diff.array() / beta.array().abs()).maxCoeff();
    int df = w2_col.array().sum() - n_parms;
    Eigen::MatrixXf temp_inv = A.inverse();
    Eigen::VectorXf diag = temp_inv.diagonal().array().sqrt();

    // convert it all back to armadillo matrix types
    arma::fcolvec temp_se = arma::fcolvec(diag.data(), diag.size(),
                                          true, false);
    arma::fcolvec temp_b = arma::fcolvec(beta.data(), beta.size(),
                                         true, false);
    beta_est.data.col(poi_col) = temp_b;
    se_beta.data.col(poi_col) = temp_se;
    arma::fcolvec neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}

void LogisticRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                             FRMatrix &interactions, arma::umat &W2,
                             FRMatrix &beta_est, FRMatrix &se_beta,
                             FRMatrix &neglog10_pvl,
                             arma::fcolvec &beta_rel_errs,
                             arma::fcolvec &beta_abs_errs, int max_iter,
                             bool is_t_dist)
{
  if (USE_BLAS)
  { // use BLAS if faster than Eigen::MatrixXf calculations run_BLAS
    run_BLAS(cov, pheno, poi_data,
             interactions, W2, beta_est, se_beta, neglog10_pvl,
             beta_rel_errs, beta_abs_errs, max_iter, is_t_dist);
  }
  else
  {
    run_EIGEN(cov, pheno, poi_data,
              interactions, W2, beta_est, se_beta, neglog10_pvl,
              beta_rel_errs, beta_abs_errs, max_iter, is_t_dist);
  }
}

void LinearRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                           FRMatrix &interactions, arma::umat &W2,
                           FRMatrix &beta_est, FRMatrix &se_beta,
                           FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
                           arma::fcolvec &beta_abs_errs, int max_iter,
                           bool is_t_dist)
{

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  arma::span col_1 = arma::span(0, 0);
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist : norm_dist;
#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++)
  {
    checkInterrupt();
    // Initialize beta
    arma::fmat beta(n_parms, 1, arma::fill::zeros);
    arma::fmat cov_w_mat = cov.data;
    arma::fcolvec poi_mat = poi_data.data.col(poi_col);
    arma::ucolvec w2_col = W2.col(poi_col);
    arma::fmat int_w_mat =
        interactions.data % arma::repmat(poi_mat, 1, interactions.data.n_cols);
    arma::uword first_chunk = cov_w_mat.n_cols - 1;
    arma::uword second_chunk = n_parms - 1;
    arma::span first = arma::span(0, first_chunk);
    arma::span second = arma::span(first_chunk + 1, second_chunk);
    // cov_w_mat.each_col([&w2_col](arma::vec &a){a%w2_col;});
    // int_w_mat.each_col([&w2_col](arma::vec &a) {a%w2_col;});
    arma::fmat temp1 = cov_w_mat % arma::repmat(w2_col, 1, cov_w_mat.n_cols);
    arma::fmat temp2 = int_w_mat % arma::repmat(w2_col, 1, int_w_mat.n_cols);

    arma::fmat A(n_parms, n_parms, arma::fill::zeros);
    A.submat(first, first) = temp1.t() * cov_w_mat;        // C'C
    A.submat(second, second) = temp2.t() * int_w_mat;      // I'I
    A.submat(first, second) = temp1.t() * int_w_mat;       // C'I
    A.submat(second, first) = A.submat(first, second).t(); // I'C

    arma::fmat z = w2_col % (pheno.data);

    arma::fmat B = arma::fmat(n_parms, 1, arma::fill::zeros);
    B.submat(first, col_1) = cov_w_mat.t() * z;
    B.submat(second, col_1) = int_w_mat.t() * z;
    // beta = arma::solve(A, B, arma::solve_opts::fast);
    arma::fmat AAinv = arma::inv_sympd(A, inv_opts::allow_approx);
    beta = AAinv * B;
    beta_est.data.col(poi_col) = beta;
    arma::fmat eta = cov_w_mat * beta.submat(first, col_1);
    // arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
    // eta += cov_w_mat * beta.submat(first, col_1);
    // eta += int_w_mat * beta.submat(second, col_1);
    // int df = arma::as_scalar(arma::sum(W2.col(poi_col), 0)) - n_parms;
    // double mse = arma::conv_to<double>::from(W2.col(poi_col).t() *
    //                                          arma::square(pheno.data - eta))
    //                                          /
    //              (double)df;

    // se_beta.data.col(poi_col) =
    //     arma::sqrt(mse * arma::abs(arma::diagvec(arma::pinv(A))));
    // arma::fmat temp_se = se_beta.data.col(poi_col);
    // arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    float df = arma::accu(W2.col(poi_col)) - n_parms;
    arma::fmat tmse = (W2.col(poi_col).t() * arma::square(pheno.data - eta));
    float mse = tmse(0, 0) / df;

    se_beta.data.col(poi_col) = arma::sqrt(mse * arma::diagvec(AAinv));
    arma::fmat temp_se = se_beta.data.col(poi_col);
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}
