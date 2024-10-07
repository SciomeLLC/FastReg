// [[Rcpp::depends(RcppArmadillo)]]
#include <chrono>
#include <iostream>
#include <regression.h>
#include <utils.h>

#include <RcppEigen.h>

template <typename T> Rcpp::NumericVector arma2vec(const T &x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

arma::fcolvec t_dist_r(arma::fcolvec abs_z, int df) {
  arma::fcolvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) = -1 * log10(R::pt(abs_z(i), df, false, false));
  }
  return ret_val;
}

float chisq(float lrs, int df) {
  return R::pchisq(lrs, df, 0, 1) / (-1.0 * log(10));
}

arma::fcolvec norm_dist_r(arma::fcolvec abs_z, int df) {
  arma::fcolvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) = -1 * R::pnorm(abs_z(i), 1.0, 1.0, false, false);
  }
  return ret_val;
}

typedef Eigen::Matrix<arma::uword, Eigen::Dynamic, Eigen::Dynamic> MatrixUd;

arma::fcolvec pmean(arma::fcolvec a, arma::fcolvec b) { return (a + b) / 2; }

//////////////////////////////////////////////////
// @brief creates a matrix that will be modified to create
//        the estimates that are based upon non-normalized data
// @mns   mean of the data
// @stds  standard deviation of the data
// @n     total number of variables to be estimated
//        which is the ultimate size of the matrix
// @start start column to consider
// @return  Matrix used to recenter values
//////////////////////////////////////////////////
arma::fmat createDelta(arma::fcolvec mns, arma::fcolvec stds, unsigned int n,
                       unsigned int start) {
  arma::fmat delta(n, n);
  delta = delta.eye(n, n); // start off with the identity matrix

  for (unsigned int i = start; i < start + mns.size(); i++) {
    if (i == start) // the start location is assumed to be the intercept
    {
      for (int j = start + 1; j < start + mns.size(); j++) {
        delta(i, j) = -mns[j] / stds[j];
      }
    } else {
      delta(i, i) = 1 / stds[i];
    }
  }

  return delta;
}

Eigen::MatrixXf createDelta_EIGEN(const Eigen::VectorXf &mns,
                                  const Eigen::VectorXf &stds, unsigned int n,
                                  unsigned int start) {
  Eigen::MatrixXf delta = Eigen::MatrixXf::Identity(n, n);

  for (unsigned int i = start; i < start + mns.size(); ++i) {
    if (i == start) // the start location is assumed to be the intercept
    {
      for (unsigned int j = start + 1; j < start + mns.size(); ++j) {
        delta(i, j) = -mns(j - start) / stds(j - start);
      }
    } else {
      delta(i, i) = 1.0f / stds(i - start);
    }
  }

  return delta;
}

/////////////////////////////////////////////////
//@brief: Map the identical interaction columns to develop the
//        delta matrix
//@X   : Covariate Matrix
//@X_I : Interaction Matrix
/////////////////////////////////////////////////
arma::umat map_is_interaction(arma::fmat X, arma::fmat X_I) {
  arma::umat interCol(2, X.n_cols + X_I.n_cols);

  for (int i = 1; i < X_I.n_cols; i++) // first col is intercept
  {
    for (int j = 1; j < X.n_cols; j++) // first col is intercept
    {
      arma::ucolvec tmp = (X.col(j) == X_I.col(i));
      if (arma::accu(tmp) == X.n_rows) {
        interCol(0, j) = X.n_cols + i;
        interCol(1, X.n_cols + i) = j;
      }
    }
  }

  return interCol;
}

Eigen::MatrixXi map_is_interaction_EIGEN(const Eigen::MatrixXf &X,
                                         const Eigen::MatrixXf &X_I) {
  Eigen::MatrixXi interCol = Eigen::MatrixXi::Zero(2, X.cols() + X_I.cols());

  for (int i = 1; i < X_I.cols(); ++i) // first col is intercept
  {
    for (int j = 1; j < X.cols(); ++j) // first col is intercept
    {
      // Compare columns element-wise
      if ((X.col(j).array() == X_I.col(i).array()).all()) {
        interCol(0, j) = X.cols() + i;
        interCol(1, X.cols() + i) = j;
      }
    }
  }

  return interCol;
}

// @breif : builds the delta matrix for just the covariates.  This
//          omits the POI delta
// @X     : Covariate X matrix
// @X_I   : Interaction X matrix
// @x_mean: means of the origional covariate matrix 1xn where n is the number
//          of columns in X alone
// @x_sd  : sd of the origional covariate matrix 1xn where n is the number
//          of columns in X alone
arma::fmat build_covariate_delta(arma::fmat X, arma::fmat X_I,
                                 arma::fcolvec x_mean, arma::fcolvec x_sd) {
  ///////////////////////////////////////////////////////////////////////////
  unsigned int n_var = X.n_cols + X_I.n_cols;
  arma::fmat tmeans(n_var, 1, arma::fill::zeros);
  arma::fmat tsds(n_var, 1, arma::fill::ones);
  tmeans.submat(0, 0, x_mean.n_rows - 1, 0) = x_mean;
  tsds.submat(0, 0, x_sd.n_rows - 1, 0) = x_sd;
  arma::fmat delta = createDelta(tmeans, tsds, n_var, 0);

  return delta;
}

Eigen::MatrixXf build_covariate_delta_EIGEN(const Eigen::MatrixXf &X,
                                            const Eigen::MatrixXf &X_I,
                                            const Eigen::VectorXf &x_mean,
                                            const Eigen::VectorXf &x_sd) {
  int n_var = X.cols() + X_I.cols();
  Eigen::VectorXf tmeans = Eigen::VectorXf::Zero(n_var);
  Eigen::VectorXf tsds = Eigen::VectorXf::Ones(n_var);

  tmeans.head(x_mean.size()) = x_mean;
  tsds.head(x_sd.size()) = x_sd;

  Eigen::MatrixXf delta = createDelta_EIGEN(tmeans, tsds, n_var, 0);

  return delta;
}

void LogisticRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                             FRMatrix &interactions, arma::umat &W2,
                             FRMatrix &beta_est, FRMatrix &se_beta,
                             FRMatrix &neglog10_pvl,
                             arma::fcolvec &beta_rel_errs,
                             arma::fcolvec &beta_abs_errs, arma::fcolvec &iters,
                             int max_iter, arma::fmat &x_mean, arma::fmat &x_sd,
                             arma::fmat &xi_mean, arma::fmat &xi_sd,
                             bool is_t_dist) {

  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  arma::fmat X_C = cov.data;
  arma::fmat Y = pheno.data;
  arma::fmat X_I = interactions.data;
  arma::fmat POI = poi_data.data;
  // set up delta-method change
  arma::fmat delta = build_covariate_delta(X_C, X_I, x_mean, x_sd);
  arma::umat int_loc = map_is_interaction(X_C, X_I);
  arma::uword n_parms = X_C.n_cols + X_I.n_cols;
  arma::fmat I(n_parms, n_parms);
  I = I.eye();
  arma::fmat delta2 = I;
  for (arma::uword poi_col = 0; poi_col < POI.n_cols; poi_col++) {
    checkInterrupt();
    arma::fmat A(n_parms, n_parms, arma::fill::zeros);
    // Initialize beta
    arma::fcolvec beta(n_parms, arma::fill::zeros);
    arma::fcolvec beta_old(n_parms, arma::fill::ones);
    arma::fmat cov_w_mat = X_C;
    arma::fmat int_w_mat = X_I;
    arma::ucolvec w2_col = W2.col(poi_col);
    arma::fcolvec POIs = POI.col(poi_col);

    float meanPOI = arma::mean(POIs);
    float sdPOI = arma::stddev(POIs);
    POIs = POIs - meanPOI;
    POIs /= sdPOI;
    // POIs = (POIs - meanPOI)/sdPOI;

    int_w_mat.each_col() %= POIs;
    // Update delta with
    delta(0, X_C.n_cols) = -meanPOI / sdPOI;
    delta(X_C.n_cols, X_C.n_cols) = 1 / sdPOI;
    arma::fmat pdelta = delta; // reset pdelta each iteration

    for (unsigned int i_idx = 1; i_idx < X_I.n_cols; i_idx++) {
      for (unsigned int j_idx = 1; j_idx < X_C.n_cols; j_idx++) {
        if (int_loc(0, j_idx) != 0) {
          // figure out indexes for delta
          unsigned int idx = X_C.n_cols + i_idx;
          unsigned int idx_X = int_loc(1, X_C.n_cols + i_idx);
          //
          float sd_both = 1 / (xi_sd(i_idx, 0) * sdPOI);
          float mean_X = x_mean(j_idx, 0);
          float mean_XI = xi_mean(i_idx, 0);
          // Covariance Estimate
          pdelta(idx_X, idx) = -meanPOI * sd_both;
          // Interaction "Intercept"
          pdelta(X_C.n_cols, idx) = -mean_XI * sd_both;
          // Main "Intercept"
          pdelta(0, idx) = meanPOI * mean_XI * sd_both;
          // covariate effect
          pdelta(idx, idx) = sd_both;
        }
      }
    }
    arma::fmat X = arma::join_rows(cov_w_mat, int_w_mat);
    // arma::span col_1 = arma::span(0,0);
    float rel_errs = 1.0;
    float abs_errs = 1.0;
    arma::fcolvec beta_diff = beta;
    int iter;
    for (iter = 0; iter < (max_iter) && rel_errs > 1e-5; iter++) {
      beta_old = beta;
      arma::fmat eta = X * beta;
      arma::fmat p = 1 / (1 + arma::exp(-eta));
      arma::fmat W1 = p % (1 - p) % w2_col;
      arma::fmat temp1 = X.each_col() % W1;

      A = temp1.t() * X;
      arma::fmat z = w2_col % (Y - p);
      arma::fmat B = X.t() * z;

      beta = beta + arma::solve(A, B, arma::solve_opts::fast);
      beta_diff = arma::abs(beta - beta_old);
      abs_errs = beta_diff.max();
      // iters.at(poi_col) = iter + 1;
      rel_errs = (beta_diff / arma::abs(beta_old)).max();
    }

    abs_errs = beta_diff.max();
    rel_errs = (beta_diff / arma::abs(beta)).max();
    int df = arma::as_scalar(arma::sum(w2_col, 0)) - n_parms;

    // variance covariance matrix
    arma::fmat rV = solve(A, I, arma::solve_opts::fast);
    rV = (pdelta * rV) * pdelta.t();
    // temp std error
    auto temp_buff = arma::conv_to<arma::fcolvec>::from(arma::diagvec(rV));
    arma::fcolvec temp_se = arma::sqrt(temp_buff);

    beta_est.data.col(poi_col) = pdelta * beta; // beta coeff
    se_beta.data.col(poi_col) = temp_se;        // temp std errors
    beta_abs_errs.at(poi_col) = abs_errs;       // convergence abs err
    beta_rel_errs.at(poi_col) = rel_errs;       // convergence rel err
    iters.at(poi_col) = iter;                   // num iterations
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}

void LinearRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                           FRMatrix &interactions, arma::umat &W2,
                           FRMatrix &beta_est, FRMatrix &se_beta,
                           FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
                           arma::fcolvec &beta_abs_errs, arma::fcolvec &iters,
                           int max_iter, arma::fmat &x_mean, arma::fmat &x_sd,
                           arma::fmat &xi_mean, arma::fmat &xi_sd,
                           bool is_t_dist) {

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  arma::span col_1 = arma::span(0, 0);
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  // #pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
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
    iters.at(poi_col) = 1;
    arma::fmat eta = cov_w_mat * beta.submat(first, col_1);
    float df = arma::accu(W2.col(poi_col)) - n_parms;
    arma::fmat tmse = (W2.col(poi_col).t() * arma::square(pheno.data - eta));
    float mse = tmse(0, 0) / df;

    se_beta.data.col(poi_col) = arma::sqrt(mse * arma::diagvec(AAinv));
    arma::fmat temp_se = se_beta.data.col(poi_col);
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}
