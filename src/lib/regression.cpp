// [[Rcpp::depends(RcppArmadillo)]]
#include <regression.h>

arma::fcolvec RegressionBase::t_dist(arma::fcolvec abs_z, int df) {

  arma::fcolvec pvalues = conv_to<fcolvec>::from(
      -1 * (stats2::pt(abs_z, df, true) + log(2)) / log(10));
  return pvalues;
}

arma::fcolvec RegressionBase::norm_dist(arma::fcolvec abs_z, int df) {
  return conv_to<fcolvec>::from(
      -1 * (stats2::pnorm(abs_z, 1.0, 1.0, true) + log(2)) / log(10));
}

arma::fcolvec pmean(arma::fcolvec a, arma::fcolvec b) { return (a + b) / 2; }

void LogisticRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                             FRMatrix &interactions, arma::umat &W2,
                             FRMatrix &beta_est, FRMatrix &se_beta,
                             FRMatrix &neglog10_pvl,
                             arma::fcolvec &beta_rel_errs,
                             arma::fcolvec &beta_abs_errs, int max_iter,
                             bool is_t_dist) {
  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  // create a pointer to the specified distribution function
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist : norm_dist;
  // arma::fcolvec (*dist_func)(arma::fcolvec, int) = t_dist;

#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
    arma::fmat A(n_parms, n_parms, arma::fill::zeros);
    // arma::fmat poi_mat =
    // arma::conv_to<arma::fmat>::from(poi_data.data.col(poi_col));

    // Initialize beta
    arma::fcolvec beta(n_parms, arma::fill::zeros);
    arma::fcolvec beta_old = arma::fcolvec(beta);
    arma::fmat cov_w_mat = cov.data;
    arma::fmat int_w_mat = interactions.data;
    arma::ucolvec w2_col = W2.col(poi_col);
    int_w_mat = interactions.data % arma::repmat(poi_data.data.col(poi_col), 1,
                                                 interactions.data.n_cols);

    arma::uword first_chunk = cov_w_mat.n_cols - 1;
    arma::uword second_chunk = n_parms - 1;
    arma::span first = arma::span(0, first_chunk);
    arma::span second = arma::span(first_chunk + 1, second_chunk);
    // arma::span col_1 = arma::span(0,0);
    for (int iter = 0; iter < max_iter; iter++) {
      arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
      eta += cov_w_mat * beta.subvec(first);
      eta += int_w_mat * beta.subvec(second);

      arma::fmat p = 1 / (1 + arma::exp(-eta));
      arma::fmat W1 = p % (1 - p) % w2_col;

      // cov_w_mat.each_col([&W1](arma::vec &a){a%W1;});
      // int_w_mat.each_col([&W1](arma::vec &a){a%W1;});
      // cov_w_mat.each_col() %= W1;
      // int_w_mat.each_col() %= W1;
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
      if (iter == max_iter - 2) {
        beta_old = beta;
      }
      beta = beta + arma::solve(A, B, arma::solve_opts::fast);
    }

    arma::fcolvec beta_diff = arma::abs(beta - beta_old);
    beta_abs_errs.at(poi_col) = beta_diff.max();
    beta_rel_errs.at(poi_col) = (beta_diff / arma::abs(beta)).max();

    int df = arma::as_scalar(arma::sum(w2_col, 0)) - n_parms;
    arma::fcolvec temp_se = arma::sqrt(arma::abs(arma::diagvec(arma::pinv(A))));

    beta_est.data.col(poi_col) = beta;
    se_beta.data.col(poi_col) =
        arma::sqrt(arma::abs(arma::diagvec(arma::pinv(A))));
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}

void LinearRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                           FRMatrix &interactions, arma::umat &W2,
                           FRMatrix &beta_est, FRMatrix &se_beta,
                           FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
                           arma::fcolvec &beta_abs_errs, int max_iter,
                           bool is_t_dist) {

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  arma::span col_1 = arma::span(0, 0);
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist : norm_dist;
#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
    // Initialize beta
    arma::fmat beta(n_parms, 1, arma::fill::zeros);
    arma::fmat cov_w_mat = cov.data;
    arma::fcolvec poi_mat = poi_data.data.col(poi_col);
    arma::ucolvec w2_col = W2.col(poi_col);
    // arma::fmat int_w_mat = interactions.data;
    // int_w_mat.each_col([&poi_mat, &w2_col](arma::vec&
    // a){(a%poi_mat)%w2_col;});
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
    beta = arma::solve(A, B, arma::solve_opts::fast);

    beta_est.data.col(poi_col) = beta;
    arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
    eta += cov_w_mat * beta.submat(first, col_1);
    eta += int_w_mat * beta.submat(second, col_1);
    int df = arma::as_scalar(arma::sum(W2.col(poi_col), 0)) - n_parms;
    double mse = arma::conv_to<double>::from(W2.col(poi_col).t() *
                                             arma::square(pheno.data - eta)) /
                 (double)df;

    se_beta.data.col(poi_col) =
        arma::sqrt(mse * arma::abs(arma::diagvec(arma::pinv(A))));
    arma::fmat temp_se = se_beta.data.col(poi_col);
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    neglog10_pvl.data.col(poi_col) = (*dist_func)(neg_abs_z, df);
  }
}
