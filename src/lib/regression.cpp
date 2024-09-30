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
    ret_val(i) = -1 * (R::pt(abs_z(i), df, true, true) + log(2)) / log(10);
  }
  // Rcpp::NumericVector dt = Rcpp::dt(arma2vec(abs_z), df, true);
  return ret_val;
}

float chisq(float lrs, int df) {
  return R::pchisq(lrs, df, 0, 1) / (-1.0 * log(10));
}

arma::fcolvec norm_dist_r(arma::fcolvec abs_z, int df) {
  arma::fcolvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) =
        -1 * (R::pnorm(abs_z(i), 1.0, 1.0, true, true) + log(2)) / log(10);
  }
  // Rcpp::NumericVector dt = Rcpp::dt(arma2vec(abs_z), df, true);
  return ret_val;
  // return Rcpp::dnorm(arma2vec(abs_z), 1.0, 1.0, true);
}

typedef Eigen::Matrix<arma::uword, Eigen::Dynamic, Eigen::Dynamic> MatrixUd;

arma::fcolvec pmean(arma::fcolvec a, arma::fcolvec b) { return (a + b) / 2; }

void run_vla_2(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
               arma::fmat &poi_sqrd, FRMatrix &interactions, arma::umat &W2,
               FRMatrix &beta_est, FRMatrix &se_beta, FRMatrix &neglog10_pvl,
               arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
               FRMatrix &beta_est2, FRMatrix &se_beta2, FRMatrix &neglog10_pvl2,
               arma::fcolvec &beta_rel_errs2, arma::fcolvec &beta_abs_errs2,
               arma::fmat &iters, arma::fmat &lls, int max_iter, bool is_t_dist,
               std::vector<int> &poi_2_idx, MatrixUd &W2_EIG,
               Eigen::MatrixXf &W2f, Eigen::MatrixXf &tphenoD,
               FRMatrix &no_interactions) {

  arma::uword n_parms = cov.data.n_cols + no_interactions.data.n_cols;
  // create a pointer to the specified distribution function
  arma::fcolvec (*dist_func_r)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  float ll1, ll2, lrs, lrs_pval;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::fcolvec)poi_data.data.col(0)).size();
  arma::fcolvec temp_se(poi_col_size, arma::fill::zeros);
  // arma::fcolvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::fcolvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::fcolvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::fcolvec poi_col(poi_col_size, arma::fill::zeros);
  Eigen::MatrixXf cov_w_mat = Eigen::Map<Eigen::MatrixXf>(
      cov.data.memptr(), cov.data.n_rows, cov.data.n_cols);
  Eigen::MatrixXf int_w_mat = Eigen::Map<Eigen::MatrixXf>(
      no_interactions.data.memptr(), no_interactions.data.n_rows,
      no_interactions.data.n_cols);
  Eigen::MatrixXf int_w_mat2 = Eigen::Map<Eigen::MatrixXf>(
      interactions.data.memptr(), interactions.data.n_rows,
      interactions.data.n_cols);
  arma::uword n_parms2 = cov.data.n_cols + interactions.data.n_cols;
#pragma omp parallel for
  for (int poi_col_idx : poi_2_idx) {
    checkInterrupt();
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(n_parms, n_parms);
    A(0, 0) = 1;
    A(n_parms, n_parms) = 10;
    Eigen::VectorXf beta = Eigen::VectorXf::Zero(n_parms, 1);
    Eigen::VectorXf beta_old = beta;
    cov_w_mat = Eigen::Map<Eigen::MatrixXf>(cov.data.memptr(), cov.data.n_rows,
                                            cov.data.n_cols);
    int_w_mat = Eigen::Map<Eigen::MatrixXf>(no_interactions.data.memptr(),
                                            no_interactions.data.n_rows,
                                            no_interactions.data.n_cols);
    Eigen::VectorXf w2_col = W2f.col(poi_col_idx);
    arma::fmat POI_ARMA = poi_data.data.col(poi_col_idx);
    Eigen::MatrixXf POI =
        Eigen::Map<Eigen::MatrixXf>(POI_ARMA.memptr(), POI_ARMA.n_rows, 1);
    Eigen::MatrixXf temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXf X =
        Eigen::MatrixXf::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXf beta_diff = (beta - beta_old).array().abs();
    beta_rel_errs.at(poi_col_idx) = 1e8;

    // Fit 2
    Eigen::MatrixXf A2 = Eigen::MatrixXf::Zero(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    Eigen::VectorXf beta2 = Eigen::VectorXf::Zero(n_parms2, 1);
    Eigen::VectorXf beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXf>(interactions.data.memptr(),
                                             interactions.data.n_rows,
                                             interactions.data.n_cols);
    Eigen::VectorXf w2_col2 = W2f.col(poi_col_idx);
    Eigen::MatrixXf POI2 =
        Eigen::Map<Eigen::MatrixXf>(POI_ARMA.memptr(), POI_ARMA.n_rows, 1);
    Eigen::MatrixXf temp_mat_e2 = int_w_mat2;
    int_w_mat2 = temp_mat_e2.array().colwise() * POI2.col(0).array();
    Eigen::MatrixXf X2 = Eigen::MatrixXf::Zero(
        int_w_mat2.rows(), cov_w_mat.cols() + int_w_mat2.cols());
    X2.topLeftCorner(int_w_mat2.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_mat2.rows(), int_w_mat2.cols()) = int_w_mat2;

    Eigen::VectorXf beta_diff2 = (beta2 - beta_old2).array().abs();
    beta_rel_errs2.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    Eigen::MatrixXf eta, eta2, p, p2;

    for (int iter = 0; iter < max_iter && beta_rel_errs.at(poi_col_idx) > 1e-4;
         iter++) {
      // Fit 1
      eta = Eigen::MatrixXf::Zero(cov.data.n_rows, 1);
      eta = X * beta;
      p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();
      Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;
      Eigen::MatrixXf z =
          w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXf B = Eigen::MatrixXf::Zero(n_parms, 1);
      B = X.transpose() * z;
      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      iters.at(poi_col_idx, 0) = iter + 1;
      Eigen::MatrixXf mTemp = (beta_diff.array() / beta_old.array());
      beta_rel_errs.at(poi_col_idx) = mTemp.maxCoeff();
      beta_old = beta;
    }

    // Fit 2
    for (int iter = 0; iter < max_iter && beta_rel_errs2.at(poi_col_idx) > 1e-4;
         iter++) {
      eta2 = Eigen::MatrixXf::Zero(cov.data.n_rows, 1);
      eta2 = X2 * beta2;
      p2 = -eta2.array();
      p2 = p2.array().exp() + 1;
      p2 = 1 / p2.array();
      Eigen::MatrixXf W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      Eigen::MatrixXf temp1_2 = X2.array().colwise() * W1_2.col(0).array();
      A2 = temp1_2.transpose() * X2;
      Eigen::MatrixXf z2 =
          w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXf B2 = Eigen::MatrixXf::Zero(n_parms2, 1);
      B2 = X2.transpose() * z2;
      beta2 = beta2 +
              A2.ldlt().solve(B2); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff2 = (beta2.array() - beta_old2.array()).abs();
      beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      Eigen::MatrixXf mTemp2 = (beta_diff2.array() / beta_old2.array());
      iters.at(poi_col_idx, 1) = iter + 1;
      beta_rel_errs2.at(poi_col_idx) = mTemp2.maxCoeff();
      beta_old2 = beta2;
    }

    beta_abs_errs.at(poi_col_idx) = beta_diff.maxCoeff();
    beta_rel_errs.at(poi_col_idx) =
        (beta_diff.array() / beta.array().abs()).maxCoeff();
    int df = w2_col.array().sum() - n_parms;
    Eigen::MatrixXf temp_inv = A.inverse();
    Eigen::VectorXf diag = temp_inv.diagonal().array().sqrt();

    beta_abs_errs2.at(poi_col_idx) = beta_diff2.maxCoeff();
    beta_rel_errs2.at(poi_col_idx) =
        (beta_diff2.array() / beta2.array().abs()).maxCoeff();
    int df2 = w2_col.array().sum() - n_parms2;
    Eigen::MatrixXf temp_inv2 = A2.inverse();
    Eigen::VectorXf diag2 = temp_inv2.diagonal().array().sqrt();

    // calculate LLs
    arma::fcolvec p_arma = arma::fcolvec(p.data(), p.size(), true, false);
    arma::fcolvec p2_arma = arma::fcolvec(p2.data(), p2.size(), true, false);
    arma::fcolvec eta_arma = arma::fcolvec(eta.data(), eta.size(), true, false);
    arma::fcolvec eta2_arma =
        arma::fcolvec(eta2.data(), eta2.size(), true, false);
    ll1 = arma::accu(log(1.0 - p_arma)) +
          arma::accu(W2.col(poi_col_idx) % eta_arma);
    ll2 = arma::accu(log(1.0 - p2_arma)) +
          arma::accu(W2.col(poi_col_idx) % eta2_arma);
    lrs = 2.0 * (ll1 - ll2);
    lrs_pval = chisq(std::abs(lrs), n_parms2);
    lls.at(poi_col_idx, 0) = ll1;
    lls.at(poi_col_idx, 1) = ll2;
    lls.at(poi_col_idx, 2) = lrs;
    lls.at(poi_col_idx, 3) = lrs_pval;
    lls.at(poi_col_idx, 4) = 2;

    // convert it all back to armadillo matrix types
    temp_se = arma::fcolvec(diag.data(), diag.size(), true, false);
    arma::fcolvec temp_b = arma::fcolvec(beta.data(), beta.size(), true, false);

    beta_est.data.col(poi_col_idx) = temp_b;
    se_beta.data.col(poi_col_idx) = temp_se;

    // convert it all back to armadillo matrix types
    arma::fcolvec temp_se2 =
        arma::fcolvec(diag2.data(), diag2.size(), true, false);
    arma::fcolvec temp_b2 =
        arma::fcolvec(beta2.data(), beta2.size(), true, false);
    beta_est2.data.col(poi_col_idx) = temp_b2;
    se_beta2.data.col(poi_col_idx) = temp_se2;

    // pval for fit1
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    neglog10_pvl.data.col(poi_col_idx) = (*dist_func_r)(neg_abs_z, df);
    // pval for fit2
    arma::fcolvec neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;
    neglog10_pvl2.data.col(poi_col_idx) = (*dist_func_r)(neg_abs_z2, df2);
  } // #pragma omp parallel for
};

// void run_vla_3(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
//                arma::fmat &poi_sqrd, FRMatrix &interactions, arma::umat &W2,
//                FRMatrix &beta_est, FRMatrix &se_beta, FRMatrix &neglog10_pvl,
//                arma::fcolvec &beta_rel_errs, arma::fcolvec &beta_abs_errs,
//                FRMatrix &beta_est2, FRMatrix &se_beta2, FRMatrix &neglog10_pvl2,
//                arma::fcolvec &beta_rel_errs2, arma::fcolvec &beta_abs_errs2,
//                arma::fmat &iters, arma::fmat &lls, int max_iter, bool is_t_dist,
//                std::vector<int> &poi_2_idx, MatrixUd &W2_EIG,
//                Eigen::MatrixXf &W2f, Eigen::MatrixXf &tphenoD,
//                FRMatrix &no_interactions) {
//   arma::fmat interactions_sqrd =
//       arma::join_rows(interactions.data, arma::square(interactions.data));
//   arma::fmat pois = arma::join_rows(poi_data.data, poi_sqrd);
//   arma::uword n_parms = cov.data.n_cols + no_interactions.data.n_cols;
//   // create a pointer to the specified distribution function
//   arma::fcolvec (*dist_func_r)(arma::fcolvec, int) =
//       is_t_dist == true ? t_dist_r : norm_dist_r;
//   float ll1, ll2, lrs, lrs_pval;
//   int df, df2, lrs_df;
//   int poi_col_size = ((arma::fcolvec)poi_data.data.col(0)).size();
//   arma::fcolvec temp_se(poi_col_size, arma::fill::zeros);
//   arma::fcolvec neg_abs_z(poi_col_size, arma::fill::zeros);
//   arma::fcolvec neg_abs_z2(poi_col_size * 2, arma::fill::zeros);
//   arma::fcolvec poi_col(poi_col_size, arma::fill::zeros);
//   arma::fcolvec poi_col2(poi_col_size * 2, arma::fill::zeros);
//   Eigen::MatrixXf cov_w_mat = Eigen::Map<Eigen::MatrixXf>(
//       cov.data.memptr(), cov.data.n_rows, cov.data.n_cols);
//   Eigen::MatrixXf int_w_mat = Eigen::Map<Eigen::MatrixXf>(
//       no_interactions.data.memptr(), no_interactions.data.n_rows,
//       no_interactions.data.n_cols);
//   Eigen::MatrixXf int_w_mat2 = Eigen::Map<Eigen::MatrixXf>(
//       interactions_sqrd.data.memptr(), interactions_sqrd.n_rows,
//       interactions_sqrd.n_cols);
//   arma::uword n_parms2 = cov.data.n_cols + interactions_sqrd.n_cols;

//   arma::uword n_parms2 = cov.data.n_cols + interactions.data.n_cols;
// #pragma omp parallel for
//   for (int poi_col_idx : poi_2_idx) {
//     checkInterrupt();
//     Eigen::MatrixXf A = Eigen::MatrixXf::Zero(n_parms, n_parms);
//     A(0, 0) = 1;
//     A(n_parms, n_parms) = 10;
//     Eigen::VectorXf beta = Eigen::VectorXf::Zero(n_parms, 1);
//     Eigen::VectorXf beta_old = beta;
//     cov_w_mat = Eigen::Map<Eigen::MatrixXf>(cov.data.memptr(), cov.data.n_rows,
//                                             cov.data.n_cols);
//     int_w_mat = Eigen::Map<Eigen::MatrixXf>(no_interactions.data.memptr(),
//                                             no_interactions.data.n_rows,
//                                             no_interactions.data.n_cols);
//     Eigen::VectorXf w2_col = W2f.col(poi_col_idx);
//     arma::fmat POI_ARMA = poi_data.data.col(poi_col_idx);
//     Eigen::MatrixXf POI =
//         Eigen::Map<Eigen::MatrixXf>(POI_ARMA.memptr(), POI_ARMA.n_rows, 1);
//     Eigen::MatrixXf temp_mat_e = int_w_mat;
//     int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
//     Eigen::MatrixXf X =
//         Eigen::MatrixXf::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
//     X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
//     X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

//     Eigen::VectorXf beta_diff = (beta - beta_old).array().abs();
//     beta_rel_errs.at(poi_col_idx) = 1e8;

//     // Fit 2
//     Eigen::MatrixXf A2 = Eigen::MatrixXf::Zero(n_parms2, n_parms2);
//     A2(0, 0) = 1;
//     A2(n_parms2, n_parms2) = 10;
//     Eigen::VectorXf beta2 = Eigen::VectorXf::Zero(n_parms2, 1);
//     Eigen::VectorXf beta_old2 = beta2;
//     int_w_mat2 = Eigen::Map<Eigen::MatrixXf>(interactions.data.memptr(),
//                                              interactions.data.n_rows,
//                                              interactions.data.n_cols);
//     Eigen::VectorXf w2_col2 = W2f.col(poi_col_idx);
//     Eigen::MatrixXf POI2 =
//         Eigen::Map<Eigen::MatrixXf>(POI_ARMA.memptr(), POI_ARMA.n_rows, 1);
//     Eigen::MatrixXf temp_mat_e2 = int_w_mat2;
//     int_w_mat2 = temp_mat_e2.array().colwise() * POI2.col(0).array();
//     Eigen::MatrixXf X2 = Eigen::MatrixXf::Zero(
//         int_w_mat2.rows(), cov_w_mat.cols() + int_w_mat2.cols());
//     X2.topLeftCorner(int_w_mat2.rows(), cov_w_mat.cols()) = cov_w_mat;
//     X2.topRightCorner(int_w_mat2.rows(), int_w_mat2.cols()) = int_w_mat2;

//     Eigen::VectorXf beta_diff2 = (beta2 - beta_old2).array().abs();
//     beta_rel_errs2.at(poi_col_idx) = 1e8;

//     // initialize eta and p for fit 1/2
//     Eigen::MatrixXf eta, eta2, p, p2;

//     for (int iter = 0; iter < max_iter && beta_rel_errs.at(poi_col_idx) > 1e-4;
//          iter++) {
//       // Fit 1
//       eta = Eigen::MatrixXf::Zero(cov.data.n_rows, 1);
//       eta = X * beta;
//       p = -eta.array();
//       p = p.array().exp() + 1;
//       p = 1 / p.array();
//       Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
//       Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
//       A = temp1.transpose() * X;
//       Eigen::MatrixXf z =
//           w2_col.array() * (tphenoD.array() - p.array()).array();
//       Eigen::MatrixXf B = Eigen::MatrixXf::Zero(n_parms, 1);
//       B = X.transpose() * z;
//       beta = beta +
//              A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

//       beta_diff = (beta.array() - beta_old.array()).abs();
//       beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
//       iters.at(poi_col_idx, 0) = iter + 1;
//       Eigen::MatrixXf mTemp = (beta_diff.array() / beta_old.array());
//       beta_rel_errs.at(poi_col_idx) = mTemp.maxCoeff();
//       beta_old = beta;
//     }

//     // Fit 2
//     for (int iter = 0; iter < max_iter && beta_rel_errs2.at(poi_col_idx) > 1e-4;
//          iter++) {
//       eta2 = Eigen::MatrixXf::Zero(cov.data.n_rows, 1);
//       eta2 = X2 * beta2;
//       p2 = -eta2.array();
//       p2 = p2.array().exp() + 1;
//       p2 = 1 / p2.array();
//       Eigen::MatrixXf W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
//       Eigen::MatrixXf temp1_2 = X2.array().colwise() * W1_2.col(0).array();
//       A2 = temp1_2.transpose() * X2;
//       Eigen::MatrixXf z2 =
//           w2_col2.array() * (tphenoD.array() - p2.array()).array();
//       Eigen::MatrixXf B2 = Eigen::MatrixXf::Zero(n_parms2, 1);
//       B2 = X2.transpose() * z2;
//       beta2 = beta2 +
//               A2.ldlt().solve(B2); // arma::solve(A, B, arma::solve_opts::fast);

//       beta_diff2 = (beta2.array() - beta_old2.array()).abs();
//       beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
//       Eigen::MatrixXf mTemp2 = (beta_diff2.array() / beta_old2.array());
//       iters.at(poi_col_idx, 1) = iter + 1;
//       beta_rel_errs2.at(poi_col_idx) = mTemp2.maxCoeff();
//       beta_old2 = beta2;
//     }

//     beta_abs_errs.at(poi_col_idx) = beta_diff.maxCoeff();
//     beta_rel_errs.at(poi_col_idx) =
//         (beta_diff.array() / beta.array().abs()).maxCoeff();
//     int df = w2_col.array().sum() - n_parms;
//     Eigen::MatrixXf temp_inv = A.inverse();
//     Eigen::VectorXf diag = temp_inv.diagonal().array().sqrt();

//     beta_abs_errs2.at(poi_col_idx) = beta_diff2.maxCoeff();
//     beta_rel_errs2.at(poi_col_idx) =
//         (beta_diff2.array() / beta2.array().abs()).maxCoeff();
//     int df2 = w2_col.array().sum() - n_parms2;
//     Eigen::MatrixXf temp_inv2 = A2.inverse();
//     Eigen::VectorXf diag2 = temp_inv2.diagonal().array().sqrt();

//     // calculate LLs
//     arma::fcolvec p_arma = arma::fcolvec(p.data(), p.size(), true, false);
//     arma::fcolvec p2_arma = arma::fcolvec(p2.data(), p2.size(), true, false);
//     arma::fcolvec eta_arma = arma::fcolvec(eta.data(), eta.size(), true, false);
//     arma::fcolvec eta2_arma =
//         arma::fcolvec(eta2.data(), eta2.size(), true, false);
//     ll1 = arma::accu(log(1.0 - p_arma)) +
//           arma::accu(W2.col(poi_col_idx) % eta_arma);
//     ll2 = arma::accu(log(1.0 - p2_arma)) +
//           arma::accu(W2.col(poi_col_idx) % eta2_arma);
//     lrs = 2.0 * (ll1 - ll2);
//     lrs_pval = chisq(std::abs(lrs), n_parms2);
//     lls.at(poi_col_idx, 0) = ll1;
//     lls.at(poi_col_idx, 1) = ll2;
//     lls.at(poi_col_idx, 2) = lrs;
//     lls.at(poi_col_idx, 3) = lrs_pval;
//     lls.at(poi_col_idx, 4) = 2;

//     // convert it all back to armadillo matrix types
//     temp_se = arma::fcolvec(diag.data(), diag.size(), true, false);
//     arma::fcolvec temp_b = arma::fcolvec(beta.data(), beta.size(), true, false);

//     beta_est.data.col(poi_col_idx) = temp_b;
//     se_beta.data.col(poi_col_idx) = temp_se;

//     // convert it all back to armadillo matrix types
//     arma::fcolvec temp_se2 =
//         arma::fcolvec(diag2.data(), diag2.size(), true, false);
//     arma::fcolvec temp_b2 =
//         arma::fcolvec(beta2.data(), beta2.size(), true, false);
//     beta_est2.data.col(poi_col_idx) = temp_b2;
//     se_beta2.data.col(poi_col_idx) = temp_se2;

//     // pval for fit1
//     neg_abs_z = arma::abs(temp_b / temp_se) * -1;
//     neglog10_pvl.data.col(poi_col_idx) = (*dist_func_r)(neg_abs_z, df);
//     // pval for fit2
//     arma::fcolvec neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;
//     neglog10_pvl2.data.col(poi_col_idx) = (*dist_func_r)(neg_abs_z2, df2);
//   } // #pragma omp parallel for
// };

void LogisticRegression::run_vla(
    FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data, arma::fmat &poi_sqrd,
    FRMatrix &interactions, arma::umat &W2, FRMatrix &beta_est,
    FRMatrix &se_beta, FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
    arma::fcolvec &beta_abs_errs, FRMatrix &beta_est2, FRMatrix &se_beta2,
    FRMatrix &neglog10_pvl2, arma::fcolvec &beta_rel_errs2,
    arma::fcolvec &beta_abs_errs2, arma::fmat &iters, arma::fmat &lls,
    int max_iter, bool is_t_dist) {

  FRMatrix no_interactions = interactions;
  no_interactions.data = no_interactions.data.col(0); // first col is poi
  no_interactions.row_names_arr = interactions.row_names_arr;
  no_interactions.row_names = interactions.row_names;
  no_interactions.col_names_arr = {"poi"};
  no_interactions.col_names = {{"poi", 0}};
  arma::fcolvec poi_col;
  std::vector<int> poi_2_idx;
  std::vector<int> poi_3_idx;

#pragma omp parallel for
  for (arma::uword poi_col_idx = 0; poi_col_idx < poi_data.data.n_cols;
       poi_col_idx++) {
    poi_col = poi_data.data.col(poi_col_idx);
    bool has_3_unique = ((arma::fcolvec)arma::unique(poi_col)).size() == 3;
    if (has_3_unique) {
      poi_3_idx.push_back(poi_col_idx);
    } else {
      poi_2_idx.push_back(poi_col_idx);
    }
  }
  // convert to Eigen matrices
  MatrixUd W2_EIG = Eigen::Map<MatrixUd>(W2.memptr(), W2.n_rows, W2.n_cols);
  Eigen::MatrixXf W2f = W2_EIG.cast<float>();
  Eigen::MatrixXf tphenoD = Eigen::Map<Eigen::MatrixXf>(
      pheno.data.memptr(), pheno.data.n_rows, pheno.data.n_cols);

  // run regression on G with 2 values
  run_vla_2(cov, pheno, poi_data, poi_sqrd, interactions, W2, beta_est, se_beta,
            neglog10_pvl, beta_rel_errs, beta_abs_errs, beta_est2, se_beta2,
            neglog10_pvl2, beta_rel_errs2, beta_abs_errs2, iters, lls, max_iter,
            is_t_dist, poi_2_idx, W2_EIG, W2f, tphenoD, no_interactions);

  // run regression on G with 3 values
  // run_vla_3(cov, pheno, poi_data, poi_sqrd, interactions, W2, beta_est, se_beta,
  //           neglog10_pvl, beta_rel_errs, beta_abs_errs, beta_est2, se_beta2,
  //           neglog10_pvl2, beta_rel_errs2, beta_abs_errs2, iters, lls, max_iter,
  //           is_t_dist, poi_2_idx, W2_EIG, W2f, tphenoD, no_interactions);
}

// void run_2() {};
// void run_3() {};

void LogisticRegression::run_BLAS(FRMatrix &cov, FRMatrix &pheno,
                                  FRMatrix &poi_data, FRMatrix &interactions,
                                  arma::umat &W2, FRMatrix &beta_est,
                                  FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                                  arma::fcolvec &beta_rel_errs,
                                  arma::fcolvec &beta_abs_errs,
                                  arma::fcolvec &iters, int max_iter,
                                  bool is_t_dist) {

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  // create a pointer to the specified distribution function
  arma::fcolvec (*dist_func_r)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  // arma::fcolvec (*dist_func)(arma::fcolvec, int) =
  //     is_t_dist == true ? t_dist : norm_dist;
#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
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

    for (int iter = 0; iter < max_iter && beta_rel_errs.at(poi_col) > 1e-4;
         iter++) {
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
      iters.at(poi_col) = iter + 1;
      beta_rel_errs.at(poi_col) = (beta_diff / arma::abs(beta_old)).max();
      beta_old = beta;
    }

    if (beta_abs_errs.at(poi_col) >
        1e-4) { // update the covariance matrix if there
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

    // arma::fcolvec temp_se =
    // arma::sqrt(arma::diagvec(arma::inv_sympd(A)));
    arma::fcolvec temp_se = arma::sqrt(arma::diagvec(arma::pinv(A)));
    beta_est.data.col(poi_col) = beta;
    se_beta.data.col(poi_col) = temp_se;
    arma::fcolvec neg_abs_z = arma::abs(beta / temp_se) * -1;
    // if (poi_col == 0) {
    //   Rcpp::Rcout << "neg_abs_z: " << std::endl;
    //   neg_abs_z.print();
    //   arma::fcolvec dist_r = (*dist_func_r)(neg_abs_z, df);
    //   Rcpp::Rcout << "dist_r: " << std::endl;
    //   dist_r.print();
    //   arma::fcolvec dist = (*dist_func)(neg_abs_z, df);
    //   Rcpp::Rcout << "dist: " << std::endl;
    //   dist.print();
    // }
    neglog10_pvl.data.col(poi_col) = (*dist_func_r)(neg_abs_z, df);
  }
}

void LogisticRegression::run_EIGEN(FRMatrix &cov, FRMatrix &pheno,
                                   FRMatrix &poi_data, FRMatrix &interactions,
                                   arma::umat &W2, FRMatrix &beta_est,
                                   FRMatrix &se_beta, FRMatrix &neglog10_pvl,
                                   arma::fcolvec &beta_rel_errs,
                                   arma::fcolvec &beta_abs_errs,
                                   arma::fcolvec &iters, int max_iter,
                                   bool is_t_dist) {

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  // create a pointer to the specified distribution function

  arma::fcolvec (*dist_func_r)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  // arma::fcolvec (*dist_func)(arma::fcolvec, int) =
  //     is_t_dist == true ? t_dist : norm_dist;
  MatrixUd W2_EIG = Eigen::Map<MatrixUd>(W2.memptr(), W2.n_rows, W2.n_cols);
  Eigen::MatrixXf W2f = W2_EIG.cast<float>();
  Eigen::MatrixXf tphenoD = Eigen::Map<Eigen::MatrixXf>(
      pheno.data.memptr(), pheno.data.n_rows, pheno.data.n_cols);

#pragma omp parallel for
  for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
    checkInterrupt();
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(n_parms, n_parms);
    A(0, 0) = 1;
    A(n_parms - 1, n_parms - 1) = 10;
    // Initialize beta

    Eigen::VectorXf beta = Eigen::VectorXf::Zero(n_parms, 1);
    Eigen::VectorXf beta_old = beta;

    Eigen::MatrixXf cov_w_mat = Eigen::Map<Eigen::MatrixXf>(
        cov.data.memptr(), cov.data.n_rows, cov.data.n_cols);
    Eigen::MatrixXf int_w_mat = Eigen::Map<Eigen::MatrixXf>(
        interactions.data.memptr(), interactions.data.n_rows,
        interactions.data.n_cols);

    Eigen::VectorXf w2_col = W2f.col(poi_col);
    arma::fmat POI_ARMA = poi_data.data.col(poi_col);
    Eigen::MatrixXf POI =
        Eigen::Map<Eigen::MatrixXf>(POI_ARMA.memptr(), POI_ARMA.n_rows, 1);
    Eigen::MatrixXf temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXf X =
        Eigen::MatrixXf::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXf beta_diff = (beta - beta_old).array().abs();
    beta_rel_errs.at(poi_col) = 1e8;

    for (int iter = 0; iter < max_iter && beta_rel_errs.at(poi_col) > 1e-4;
         iter++) {
      Eigen::MatrixXf eta = Eigen::MatrixXf::Zero(cov.data.n_rows, 1);
      // arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
      eta = X * beta;
      Eigen::MatrixXf p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();

      Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;

      Eigen::MatrixXf z =
          w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXf B = Eigen::MatrixXf::Zero(n_parms, 1);

      B = X.transpose() * z;

      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      beta_abs_errs.at(poi_col) = beta_diff.array().maxCoeff();
      iters.at(poi_col) = iter + 1;
      Eigen::MatrixXf mTemp = (beta_diff.array() / beta_old.array());
      beta_rel_errs.at(poi_col) = mTemp.maxCoeff();
      beta_old = beta;
    }

    if (beta_abs_errs.at(poi_col) >
        1e-4) { // update the covariance matrix if there
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
    beta_rel_errs.at(poi_col) =
        (beta_diff.array() / beta.array().abs()).maxCoeff();
    int df = w2_col.array().sum() - n_parms;
    Eigen::MatrixXf temp_inv = A.inverse();
    Eigen::VectorXf diag = temp_inv.diagonal().array().sqrt();

    // convert it all back to armadillo matrix types
    arma::fcolvec temp_se =
        arma::fcolvec(diag.data(), diag.size(), true, false);
    arma::fcolvec temp_b = arma::fcolvec(beta.data(), beta.size(), true, false);
    beta_est.data.col(poi_col) = temp_b;
    se_beta.data.col(poi_col) = temp_se;
    arma::fcolvec neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    // if (poi_col == 0) {
    //   Rcpp::Rcout << "neg_abs_z: " << std::endl;
    //   neg_abs_z.print();
    //   arma::fcolvec dist_r = (*dist_func_r)(neg_abs_z, df);
    //   Rcpp::Rcout << "dist_r: " << std::endl;
    //   dist_r.print();
    //   arma::fcolvec dist = (*dist_func)(neg_abs_z, df);
    //   Rcpp::Rcout << "dist: " << std::endl;
    //   dist.print();
    // }
    neglog10_pvl.data.col(poi_col) = (*dist_func_r)(neg_abs_z, df);
  }
}

void LogisticRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                             FRMatrix &interactions, arma::umat &W2,
                             FRMatrix &beta_est, FRMatrix &se_beta,
                             FRMatrix &neglog10_pvl,
                             arma::fcolvec &beta_rel_errs,
                             arma::fcolvec &beta_abs_errs, arma::fcolvec &iters,
                             int max_iter, bool is_t_dist, bool use_blas) {
  if (use_blas) { // use BLAS if faster than Eigen::MatrixXf
                  // calculations
    // run_BLAS
    run_BLAS(cov, pheno, poi_data, interactions, W2, beta_est, se_beta,
             neglog10_pvl, beta_rel_errs, beta_abs_errs, iters, max_iter,
             is_t_dist);
  } else {
    run_EIGEN(cov, pheno, poi_data, interactions, W2, beta_est, se_beta,
              neglog10_pvl, beta_rel_errs, beta_abs_errs, iters, max_iter,
              is_t_dist);
  }
}

void LinearRegression::run(FRMatrix &cov, FRMatrix &pheno, FRMatrix &poi_data,
                           FRMatrix &interactions, arma::umat &W2,
                           FRMatrix &beta_est, FRMatrix &se_beta,
                           FRMatrix &neglog10_pvl, arma::fcolvec &beta_rel_errs,
                           arma::fcolvec &beta_abs_errs, arma::fcolvec &iters,
                           int max_iter, bool is_t_dist, bool use_blas) {

  arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
  arma::span col_1 = arma::span(0, 0);
  arma::fcolvec (*dist_func)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
#pragma omp parallel for
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
    // arma::fmat eta(cov.data.n_rows, 1, arma::fill::zeros);
    // eta += cov_w_mat * beta.submat(first, col_1);
    // eta += int_w_mat * beta.submat(second, col_1);
    // int df = arma::as_scalar(arma::sum(W2.col(poi_col), 0)) -
    // n_parms; double mse =
    // arma::conv_to<double>::from(W2.col(poi_col).t() *
    //                                          arma::square(pheno.data
    //                                          - eta))
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
