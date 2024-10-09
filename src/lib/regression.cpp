// [[Rcpp::depends(RcppArmadillo)]]
#include <chrono>
#include <iostream>
#include <regression.h>

#include <RcppEigen.h>

static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

void checkInterrupt() {
  if (R_ToplevelExec(chkIntFn, NULL) == FALSE) {
    Rcpp::stop("Received user interrupt. Stopping FastReg...");
  }
}

arma::colvec t_dist_r(arma::colvec abs_z, int df) {
  arma::colvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) = -1 * (R::pt(abs_z(i), df, true, true) + log(2)) / log(10);
  }

  return ret_val;
}

double chisq(double lrs, int df) {
  return R::pchisq(lrs, df, 0, 1) / (-1.0 * log(10));
}

arma::colvec norm_dist_r(arma::colvec abs_z, int df) {
  arma::colvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) =
        -1 * (R::pnorm(abs_z(i), 1.0, 1.0, true, true) + log(2)) / log(10);
  }
  return ret_val;
}

void run_vla_2(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
               VLAResult &result, int max_iter, bool is_t_dist,
               std::vector<int> &poi_2_idx, Eigen::MatrixXd &W2f,
               Eigen::MatrixXd &tphenoD) {
  arma::colvec (*dist_func_r)(arma::colvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  float ll1, ll2, lrs, lrs_pval;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::colvec)poi_data.col(0)).size();
  arma::colvec temp_se(poi_col_size, arma::fill::zeros);
  // arma::colvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::mat poi_col(poi_col_size, 1, arma::fill::zeros);
  Eigen::MatrixXd cov_w_mat =
      Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
  Eigen::MatrixXd int_w_mat, int_w_mat2;
  arma::uword n_parms2 = result.cov_int_names.size();
// #pragma omp parallel for
  for (int poi_col_idx : poi_2_idx) {
    checkInterrupt();
    Eigen::MatrixXd A =
        Eigen::MatrixXd::Zero(result.num_parms, result.num_parms);
    A(0, 0) = 1;
    A(result.num_parms, result.num_parms) = 10;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(result.num_parms, 1);
    Eigen::VectorXd beta_old = beta;
    cov_w_mat =
        Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
    int_w_mat = Eigen::Map<Eigen::MatrixXd>(result.no_interactions->memptr(),
                                            result.no_interactions->n_rows,
                                            result.no_interactions->n_cols);
    Eigen::VectorXd w2_col = W2f.col(poi_col_idx);
    poi_col = poi_data.col(poi_col_idx);
    Eigen::MatrixXd POI =
        Eigen::Map<Eigen::MatrixXd>(poi_col.memptr(), poi_col.n_rows, 1);
    Eigen::MatrixXd temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXd X =
        Eigen::MatrixXd::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXd beta_diff = (beta - beta_old).array().abs();
    result.beta_rel_errs.at(poi_col_idx) = 1e8;

    // Fit 2
    Eigen::MatrixXd A2 = Eigen::MatrixXd::Zero(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    Eigen::VectorXd beta2 = Eigen::VectorXd::Zero(n_parms2, 1);
    Eigen::VectorXd beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXd>(result.interactions->memptr(),
                                             result.interactions->n_rows,
                                             result.interactions->n_cols);
    Eigen::VectorXd w2_col2 = W2f.col(poi_col_idx);
    Eigen::MatrixXd temp_mat_e2 = int_w_mat2;
    int_w_mat2 = temp_mat_e2.array().colwise() * POI.col(0).array();

    Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero(
        int_w_mat2.rows(), cov_w_mat.cols() + int_w_mat2.cols());
    X2.topLeftCorner(int_w_mat2.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_mat2.rows(), int_w_mat2.cols()) = int_w_mat2;

    Eigen::VectorXd beta_diff2 = (beta2 - beta_old2).array().abs();
    result.beta_rel_errs2.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    Eigen::MatrixXd eta, eta2, p, p2;

    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs.at(poi_col_idx) > 1e-4;
         iter++) {
      // Fit 1
      eta = Eigen::MatrixXd::Zero(cov.n_rows, 1);
      eta = X * beta;
      p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();
      Eigen::MatrixXd W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXd temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;
      Eigen::MatrixXd z =
          w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXd B = X.transpose() * z;
      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      result.beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      result.iters.at(poi_col_idx, 0) = iter + 1;
      Eigen::MatrixXd mTemp = (beta_diff.array() / beta_old.array());
      result.beta_rel_errs.at(poi_col_idx) = mTemp.maxCoeff();
      beta_old = beta;
    }

    // Fit 2
    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs2.at(poi_col_idx) > 1e-10;
         iter++) {
      eta2 = Eigen::MatrixXd::Zero(cov.n_rows, 1);
      eta2 = X2 * beta2;
      p2 = -eta2.array();
      p2 = p2.array().exp() + 1;
      p2 = 1 / p2.array();
      Eigen::MatrixXd W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      Eigen::MatrixXd temp1_2 = X2.array().colwise() * W1_2.col(0).array();
      A2 = temp1_2.transpose() * X2;
      Eigen::MatrixXd z2 =
          w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXd B2 = X2.transpose() * z2;
      beta2 = beta2 +
              A2.ldlt().solve(B2); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff2 = (beta2.array() - beta_old2.array()).abs();
      result.beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      Eigen::MatrixXd mTemp2 = (beta_diff2.array() / beta_old2.array());
      result.iters.at(poi_col_idx, 1) = iter + 1;
      result.beta_rel_errs2.at(poi_col_idx) = mTemp2.maxCoeff();
      beta_old2 = beta2;
    }

    int df = w2_col.array().sum() - result.num_parms;
    Eigen::VectorXd diag = A.inverse().diagonal().array().sqrt();

    int df2 = w2_col.array().sum() - n_parms2;
    Eigen::VectorXd diag2 = A2.inverse().diagonal().array().sqrt();

    // calculate LLs
    eta = X * beta;
    p = -eta.array();
    p = p.array().exp() + 1;
    p = 1 / p.array();
    eta2 = X2 * beta2;
    p2 = -eta2.array();
    p2 = p2.array().exp() + 1;
    p2 = 1 / p2.array();
    arma::colvec p_arma = arma::colvec(p.data(), p.size(), true, false);
    arma::colvec p2_arma = arma::colvec(p2.data(), p2.size(), true, false);
    arma::colvec eta_arma = arma::colvec(eta.data(), eta.size(), true, false);
    arma::colvec eta2_arma =
        arma::colvec(eta2.data(), eta2.size(), true, false);

    ll1 = arma::accu((log(1.0 - p_arma) + pheno.col(0) % eta_arma) %
                     result.W2.col(poi_col_idx));

    p2_arma.elem(arma::find_nonfinite(p2_arma)).zeros();
    p2_arma.elem(arma::find(p2_arma >= 1.0)).fill(1.0 - 1e-4);
    eta2_arma.elem(arma::find_nonfinite(eta2_arma)).zeros();
    ll2 = arma::accu((log(1.0 - p2_arma) + pheno.col(0) % eta2_arma) %
                     result.W2.col(poi_col_idx));
    lrs = 2.0 * (ll2 - ll1);
    lrs_pval = chisq(lrs, n_parms2 - result.num_parms);

    // calc and set betas
    temp_se = arma::colvec(diag.data(), diag.size(), true, false);
    arma::colvec temp_b = arma::colvec(beta.data(), beta.size(), true, false);
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    arma::colvec temp_se2 =
        arma::colvec(diag2.data(), diag2.size(), true, false);
    arma::colvec temp_b2 =
        arma::colvec(beta2.data(), beta2.size(), true, false);
    neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(X2);
    auto rank = lu_decomp.rank();
    result.set_lls(ll1, ll2, lrs, lrs_pval, 2, poi_col_idx, rank);
    arma::colvec pval = (*dist_func_r)(neg_abs_z, df);
    arma::colvec pval2 = (*dist_func_r)(neg_abs_z2, df2);
    result.set_betas_fit1(temp_b, temp_se, pval, poi_col_idx);
    result.set_betas_fit2(temp_b2, temp_se2, pval2, poi_col_idx);
  } // #pragma omp parallel for
};

void run_vla_3(arma::mat &cov, arma::mat &pheno, arma::mat &poi_data,
               VLAResult &result, int max_iter, bool is_t_dist,
               std::vector<int> &poi_3_idx, Eigen::MatrixXd &W2f,
               Eigen::MatrixXd &tphenoD) {
  arma::colvec (*dist_func_r)(arma::colvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  double ll1, ll2, lrs, lrs_pval;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::colvec)poi_data.col(0)).size();
  arma::colvec temp_se(poi_col_size, arma::fill::zeros);
  // arma::colvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::mat poi_col(poi_col_size, 1, arma::fill::zeros);
  arma::mat poi_col_sqrd(poi_col_size, 1, arma::fill::zeros);
  Eigen::MatrixXd cov_w_mat =
      Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
  Eigen::MatrixXd int_w_mat, int_w_mat2, int_w_mat2_sqrd;
  arma::uword n_parms2 = result.cov_int_names_sqrd.size();
  // #pragma omp parallel for
  for (int poi_col_idx : poi_3_idx) {
    checkInterrupt();
    Eigen::MatrixXd A =
        Eigen::MatrixXd::Zero(result.num_parms, result.num_parms);
    A(0, 0) = 1;
    A(result.num_parms, result.num_parms) = 10;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(result.num_parms, 1);
    Eigen::VectorXd beta_old = beta;
    cov_w_mat =
        Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
    int_w_mat = Eigen::Map<Eigen::MatrixXd>(result.no_interactions->memptr(),
                                            result.no_interactions->n_rows,
                                            result.no_interactions->n_cols);
    Eigen::VectorXd w2_col = W2f.col(poi_col_idx);
    poi_col = poi_data.col(poi_col_idx);

    Eigen::MatrixXd POI =
        Eigen::Map<Eigen::MatrixXd>(poi_col.memptr(), poi_col.n_rows, 1);
    Eigen::MatrixXd temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXd X =
        Eigen::MatrixXd::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXd beta_diff = (beta - beta_old).array().abs();
    result.beta_rel_errs.at(poi_col_idx) = 1e8;

    // Fit 2
    poi_col_sqrd = result.poi_sqrd_mat.col(poi_col_idx);
    Eigen::MatrixXd POI_sqrd = Eigen::Map<Eigen::MatrixXd>(
        poi_col_sqrd.memptr(), poi_col_sqrd.n_rows, 1);
    Eigen::MatrixXd A2 = Eigen::MatrixXd::Zero(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    Eigen::VectorXd beta2 = Eigen::VectorXd::Zero(n_parms2, 1);
    Eigen::VectorXd beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXd>(result.interactions->memptr(),
                                             result.interactions->n_rows,
                                             result.interactions->n_cols);
    int_w_mat2_sqrd = Eigen::Map<Eigen::MatrixXd>(
        result.interactions_sqrd->memptr(), result.interactions_sqrd->n_rows,
        result.interactions_sqrd->n_cols);
    Eigen::VectorXd w2_col2 = W2f.col(poi_col_idx);
    Eigen::MatrixXd temp_mat_e2 = int_w_mat2;
    Eigen::MatrixXd temp_mat_e2_sqrd = int_w_mat2_sqrd;
    int_w_mat2 = temp_mat_e2.array().colwise() * POI.col(0).array();
    int_w_mat2_sqrd =
        temp_mat_e2_sqrd.array().colwise() * POI_sqrd.col(0).array();

    Eigen::MatrixXd int_w_sqrd(int_w_mat2.rows(),
                               int_w_mat2.cols() + int_w_mat2_sqrd.cols());
    int_w_sqrd << int_w_mat2, int_w_mat2_sqrd;
    Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero(
        int_w_mat2.rows(), cov_w_mat.cols() + int_w_sqrd.cols());
    X2.topLeftCorner(int_w_sqrd.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_sqrd.rows(), int_w_sqrd.cols()) = int_w_sqrd;

    Eigen::VectorXd beta_diff2 = (beta2 - beta_old2).array().abs();
    result.beta_rel_errs2.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    Eigen::MatrixXd eta, eta2, p, p2;

    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs.at(poi_col_idx) > 1e-4;
         iter++) {
      // Fit 1
      eta = Eigen::MatrixXd::Zero(cov.n_rows, 1);
      eta = X * beta;
      p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();
      Eigen::MatrixXd W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXd temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;
      Eigen::MatrixXd z =
          w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXd B = X.transpose() * z;
      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      result.beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      result.iters.at(poi_col_idx, 0) = iter + 1;
      Eigen::MatrixXd mTemp = (beta_diff.array() / beta_old.array());
      result.beta_rel_errs.at(poi_col_idx) = mTemp.maxCoeff();
      beta_old = beta;
    }

    // Fit 2
    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs2.at(poi_col_idx) > 1e-10;
         iter++) {
      eta2 = Eigen::MatrixXd::Zero(cov.n_rows, 1);
      eta2 = X2 * beta2;
      p2 = -eta2.array();
      p2 = p2.array().exp() + 1;
      p2 = 1 / p2.array();
      Eigen::MatrixXd W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      Eigen::MatrixXd temp1_2 = X2.array().colwise() * W1_2.col(0).array();
      A2 = temp1_2.transpose() * X2;
      Eigen::MatrixXd z2 =
          w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXd B2 = X2.transpose() * z2;
      beta2 = beta2 +
              A2.ldlt().solve(B2); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff2 = (beta2.array() - beta_old2.array()).abs();
      result.beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      Eigen::MatrixXd mTemp2 = (beta_diff2.array() / beta_old2.array());
      result.iters.at(poi_col_idx, 1) = iter + 1;
      result.beta_rel_errs2.at(poi_col_idx) = mTemp2.maxCoeff();
      beta_old2 = beta2;
    }

    int df = w2_col.array().sum() - result.num_parms;
    Eigen::VectorXd diag = A.inverse().diagonal().array().sqrt();

    int df2 = w2_col.array().sum() - n_parms2;
    Eigen::VectorXd diag2 = A2.inverse().diagonal().array().sqrt();

    // calculate LLs
    eta = X * beta;
    p = -eta.array();
    p = p.array().exp() + 1;
    p = 1 / p.array();
    eta2 = X2 * beta2;

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(X2);
    auto rank = lu_decomp.rank();
    p2 = -eta2.array();
    p2 = p2.array().exp() + 1;
    p2 = 1 / p2.array();
    arma::colvec p_arma = arma::colvec(p.data(), p.size(), true, false);
    arma::colvec p2_arma = arma::colvec(p2.data(), p2.size(), true, false);
    arma::colvec eta_arma = arma::colvec(eta.data(), eta.size(), true, false);
    arma::colvec eta2_arma =
        arma::colvec(eta2.data(), eta2.size(), true, false);

    ll1 = arma::accu((log(1.0 - p_arma) + pheno.col(0) % eta_arma) %
                     result.W2.col(poi_col_idx));
    p2_arma.elem(arma::find_nonfinite(p2_arma)).zeros();
    p2_arma.elem(arma::find(p2_arma >= 1.0)).fill(1.0 - 1e-4);
    eta2_arma.elem(arma::find_nonfinite(eta2_arma)).zeros();
    ll2 = arma::accu((log(1.0 - p2_arma) + pheno.col(0) % eta2_arma) %
                     result.W2.col(poi_col_idx));
    lrs = 2.0 * (ll2 - ll1);
    lrs_pval = chisq(lrs, n_parms2 - result.num_parms);
    // calc and set betas
    temp_se = arma::colvec(diag.data(), diag.size(), true, false);
    arma::colvec temp_b = arma::colvec(beta.data(), beta.size(), true, false);
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    arma::colvec temp_se2 =
        arma::colvec(diag2.data(), diag2.size(), true, false);
    arma::colvec temp_b2 =
        arma::colvec(beta2.data(), beta2.size(), true, false);
    neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;

    result.set_lls(ll1, ll2, lrs, lrs_pval, 3, poi_col_idx, rank);
    arma::colvec pval = (*dist_func_r)(neg_abs_z, df);
    arma::colvec pval2 = (*dist_func_r)(neg_abs_z2, df2);
    result.set_betas_fit1(temp_b, temp_se, pval, poi_col_idx);
    result.set_betas_fit2_sqrd(temp_b2, temp_se2, pval2, poi_col_idx);
  } // #pragma omp parallel for
};

void LogisticRegression::run_vla(arma::mat &cov, arma::mat &pheno,
                                 arma::mat &poi_data, VLAResult &result,
                                 int max_iter, bool is_t_dist) {
  arma::colvec poi_col;
  std::vector<int> poi_2_idx;
  std::vector<int> poi_3_idx;

#pragma omp parallel for
  for (arma::uword poi_col_idx = 0; poi_col_idx < poi_data.n_cols;
       poi_col_idx++) {
    poi_col = poi_data.col(poi_col_idx);
    bool has_3_unique = ((arma::colvec)arma::unique(poi_col)).size() == 3;
    if (has_3_unique) {
      poi_3_idx.push_back(poi_col_idx);
    } else {
      poi_2_idx.push_back(poi_col_idx);
    }
  }

  Eigen::MatrixXd tphenoD =
      Eigen::Map<Eigen::MatrixXd>(pheno.memptr(), pheno.n_rows, pheno.n_cols);
  Eigen::MatrixXd W2f = Eigen::Map<Eigen::MatrixXd>(
      result.W2.memptr(), result.W2.n_rows, result.W2.n_cols);

  // run regression on G with 2 values
  run_vla_2(cov, pheno, poi_data, result, max_iter, is_t_dist, poi_2_idx, W2f,
            tphenoD);
  run_vla_3(cov, pheno, poi_data, result, max_iter, is_t_dist, poi_3_idx, W2f,
            tphenoD);
}
