#include <regression.h>

arma::colvec t_dist_r(arma::colvec abs_z, int df) {
  arma::colvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val[i] = -1 * (R::pt(abs_z[i], df, true, true) + log(2)) / log(10);
  }
  return ret_val;
}

double chisq(double lrs, int df) {
  return R::pchisq(lrs, df, 0, 1) / (-1.0 * log(10));
}

arma::colvec norm_dist_r(arma::colvec abs_z, int df) {
  arma::colvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val[i] =
        -1 * (R::pnorm(abs_z[i], 1.0, 1.0, true, true) + log(2)) / log(10);
  }
  return ret_val;
}

void LogisticRegression::run_vla_2(arma::mat &cov, arma::mat &pheno,
                                   arma::mat &poi_data, VLAResult &result,
                                   int max_iter, bool is_t_dist,
                                   std::vector<int> &poi_2_idx,
                                   Eigen::MatrixXd &W2f,
                                   Eigen::MatrixXd &tphenoD) {
  arma::colvec (*dist_func_r)(arma::colvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  Eigen::MatrixXd X, X2, A, A2, temp1, temp1_2, z, z2, eta, eta2, p, p2,
          temp_eta, temp_eta2, temp_p, temp_p2;
  Eigen::VectorXd beta, beta_old, beta2, beta_old2, beta_diff, beta_diff2,
      w2_col, w2_col2, diag, diag2;
  arma::colvec pval, pval2;
  double ll1, ll1_old, ll2, ll2_old, lrs, lrs_pval, ll1_diff, ll2_diff,
      ll1_rel_err, ll2_rel_err;
  ll1 = ll2 = ll1_old = ll2_old = 1.0;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::colvec)poi_data.col(0)).size();
  arma::colvec temp_se(poi_col_size, arma::fill::zeros);
  arma::colvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::mat poi_col(poi_col_size, 1, arma::fill::zeros);
  arma::mat poi_col_sqrd(poi_col_size, 1, arma::fill::zeros);
  Eigen::MatrixXd cov_w_mat =
      Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
  Eigen::MatrixXd int_w_mat, int_w_mat2, int_w_mat2_sqrd;
  arma::uword n_parms2 = result.cov_int_names.size();
  Eigen::MatrixXd POI;
  X.resize(0, 0);
  X2.resize(0, 0);
  POI.resize(poi_data.n_rows, 1);
  // #pragma omp parallel for
  for (int poi_col_idx : poi_2_idx) {
    checkInterrupt();
    A.resize(result.num_parms, result.num_parms);
    A(0, 0) = 1;
    A(result.num_parms, result.num_parms) = 10;
    beta = Eigen::VectorXd::Zero(result.num_parms, 1);
    beta_old = beta;
    cov_w_mat =
        Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
    int_w_mat = Eigen::Map<Eigen::MatrixXd>(result.no_interactions->memptr(),
                                            result.no_interactions->n_rows,
                                            result.no_interactions->n_cols);
    w2_col = W2f.col(poi_col_idx);
    poi_col = poi_data.col(poi_col_idx);
    POI = Eigen::Map<Eigen::MatrixXd>(poi_col.memptr(), poi_col.n_rows, 1);
    Eigen::MatrixXd temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    X = Eigen::MatrixXd::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    beta_diff = (beta - beta_old).array().abs();
    result.beta_rel_errs.at(poi_col_idx) = 1e8;

    // Fit 2
    A2.resize(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    beta2 = Eigen::VectorXd::Zero(n_parms2, 1);
    beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXd>(result.interactions->memptr(),
                                             result.interactions->n_rows,
                                             result.interactions->n_cols);
    w2_col2 = W2f.col(poi_col_idx);
    Eigen::MatrixXd temp_mat_e2 = int_w_mat2;
    int_w_mat2 = temp_mat_e2.array().colwise() * POI.col(0).array();

    X2 = Eigen::MatrixXd::Zero(int_w_mat2.rows(),
                               cov_w_mat.cols() + int_w_mat2.cols());
    X2.topLeftCorner(int_w_mat2.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_mat2.rows(), int_w_mat2.cols()) = int_w_mat2;

    beta_diff2 = (beta2 - beta_old2).array().abs();
    result.beta_rel_errs2.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    eta.noalias() = X * beta;
    p = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
    ll1_rel_err = 1.0;
    for (int iter = 0;
         iter < max_iter && ll1_rel_err > 1e-7;
         iter++) {
      // Fit 1
      Eigen::MatrixXd W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXd temp1 = X.array().colwise() * W1.col(0).array();
      A.noalias() = temp1.transpose() * X;
      z = w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXd B = X.transpose() * z;
      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      result.beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      result.iters.at(poi_col_idx, 0) = iter + 1;
      result.beta_rel_errs.at(poi_col_idx) =
          (beta_diff.array() / beta_old.array()).maxCoeff();
      beta_old = beta;
      eta.noalias() = X * beta;
      p = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
      ll1_old = ll1;
      temp_p = (p.array().isFinite()).select(p.array(), 0.0);
      temp_p = temp_p.cwiseMin(1.0 - 1e-7);
      temp_eta = (eta.array().isFinite()).select(eta.array(), 0.0);
      ll1 = (((1.0 - temp_p.array()).array().log() +
              tphenoD.col(0).array() * temp_eta.array())
                 .array() *
             w2_col.array())
                .sum();
      ll1_diff = std::abs(ll1 - ll1_old);
      ll1_rel_err = std::abs(ll1_diff / (std::abs(ll1) + 0.05));
    }

    eta2.noalias() = X2 * beta2;
    p2 = (1.0 / (1.0 + (-eta2.array()).exp())).matrix();
    // Fit 2
    for (int iter = 0;
         iter < max_iter && ll2_rel_err > 1e-7;
         iter++) {
      Eigen::MatrixXd W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      temp1_2.noalias() = (X2.array().colwise() * W1_2.col(0).array()).matrix();
      A2.noalias() = temp1_2.transpose() * X2;
      z2 = w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXd B2 = X2.transpose() * z2;
      beta2.noalias() += A2.ldlt().solve(B2);

      beta_diff2 = (beta2 - beta_old2).array().abs();
      result.beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      result.iters.at(poi_col_idx, 1) = iter + 1;
      result.beta_rel_errs2.at(poi_col_idx) =
          (beta_diff2.array() / beta_old2.array()).maxCoeff();
      beta_old2 = beta2;
      eta2.noalias() = X2 * beta2;
      p2 = (1.0 / (1.0 + (-eta2.array()).exp())).matrix();
      ll2_old = ll2;
      temp_p2 = (p2.array().isFinite()).select(p2.array(), 0.0);
      temp_p2 = temp_p2.cwiseMin(1.0 - 1e-7);
      temp_eta2 = (eta2.array().isFinite()).select(eta2.array(), 0.0);
      ll2 = (((1.0 - temp_p2.array()).array().log() +
              tphenoD.col(0).array() * temp_eta2.array())
                 .array() *
             w2_col2.array())
                .sum();
      ll2_diff = std::abs(ll2_old - ll2);
      ll2_rel_err = std::abs(ll2_diff / (std::abs(ll2) + 0.05));
    }
    
    result.lls.at(poi_col_idx, 8) = ll1_diff;
    result.lls.at(poi_col_idx, 9) = ll1_rel_err;
    result.lls.at(poi_col_idx, 10) = ll2_diff;
    result.lls.at(poi_col_idx, 11) = ll2_rel_err;

    int df = w2_col.array().sum() - result.num_parms;
    diag = A.inverse().diagonal().array().sqrt();

    int df2 = w2_col.array().sum() - n_parms2;
    diag2 = A2.inverse().diagonal().array().sqrt();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(X2);
    auto rank = lu_decomp.rank();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp2(X);
    auto rank1 = lu_decomp2.rank();
    result.lls.at(poi_col_idx, 12) = rank1;

    lrs = 2.0 * (ll2 - ll1);
    // lrs_pval = std::abs(chisq(lrs, n_parms2 - result.num_parms));
    lrs_pval = std::abs(chisq(lrs, rank - rank1));

    // calc and set betas
    temp_se = arma::colvec(diag.data(), diag.size(), true, false);
    arma::colvec temp_b = arma::colvec(beta.data(), beta.size(), true, false);
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    temp_se2 = arma::colvec(diag2.data(), diag2.size(), true, false);
    arma::colvec temp_b2 =
        arma::colvec(beta2.data(), beta2.size(), true, false);
    neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;

    result.set_lls(ll1, ll2, lrs, lrs_pval, 2, poi_col_idx, rank);
    arma::colvec pval = (*dist_func_r)(neg_abs_z, df);
    arma::colvec pval2 = (*dist_func_r)(neg_abs_z2, df2);

    result.set_betas_fit1(temp_b, temp_se, pval, poi_col_idx);
    result.set_betas_fit2(temp_b2, temp_se2, pval2, poi_col_idx);
  } // #pragma omp parallel for
};

void LogisticRegression::run_vla_3(arma::mat &cov, arma::mat &pheno,
                                   arma::mat &poi_data, VLAResult &result,
                                   int max_iter, bool is_t_dist,
                                   std::vector<int> &poi_3_idx,
                                   Eigen::MatrixXd &W2f,
                                   Eigen::MatrixXd &tphenoD) {

  arma::colvec (*dist_func_r)(arma::colvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;

  Eigen::MatrixXd X, X2, A, A2, temp1, temp1_2, z, z2, eta, eta2, p, p2,
      temp_eta, temp_eta2, temp_p, temp_p2;
  Eigen::VectorXd beta, beta_old, beta2, beta_old2, beta_diff, beta_diff2,
      w2_col, w2_col2, diag, diag2;
  arma::colvec pval, pval2;
  double ll1, ll1_old, ll2, ll2_old, lrs, lrs_pval, ll1_diff, ll2_diff,
      ll1_rel_err, ll2_rel_err;
  ll1 = ll2 = ll1_old = ll2_old = 1.0;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::colvec)poi_data.col(0)).size();
  arma::colvec temp_se(poi_col_size, arma::fill::zeros);
  arma::colvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::colvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::mat poi_col(poi_col_size, 1, arma::fill::zeros);
  arma::mat poi_col_sqrd(poi_col_size, 1, arma::fill::zeros);
  Eigen::MatrixXd cov_w_mat =
      Eigen::Map<Eigen::MatrixXd>(cov.memptr(), cov.n_rows, cov.n_cols);
  Eigen::MatrixXd int_w_mat, int_w_mat2, int_w_mat2_sqrd;
  arma::uword n_parms2 = result.cov_int_names_sqrd.size();
  Eigen::MatrixXd POI, POI_sqrd;
  X.resize(0, 0);
  X2.resize(0, 0);
  POI.resize(poi_data.n_rows, 1);
  POI_sqrd.resize(result.poi_sqrd_mat.n_rows, 1);
  // #pragma omp parallel for private(                                              \
//     X, X2, A, A2, temp1, temp1_2, z, z2, eta, eta2, p, p2, beta, beta_old,     \
//     beta2, beta_old2, beta_diff, beta_diff2, w2_col, w2_col2, diag, diag2,     \
//     temp_se, neg_abs_z, temp_se2, neg_abs_z2, poi_col, poi_col_sqrd, pval,     \
//     pval2)
  for (int poi_col_idx : poi_3_idx) {
    checkInterrupt();
    A = Eigen::MatrixXd::Zero(result.num_parms, result.num_parms);
    A(0, 0) = 1;
    A(result.num_parms, result.num_parms) = 10;
    beta = Eigen::VectorXd::Zero(result.num_parms, 1);
    beta_old = beta;
    int_w_mat = Eigen::Map<Eigen::MatrixXd>(result.no_interactions->memptr(),
                                            result.no_interactions->n_rows,
                                            result.no_interactions->n_cols);
    w2_col = W2f.col(poi_col_idx);
    poi_col = poi_data.col(poi_col_idx);

    POI = Eigen::Map<Eigen::MatrixXd>(poi_col.memptr(), poi_col.n_rows, 1);
    int_w_mat = int_w_mat.array().colwise() * POI.col(0).array();

    // Construct the design matrix X
    X.resize(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    result.beta_rel_errs.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    eta.noalias() = X * beta;
    p = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
    ll1_rel_err = 1.0;

    for (int iter = 0;
         iter < max_iter && ll1_rel_err > 1e-7;
         iter++) {
      eta.noalias() = X * beta;
      p = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
      Eigen::MatrixXd W1 = p.array() * (1 - p.array()) * w2_col.array();
      temp1.noalias() = (X.array().colwise() * W1.col(0).array()).matrix();
      A.noalias() = temp1.transpose() * X;
      z = w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXd B = X.transpose() * z;
      // arma::solve(A, B, arma::solve_opts::fast);
      beta.noalias() += A.ldlt().solve(B);

      beta_diff = (beta - beta_old).array().abs();
      result.beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      result.iters.at(poi_col_idx, 0) = iter + 1;
      result.beta_rel_errs.at(poi_col_idx) =
          (beta_diff.array() / beta_old.array()).maxCoeff();
      beta_old = beta;
      eta.noalias() = X * beta;
      p = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
      ll1_old = ll1;
      temp_p = (p.array().isFinite()).select(p.array(), 0.0);
      temp_p = temp_p.cwiseMin(1.0 - 1e-7);
      temp_eta = (eta.array().isFinite()).select(eta.array(), 0.0);
      ll1 = (((1.0 - temp_p.array()).array().log() +
              tphenoD.col(0).array() * temp_eta.array())
                 .array() *
             w2_col.array())
                .sum();
      ll1_diff = std::abs(ll1_old - ll1);
      // ll1_rel_err = std::abs(ll1_diff / ll1_old);
      ll1_rel_err = std::abs(ll1_diff / (std::abs(ll1) + 0.05));
    }

    // Fit 2
    poi_col_sqrd = result.poi_sqrd_mat.col(poi_col_idx);
    POI_sqrd = Eigen::Map<Eigen::MatrixXd>(poi_col_sqrd.memptr(),
                                           poi_col_sqrd.n_rows, 1);
    A2 = Eigen::MatrixXd::Zero(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    beta2 = Eigen::VectorXd::Zero(n_parms2, 1);
    beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXd>(result.interactions->memptr(),
                                             result.interactions->n_rows,
                                             result.interactions->n_cols);
    int_w_mat2_sqrd = Eigen::Map<Eigen::MatrixXd>(
        result.interactions_sqrd->memptr(), result.interactions_sqrd->n_rows,
        result.interactions_sqrd->n_cols);
    w2_col2 = W2f.col(poi_col_idx);
    int_w_mat2 = int_w_mat2.array().colwise() * POI.col(0).array();
    int_w_mat2_sqrd =
        int_w_mat2_sqrd.array().colwise() * POI_sqrd.col(0).array();

    Eigen::MatrixXd int_w_sqrd(int_w_mat2.rows(),
                               int_w_mat2.cols() + int_w_mat2_sqrd.cols());
    int_w_sqrd << int_w_mat2, int_w_mat2_sqrd;

    X2.resize(int_w_mat2.rows(), cov_w_mat.cols() + int_w_sqrd.cols());
    X2.topLeftCorner(int_w_sqrd.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_sqrd.rows(), int_w_sqrd.cols()) = int_w_sqrd;

    result.beta_rel_errs2.at(poi_col_idx) = 1e8;

    eta2.noalias() = X2 * beta2;
    p2 = (1.0 / (1.0 + (-eta2.array()).exp())).matrix();
    ll2_rel_err = 1.0;
    // Fit 2
    for (int iter = 0;
         iter < max_iter && ll2_rel_err > 1e-7;
         iter++) {
      Eigen::MatrixXd W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      temp1_2.noalias() = (X2.array().colwise() * W1_2.col(0).array()).matrix();
      A2.noalias() = temp1_2.transpose() * X2;
      z2 = w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXd B2 = X2.transpose() * z2;
      // arma::solve(A, B, arma::solve_opts::fast);
      beta2.noalias() += A2.ldlt().solve(B2);

      beta_diff2 = (beta2 - beta_old2).array().abs();
      result.beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      result.iters.at(poi_col_idx, 1) = iter + 1;
      result.beta_rel_errs2.at(poi_col_idx) =
          (beta_diff2.array() / beta_old2.array()).maxCoeff();
      beta_old2 = beta2;
      eta2.noalias() = X2 * beta2;
      p2 = (1.0 / (1.0 + (-eta2.array()).exp())).matrix();

      // WARNING: Re-using p/eta variables as temporary to avoid memory
      // reallocations
      ll2_old = ll2;
      temp_p2 = (p2.array().isFinite()).select(p2.array(), 0.0);
      temp_p2 = temp_p2.cwiseMin(1.0 - 1e-7);
      temp_eta2 = (eta2.array().isFinite()).select(eta2.array(), 0.0);
      ll2 = (((1.0 - temp_p2.array()).array().log() +
              tphenoD.col(0).array() * temp_eta2.array())
                 .array() *
             w2_col2.array())
                .sum();
      ll2_diff = std::abs(ll2_old - ll2);
      ll2_rel_err = std::abs(ll2_diff / (std::abs(ll2) + 0.05));
    }

    result.lls.at(poi_col_idx, 8) = ll1_diff;
    result.lls.at(poi_col_idx, 9) = ll1_rel_err;
    result.lls.at(poi_col_idx, 10) = ll2_diff;
    result.lls.at(poi_col_idx, 11) = ll2_rel_err;

    int df = w2_col.array().sum() - result.num_parms;
    diag = A.inverse().diagonal().array().sqrt();

    int df2 = w2_col.array().sum() - n_parms2;
    diag2 = A2.inverse().diagonal().array().sqrt();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(X2);
    auto rank = lu_decomp.rank();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp2(X);
    auto rank1 = lu_decomp2.rank();
    result.lls.at(poi_col_idx, 12) = rank1;

    lrs = 2.0 * (ll2 - ll1);
    // lrs_pval = std::abs(chisq(lrs, n_parms2 - result.num_parms));
    lrs_pval = std::abs(chisq(lrs, rank - rank1));


    // calc and set betas
    arma::colvec temp_b = arma::colvec(beta.data(), beta.size(), true, false);
    arma::colvec temp_b2 =
        arma::colvec(beta2.data(), beta2.size(), true, false);
    temp_se = arma::colvec(diag.data(), diag.size(), true, false);
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    temp_se2 = arma::colvec(diag2.data(), diag2.size(), true, false);
    neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;

    result.set_lls(ll1, ll2, lrs, lrs_pval, 3, poi_col_idx, rank);
    pval = (*dist_func_r)(neg_abs_z, df);
    pval2 = (*dist_func_r)(neg_abs_z2, df2);

    result.set_betas_fit1(temp_b, temp_se, pval, poi_col_idx);
    result.set_betas_fit2_sqrd(temp_b2, temp_se2, pval2, poi_col_idx);
  } // #pragma omp parallel for
};

void LogisticRegression::run_vla(arma::mat &cov, arma::mat &pheno,
                                 arma::mat &poi_data, VLAResult &result,
                                 int max_iter, bool is_t_dist,
                                 double maf_thresh) {
  arma::colvec poi_col;
  std::vector<int> poi_2_idx;
  std::vector<int> poi_3_idx;

  // Initialize necessary variables
  double w, maf, one_counter;
  bool ones, twos, temp;

// Parallel region for classifying poi_data into poi_2_idx and poi_3_idx
#pragma omp parallel
  {
    // Thread-local buffers to avoid data races
    std::vector<int> local_poi_2_idx;
    std::vector<int> local_poi_3_idx;

#pragma omp for nowait
    for (arma::uword poi_col_idx = 0; poi_col_idx < poi_data.n_cols;
         poi_col_idx++) {
      arma::colvec w2_col = result.W2.col(poi_col_idx);
      w = arma::accu(w2_col);
      poi_col = poi_data.col(poi_col_idx);
      maf = arma::dot(poi_col.t(), w2_col) / (2.0 * w);
      result.lls.at(poi_col_idx, 6) = maf;

      if ((maf > maf_thresh) & (maf < (1.0 - maf_thresh))) {
        ones = arma::any(poi_col == 1.0);
        twos = arma::any(poi_col == 2.0);
        one_counter = arma::accu(poi_col == 1.0); // Count occurrences of 1.0

        if (ones && twos) {
          local_poi_3_idx.push_back(poi_col_idx);
        } else {
          local_poi_2_idx.push_back(poi_col_idx);
        }

        result.lls.at(poi_col_idx, 7) = one_counter;
      }
    }

// Combine thread-local buffers into the shared vectors
#pragma omp critical
    {
      poi_2_idx.insert(poi_2_idx.end(), local_poi_2_idx.begin(),
                       local_poi_2_idx.end());
      poi_3_idx.insert(poi_3_idx.end(), local_poi_3_idx.begin(),
                       local_poi_3_idx.end());
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
