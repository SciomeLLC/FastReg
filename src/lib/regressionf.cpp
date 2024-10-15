#include <regression.h>

arma::fcolvec t_dist_r(arma::fcolvec abs_z, int df) {
  arma::fcolvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) = -1 * (R::pt(abs_z(i), df, true, true) + log(2)) / log(10);
  }

  return ret_val;
}

float chisq(float lrs, int df) {
  return (float)R::pchisq(lrs, df, 0, 1) / (-1.0 * log(10));
}

arma::fcolvec norm_dist_r(arma::fcolvec abs_z, int df) {
  arma::fcolvec ret_val(abs_z.size());
  for (size_t i = 0; i < abs_z.size(); i++) {
    ret_val(i) =
        -1 * (R::pnorm(abs_z(i), 1.0, 1.0, true, true) + log(2)) / log(10);
  }
  return ret_val;
}

void LogisticRegression::run_vla_2(arma::fmat &cov, arma::fmat &pheno, arma::fmat &poi_data,
               VLAResultf &result, int max_iter, bool is_t_dist,
               std::vector<int> &poi_2_idx, Eigen::MatrixXf &W2f,
               Eigen::MatrixXf &tphenoD) {
  arma::fcolvec (*dist_func_r)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  double ll1, ll2, lrs, lrs_pval;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::fcolvec)poi_data.col(0)).size();
  arma::fcolvec temp_se(poi_col_size, arma::fill::zeros);
  // arma::fcolvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::fcolvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::fcolvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::fmat poi_col(poi_col_size, 1, arma::fill::zeros);
  Eigen::MatrixXf cov_w_mat =
      Eigen::Map<Eigen::MatrixXf>(cov.memptr(), cov.n_rows, cov.n_cols);
  Eigen::MatrixXf int_w_mat, int_w_mat2;
  arma::uword n_parms2 = result.cov_int_names.size();
  // #pragma omp parallel for
  for (int poi_col_idx : poi_2_idx) {
    checkInterrupt();
    Eigen::MatrixXf A =
        Eigen::MatrixXf::Zero(result.num_parms, result.num_parms);
    A(0, 0) = 1;
    A(result.num_parms, result.num_parms) = 10;
    Eigen::VectorXf beta = Eigen::VectorXf::Zero(result.num_parms, 1);
    Eigen::VectorXf beta_old = beta;
    cov_w_mat =
        Eigen::Map<Eigen::MatrixXf>(cov.memptr(), cov.n_rows, cov.n_cols);
    int_w_mat = Eigen::Map<Eigen::MatrixXf>(result.no_interactions->memptr(),
                                            result.no_interactions->n_rows,
                                            result.no_interactions->n_cols);
    Eigen::VectorXf w2_col = W2f.col(poi_col_idx);
    poi_col = poi_data.col(poi_col_idx);
    Eigen::MatrixXf POI =
        Eigen::Map<Eigen::MatrixXf>(poi_col.memptr(), poi_col.n_rows, 1);
    Eigen::MatrixXf temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXf X =
        Eigen::MatrixXf::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXf beta_diff = (beta - beta_old).array().abs();
    result.beta_rel_errs.at(poi_col_idx) = 1e8;

    // Fit 2
    Eigen::MatrixXf A2 = Eigen::MatrixXf::Zero(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    Eigen::VectorXf beta2 = Eigen::VectorXf::Zero(n_parms2, 1);
    Eigen::VectorXf beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXf>(result.interactions->memptr(),
                                             result.interactions->n_rows,
                                             result.interactions->n_cols);
    Eigen::VectorXf w2_col2 = W2f.col(poi_col_idx);
    Eigen::MatrixXf temp_mat_e2 = int_w_mat2;
    int_w_mat2 = temp_mat_e2.array().colwise() * POI.col(0).array();

    Eigen::MatrixXf X2 = Eigen::MatrixXf::Zero(
        int_w_mat2.rows(), cov_w_mat.cols() + int_w_mat2.cols());
    X2.topLeftCorner(int_w_mat2.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_mat2.rows(), int_w_mat2.cols()) = int_w_mat2;

    Eigen::VectorXf beta_diff2 = (beta2 - beta_old2).array().abs();
    result.beta_rel_errs2.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    Eigen::MatrixXf eta, eta2, p, p2;

    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs.at(poi_col_idx) > 1e-4;
         iter++) {
      // Fit 1
      eta = Eigen::MatrixXf::Zero(cov.n_rows, 1);
      eta = X * beta;
      p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();
      Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;
      Eigen::MatrixXf z =
          w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXf B = X.transpose() * z;
      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      result.beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      result.iters.at(poi_col_idx, 0) = iter + 1;
      Eigen::MatrixXf mTemp = (beta_diff.array() / beta_old.array());
      result.beta_rel_errs.at(poi_col_idx) = mTemp.maxCoeff();
      beta_old = beta;
    }

    // Fit 2
    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs2.at(poi_col_idx) > 1e-10;
         iter++) {
      eta2 = Eigen::MatrixXf::Zero(cov.n_rows, 1);
      eta2 = X2 * beta2;
      p2 = -eta2.array();
      p2 = p2.array().exp() + 1;
      p2 = 1 / p2.array();
      Eigen::MatrixXf W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      Eigen::MatrixXf temp1_2 = X2.array().colwise() * W1_2.col(0).array();
      A2 = temp1_2.transpose() * X2;
      Eigen::MatrixXf z2 =
          w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXf B2 = X2.transpose() * z2;
      beta2 = beta2 +
              A2.ldlt().solve(B2); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff2 = (beta2.array() - beta_old2.array()).abs();
      result.beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      Eigen::MatrixXf mTemp2 = (beta_diff2.array() / beta_old2.array());
      result.iters.at(poi_col_idx, 1) = iter + 1;
      result.beta_rel_errs2.at(poi_col_idx) = mTemp2.maxCoeff();
      beta_old2 = beta2;
    }

    int df = w2_col.array().sum() - result.num_parms;
    Eigen::VectorXf diag = A.inverse().diagonal().array().sqrt();

    int df2 = w2_col.array().sum() - n_parms2;
    Eigen::VectorXf diag2 = A2.inverse().diagonal().array().sqrt();

    // calculate LLs
    eta = X * beta;
    p = -eta.array();
    p = p.array().exp() + 1;
    p = 1 / p.array();
    eta2 = X2 * beta2;
    p2 = -eta2.array();
    p2 = p2.array().exp() + 1;
    p2 = 1 / p2.array();
    arma::fcolvec p_arma = arma::fcolvec(p.data(), p.size(), true, false);
    arma::fcolvec p2_arma = arma::fcolvec(p2.data(), p2.size(), true, false);
    arma::fcolvec eta_arma = arma::fcolvec(eta.data(), eta.size(), true, false);
    arma::fcolvec eta2_arma =
        arma::fcolvec(eta2.data(), eta2.size(), true, false);

    p_arma.elem(arma::find_nonfinite(p_arma)).zeros();
    p_arma.elem(arma::find(p_arma >= 1.0)).fill(1.0 - 1e-4);
    ll1 = arma::accu((log(1.0 - p_arma) + pheno.col(0) % eta_arma) %
                     result.W2.col(poi_col_idx));

    p2_arma.elem(arma::find_nonfinite(p2_arma)).zeros();
    p2_arma.elem(arma::find(p2_arma >= 1.0)).fill(1.0 - 1e-4);
    eta2_arma.elem(arma::find_nonfinite(eta2_arma)).zeros();
    ll2 = arma::accu((log(1.0 - p2_arma) + pheno.col(0) % eta2_arma) %
                     result.W2.col(poi_col_idx));
    lrs = 2.0 * (ll2 - ll1);
    lrs_pval = std::abs(chisq(lrs, n_parms2 - result.num_parms));

    // calc and set betas
    temp_se = arma::fcolvec(diag.data(), diag.size(), true, false);
    arma::fcolvec temp_b = arma::fcolvec(beta.data(), beta.size(), true, false);
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    arma::fcolvec temp_se2 =
        arma::fcolvec(diag2.data(), diag2.size(), true, false);
    arma::fcolvec temp_b2 =
        arma::fcolvec(beta2.data(), beta2.size(), true, false);
    neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;

    Eigen::FullPivLU<Eigen::MatrixXf> lu_decomp(X2);
    auto rank = lu_decomp.rank();
    result.set_lls(ll1, ll2, lrs, lrs_pval, 2, poi_col_idx, rank);
    arma::fcolvec pval = (*dist_func_r)(neg_abs_z, df);
    arma::fcolvec pval2 = (*dist_func_r)(neg_abs_z2, df2);
    result.set_betas_fit1(temp_b, temp_se, pval, poi_col_idx);
    result.set_betas_fit2(temp_b2, temp_se2, pval2, poi_col_idx);
  } // #pragma omp parallel for
};

void LogisticRegression::run_vla_3(arma::fmat &cov, arma::fmat &pheno, arma::fmat &poi_data,
               VLAResultf &result, int max_iter, bool is_t_dist,
               std::vector<int> &poi_3_idx, Eigen::MatrixXf &W2f,
               Eigen::MatrixXf &tphenoD) {
  arma::fcolvec (*dist_func_r)(arma::fcolvec, int) =
      is_t_dist == true ? t_dist_r : norm_dist_r;
  double ll1, ll2, lrs, lrs_pval;
  int df, df2, lrs_df;
  int poi_col_size = ((arma::fcolvec)poi_data.col(0)).size();
  arma::fcolvec temp_se(poi_col_size, arma::fill::zeros);
  // arma::fcolvec temp_se2(poi_col_size, arma::fill::zeros);
  arma::fcolvec neg_abs_z(poi_col_size, arma::fill::zeros);
  arma::fcolvec neg_abs_z2(poi_col_size, arma::fill::zeros);
  arma::fmat poi_col(poi_col_size, 1, arma::fill::zeros);
  arma::fmat poi_col_sqrd(poi_col_size, 1, arma::fill::zeros);
  Eigen::MatrixXf cov_w_mat =
      Eigen::Map<Eigen::MatrixXf>(cov.memptr(), cov.n_rows, cov.n_cols);
  Eigen::MatrixXf int_w_mat, int_w_mat2, int_w_mat2_sqrd;
  arma::uword n_parms2 = result.cov_int_names_sqrd.size();
  // #pragma omp parallel for
  for (int poi_col_idx : poi_3_idx) {
    checkInterrupt();
    Eigen::MatrixXf A =
        Eigen::MatrixXf::Zero(result.num_parms, result.num_parms);
    A(0, 0) = 1;
    A(result.num_parms, result.num_parms) = 10;
    Eigen::VectorXf beta = Eigen::VectorXf::Zero(result.num_parms, 1);
    Eigen::VectorXf beta_old = beta;
    cov_w_mat =
        Eigen::Map<Eigen::MatrixXf>(cov.memptr(), cov.n_rows, cov.n_cols);
    int_w_mat = Eigen::Map<Eigen::MatrixXf>(result.no_interactions->memptr(),
                                            result.no_interactions->n_rows,
                                            result.no_interactions->n_cols);
    Eigen::VectorXf w2_col = W2f.col(poi_col_idx);
    poi_col = poi_data.col(poi_col_idx);

    Eigen::MatrixXf POI =
        Eigen::Map<Eigen::MatrixXf>(poi_col.memptr(), poi_col.n_rows, 1);
    Eigen::MatrixXf temp_mat_e = int_w_mat;
    int_w_mat = temp_mat_e.array().colwise() * POI.col(0).array();
    Eigen::MatrixXf X =
        Eigen::MatrixXf::Zero(int_w_mat.rows(), cov_w_mat.cols() + 1);
    X.topLeftCorner(int_w_mat.rows(), cov_w_mat.cols()) = cov_w_mat;
    X.topRightCorner(int_w_mat.rows(), int_w_mat.cols()) = int_w_mat;

    Eigen::VectorXf beta_diff = (beta - beta_old).array().abs();
    result.beta_rel_errs.at(poi_col_idx) = 1e8;

    // Fit 2
    poi_col_sqrd = result.poi_sqrd_mat.col(poi_col_idx);
    Eigen::MatrixXf POI_sqrd = Eigen::Map<Eigen::MatrixXf>(
        poi_col_sqrd.memptr(), poi_col_sqrd.n_rows, 1);
    Eigen::MatrixXf A2 = Eigen::MatrixXf::Zero(n_parms2, n_parms2);
    A2(0, 0) = 1;
    A2(n_parms2, n_parms2) = 10;
    Eigen::VectorXf beta2 = Eigen::VectorXf::Zero(n_parms2, 1);
    Eigen::VectorXf beta_old2 = beta2;
    int_w_mat2 = Eigen::Map<Eigen::MatrixXf>(result.interactions->memptr(),
                                             result.interactions->n_rows,
                                             result.interactions->n_cols);
    int_w_mat2_sqrd = Eigen::Map<Eigen::MatrixXf>(
        result.interactions_sqrd->memptr(), result.interactions_sqrd->n_rows,
        result.interactions_sqrd->n_cols);
    Eigen::VectorXf w2_col2 = W2f.col(poi_col_idx);
    Eigen::MatrixXf temp_mat_e2 = int_w_mat2;
    Eigen::MatrixXf temp_mat_e2_sqrd = int_w_mat2_sqrd;
    int_w_mat2 = temp_mat_e2.array().colwise() * POI.col(0).array();
    int_w_mat2_sqrd =
        temp_mat_e2_sqrd.array().colwise() * POI_sqrd.col(0).array();

    Eigen::MatrixXf int_w_sqrd(int_w_mat2.rows(),
                               int_w_mat2.cols() + int_w_mat2_sqrd.cols());
    int_w_sqrd << int_w_mat2, int_w_mat2_sqrd;
    Eigen::MatrixXf X2 = Eigen::MatrixXf::Zero(
        int_w_mat2.rows(), cov_w_mat.cols() + int_w_sqrd.cols());
    X2.topLeftCorner(int_w_sqrd.rows(), cov_w_mat.cols()) = cov_w_mat;
    X2.topRightCorner(int_w_sqrd.rows(), int_w_sqrd.cols()) = int_w_sqrd;

    Eigen::VectorXf beta_diff2 = (beta2 - beta_old2).array().abs();
    result.beta_rel_errs2.at(poi_col_idx) = 1e8;

    // initialize eta and p for fit 1/2
    Eigen::MatrixXf eta, eta2, p, p2;

    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs.at(poi_col_idx) > 1e-4;
         iter++) {
      // Fit 1
      eta = Eigen::MatrixXf::Zero(cov.n_rows, 1);
      eta = X * beta;
      p = -eta.array();
      p = p.array().exp() + 1;
      p = 1 / p.array();
      Eigen::MatrixXf W1 = p.array() * (1 - p.array()) * w2_col.array();
      Eigen::MatrixXf temp1 = X.array().colwise() * W1.col(0).array();
      A = temp1.transpose() * X;
      Eigen::MatrixXf z =
          w2_col.array() * (tphenoD.array() - p.array()).array();
      Eigen::MatrixXf B = X.transpose() * z;
      beta = beta +
             A.ldlt().solve(B); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff = (beta.array() - beta_old.array()).abs();
      result.beta_abs_errs.at(poi_col_idx) = beta_diff.array().maxCoeff();
      result.iters.at(poi_col_idx, 0) = iter + 1;
      Eigen::MatrixXf mTemp = (beta_diff.array() / beta_old.array());
      result.beta_rel_errs.at(poi_col_idx) = mTemp.maxCoeff();
      beta_old = beta;
    }

    // Fit 2
    for (int iter = 0;
         iter < max_iter && result.beta_rel_errs2.at(poi_col_idx) > 1e-10;
         iter++) {
      eta2 = Eigen::MatrixXf::Zero(cov.n_rows, 1);
      eta2 = X2 * beta2;
      p2 = -eta2.array();
      p2 = p2.array().exp() + 1;
      p2 = 1 / p2.array();
      Eigen::MatrixXf W1_2 = p2.array() * (1 - p2.array()) * w2_col2.array();
      Eigen::MatrixXf temp1_2 = X2.array().colwise() * W1_2.col(0).array();
      A2 = temp1_2.transpose() * X2;
      Eigen::MatrixXf z2 =
          w2_col2.array() * (tphenoD.array() - p2.array()).array();
      Eigen::MatrixXf B2 = X2.transpose() * z2;
      beta2 = beta2 +
              A2.ldlt().solve(B2); // arma::solve(A, B, arma::solve_opts::fast);

      beta_diff2 = (beta2.array() - beta_old2.array()).abs();
      result.beta_abs_errs2.at(poi_col_idx) = beta_diff2.array().maxCoeff();
      Eigen::MatrixXf mTemp2 = (beta_diff2.array() / beta_old2.array());
      result.iters.at(poi_col_idx, 1) = iter + 1;
      result.beta_rel_errs2.at(poi_col_idx) = mTemp2.maxCoeff();
      beta_old2 = beta2;
    }

    int df = w2_col.array().sum() - result.num_parms;
    Eigen::VectorXf diag = A.inverse().diagonal().array().sqrt();

    int df2 = w2_col.array().sum() - n_parms2;
    Eigen::VectorXf diag2 = A2.inverse().diagonal().array().sqrt();

    // calculate LLs
    eta = X * beta;
    p = -eta.array();
    p = p.array().exp() + 1;
    p = 1 / p.array();
    eta2 = X2 * beta2;

    Eigen::FullPivLU<Eigen::MatrixXf> lu_decomp(X2);
    auto rank = lu_decomp.rank();
    p2 = -eta2.array();
    p2 = p2.array().exp() + 1;
    p2 = 1 / p2.array();
    arma::fcolvec p_arma = arma::fcolvec(p.data(), p.size(), true, false);
    arma::fcolvec p2_arma = arma::fcolvec(p2.data(), p2.size(), true, false);
    arma::fcolvec eta_arma = arma::fcolvec(eta.data(), eta.size(), true, false);
    arma::fcolvec eta2_arma =
        arma::fcolvec(eta2.data(), eta2.size(), true, false);

    p_arma.elem(arma::find_nonfinite(p_arma)).zeros();
    p_arma.elem(arma::find(p_arma >= 1.0)).fill(1.0 - 1e-4);
    ll1 = arma::accu((log(1.0 - p_arma) + pheno.col(0) % eta_arma) %
                     result.W2.col(poi_col_idx));

    p2_arma.elem(arma::find_nonfinite(p2_arma)).zeros();
    p2_arma.elem(arma::find(p2_arma >= 1.0)).fill(1.0 - 1e-4);
    eta2_arma.elem(arma::find_nonfinite(eta2_arma)).zeros();
    ll2 = arma::accu((log(1.0 - p2_arma) + pheno.col(0) % eta2_arma) %
                     result.W2.col(poi_col_idx));
    lrs = 2.0 * (ll2 - ll1);
    lrs_pval = std::abs(chisq(lrs, n_parms2 - result.num_parms));
    // calc and set betas
    temp_se = arma::fcolvec(diag.data(), diag.size(), true, false);
    arma::fcolvec temp_b = arma::fcolvec(beta.data(), beta.size(), true, false);
    neg_abs_z = arma::abs(temp_b / temp_se) * -1;
    arma::fcolvec temp_se2 =
        arma::fcolvec(diag2.data(), diag2.size(), true, false);
    arma::fcolvec temp_b2 =
        arma::fcolvec(beta2.data(), beta2.size(), true, false);
    neg_abs_z2 = arma::abs(temp_b2 / temp_se2) * -1;

    result.set_lls(ll1, ll2, lrs, lrs_pval, 3, poi_col_idx, rank);
    arma::fcolvec pval = (*dist_func_r)(neg_abs_z, df);
    arma::fcolvec pval2 = (*dist_func_r)(neg_abs_z2, df2);
    result.set_betas_fit1(temp_b, temp_se, pval, poi_col_idx);
    result.set_betas_fit2_sqrd(temp_b2, temp_se2, pval2, poi_col_idx);
  } // #pragma omp parallel for
};

void LogisticRegression::run_vla(arma::fmat &cov, arma::fmat &pheno,
                                 arma::fmat &poi_data, VLAResultf &result,
                                 int max_iter, bool is_t_dist,
                                 double maf_thresh) {
  arma::fcolvec poi_col;
  std::vector<int> poi_2_idx;
  std::vector<int> poi_3_idx;

  double w, maf, one_counter;
  bool ones, twos, temp;

  // #pragma omp parallel for
  for (arma::uword poi_col_idx = 0; poi_col_idx < poi_data.n_cols;
       poi_col_idx++) {
    arma::fcolvec w2_col = result.W2.col(poi_col_idx);
    w = arma::accu(w2_col);
    poi_col = poi_data.col(poi_col_idx);
    maf = arma::dot(poi_col.t(), w2_col) / (2.0 * w);
    // Changed 10/10 to CAF from MAF
    // maf = std::min(maf, 1.0 - maf);
    result.lls.at(poi_col_idx, 6) = maf;
    if ((maf > maf_thresh) & (maf < (1.0 - maf_thresh))) {
      // bool has_3_unique = ((arma::fcolvec) arma::unique(poi_col)).size() ==
      // 3;
      ones = false;
      twos = false;
      one_counter = 0.0;
      for (const auto &value : poi_col) {
        temp = (value == 1.0);
        ones = ones | temp;
        one_counter += 1.0 * (temp);
        twos = twos | (value == 2.0);
      }
      if (ones & twos) {
        poi_3_idx.push_back(poi_col_idx);
      } else {
        poi_2_idx.push_back(poi_col_idx);
      }
      result.lls.at(poi_col_idx, 7) = one_counter;
    }
  }

  Eigen::MatrixXf tphenoD =
      Eigen::Map<Eigen::MatrixXf>(pheno.memptr(), pheno.n_rows, pheno.n_cols);
  Eigen::MatrixXf W2f = Eigen::Map<Eigen::MatrixXf>(
      result.W2.memptr(), result.W2.n_rows, result.W2.n_cols);

  // run regression on G with 2 values
  run_vla_2(cov, pheno, poi_data, result, max_iter, is_t_dist, poi_2_idx, W2f,
            tphenoD);
  run_vla_3(cov, pheno, poi_data, result, max_iter, is_t_dist, poi_3_idx, W2f,
            tphenoD);
}