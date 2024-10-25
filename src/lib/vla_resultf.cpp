#include <vla_result.h>

VLAResultf::VLAResultf(arma::fmat &covar_matrix, arma::fmat &poi_matrix,
                       arma::fmat &no_interactions, arma::fmat &interactions,
                       arma::fmat &interactions_sqrd, std::string lt) {
  local_time = lt;
  cov_n_cols = covar_matrix.n_cols;
  num_parms = no_interactions.n_cols + covar_matrix.n_cols;
  num_parms2 = interactions.n_cols + covar_matrix.n_cols;
  num_parms2_sqrd = 2 * interactions.n_cols + covar_matrix.n_cols;

  beta_est = arma::fmat(num_parms, poi_matrix.n_cols, arma::fill::zeros);
  se_beta = arma::fmat(num_parms, poi_matrix.n_cols, arma::fill::zeros);
  neglog10_pvl = arma::fmat(num_parms, poi_matrix.n_cols, arma::fill::zeros);

  beta_est2 = arma::fmat(num_parms2, poi_matrix.n_cols, arma::fill::zeros);
  se_beta2 = arma::fmat(num_parms2, poi_matrix.n_cols, arma::fill::zeros);
  beta_est2_sqrd =
      arma::fmat(num_parms2_sqrd, poi_matrix.n_cols, arma::fill::zeros);
  se_beta2_sqrd =
      arma::fmat(num_parms2_sqrd, poi_matrix.n_cols, arma::fill::zeros);
  neglog10_pvl2 = arma::fmat(num_parms2, poi_matrix.n_cols, arma::fill::zeros);
  neglog10_pvl2_sqrd =
      arma::fmat(num_parms2_sqrd, poi_matrix.n_cols, arma::fill::zeros);

  beta_rel_errs = arma::fcolvec(poi_matrix.n_cols, arma::fill::zeros);
  beta_abs_errs = arma::fcolvec(poi_matrix.n_cols, arma::fill::zeros);
  beta_rel_errs2 = arma::fcolvec(poi_matrix.n_cols, arma::fill::zeros);
  beta_abs_errs2 = arma::fcolvec(poi_matrix.n_cols, arma::fill::zeros);

  beta_est.fill(arma::datum::nan);
  beta_est2.fill(arma::datum::nan);
  se_beta.fill(arma::datum::nan);
  se_beta2.fill(arma::datum::nan);
  beta_est2_sqrd.fill(arma::datum::nan);
  se_beta2_sqrd.fill(arma::datum::nan);
  neglog10_pvl2.fill(arma::datum::nan);
  neglog10_pvl2_sqrd.fill(arma::datum::nan);
  beta_rel_errs.fill(arma::datum::nan);
  beta_rel_errs2.fill(arma::datum::nan);
  beta_abs_errs.fill(arma::datum::nan);
  beta_abs_errs2.fill(arma::datum::nan);

  iters = arma::fmat(poi_matrix.n_cols, 2, arma::fill::zeros);
  lls = arma::fmat(poi_matrix.n_cols, 8, arma::fill::zeros);
  iters.fill(arma::datum::nan);
  lls.fill(arma::datum::nan);
  W2 = arma::fmat(poi_matrix.n_rows, poi_matrix.n_cols, arma::fill::ones);

  cov_no_int_names.resize(num_parms);
  cov_int_names.resize(num_parms2);
  cov_int_names_sqrd.resize(num_parms2_sqrd);

  // TODO: Create and set column names and row names for results
  //  write_to_file() will require the col/row names to output to file correctly

  for (arma::uword v = 0; v < poi_matrix.n_cols; v++) {
    arma::uvec G_na = arma::find_nonfinite(poi_matrix.col(v));
    for (arma::uword i = 0; i < G_na.n_elem; i++) {
      W2(G_na(i), v) = 0.0;
      poi_matrix(G_na(i), v) = 0.0;
    }
  }

  poi_sqrd_mat = arma::square(poi_matrix);
  this->interactions = &interactions;
  this->no_interactions = &no_interactions;
  this->interactions_sqrd = &interactions_sqrd;
}

void VLAResultf::set_lls(float ll1, float ll2, float lrs, float lrs_pval,
                         int num_g, int idx, int rank) {
  lls.at(idx, 0) = ll1;
  lls.at(idx, 1) = ll2;
  lls.at(idx, 2) = lrs;
  lls.at(idx, 3) = lrs_pval;
  lls.at(idx, 4) = num_g;
  lls.at(idx, 5) = rank;
}

void VLAResultf::set_betas_fit1(arma::fcolvec &beta, arma::fcolvec &se,
                                arma::fcolvec &pval, int idx) {
  beta_est.col(idx) = beta;
  se_beta.col(idx) = se;
  neglog10_pvl.col(idx) = pval;
}

void VLAResultf::set_betas_fit2(arma::fcolvec &beta, arma::fcolvec &se,
                                arma::fcolvec &pval, int idx) {
  beta_est2.col(idx) = beta;
  se_beta2.col(idx) = se;
  neglog10_pvl2.col(idx) = pval;
}

void VLAResultf::set_betas_fit2_sqrd(arma::fcolvec &beta, arma::fcolvec &se,
                                     arma::fcolvec &pval, int idx) {
  beta_est2_sqrd.col(idx) = beta;
  se_beta2_sqrd.col(idx) = se;
  neglog10_pvl2_sqrd.col(idx) = pval;
}

void VLAResultf::write_to_file(std::string dir, std::string file_name,
                               std::string pheno_name,
                               std::vector<std::string> row_names) {
  // create dir if it doesn't exist
  fs::create_directory(dir);
  fs::create_directory(dir + "/" + pheno_name);

  int n_parms = beta_est.n_rows;
  int n_parms2 = beta_est2.n_rows;

  std::stringstream ss;
  ss << dir << "/" << pheno_name << "/" << file_name << "_results_"
     << getLocalTime() << ".tsv";
  std::string result_file = ss.str();
  std::ofstream outfile;

  int idx = cov_n_cols + 1;
  int G_idx = cov_n_cols;
  int Gsq_idx = beta_est2.n_rows;
  int poi_idx = beta_est.n_rows - 1;
  if (fs::exists(result_file)) {
    outfile.open(result_file, std::ios::app);
    if (!outfile.is_open()) {
      Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
    }
  } else {
    outfile.open(result_file);
    if (!outfile.is_open()) {
      Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
    }
  }

  outfile << std::fixed << std::setprecision(17);

  std::stringstream buffer;
  float abs_err_val, rel_err_val, iter1, iter2, abs_err_val2, rel_err_val2, ll1,
      ll2, lrs, lrs_pval, num_G, rank;
  int N, df, df2, row;
  for (int col = 0; col < (int)beta_est.n_cols; col++) {
    abs_err_val = beta_abs_errs.at(col);
    rel_err_val = beta_rel_errs.at(col);
    iter1 = iters.at(col, 0);
    iter2 = iters.at(col, 1);
    abs_err_val2 = beta_abs_errs2.at(col);
    rel_err_val2 = beta_rel_errs2.at(col);
    ll1 = lls.at(col, 0);
    ll2 = lls.at(col, 1);
    lrs = lls.at(col, 2);
    lrs_pval = lls.at(col, 3);
    num_G = lls.at(col, 4);
    rank = lls.at(col, 5);

    N = arma::as_scalar(arma::sum(W2.col(col), 0));
    df = N - n_parms;
    n_parms2 = (num_G == 3) ? beta_est2_sqrd.n_rows : beta_est2.n_rows;
    df2 = N - n_parms2;
    buffer << row_names[col] << "\t" << N << "\t" << df << "\t"
           << beta_est.at(poi_idx, col) << "\t" << se_beta.at(poi_idx, col) << "\t"
           << neglog10_pvl.at(poi_idx, col) << "\t" << abs_err_val << "\t"
           << rel_err_val << "\t" << iter1 << "\t" << df2;

    row = idx;
    if (num_G == 2) {
      buffer << "\t" << beta_est2.at(row - 1, col) << "\t"
             << se_beta2.at(row - 1, col) << "\t"
             << neglog10_pvl2.at(row - 1, col);
      while (row < (int)beta_est2.n_rows) {
        buffer << "\t" << beta_est2.at(row, col) << "\t"
               << se_beta2.at(row, col) << "\t" << neglog10_pvl2.at(row, col);
        row++;
      }
    } else {
      buffer << "\t" << beta_est2_sqrd.at(row - 1, col) << "\t"
             << se_beta2_sqrd.at(row - 1, col) << "\t"
             << neglog10_pvl2_sqrd.at(row - 1, col);
      while (row < (int)beta_est2.n_rows) {
        buffer << "\t" << beta_est2_sqrd.at(row, col) << "\t"
               << se_beta2_sqrd.at(row, col) << "\t"
               << neglog10_pvl2_sqrd.at(row, col);
        row++;
      }
    }

    while (row < (int)beta_est2_sqrd.n_rows) {
      buffer << "\t" << beta_est2_sqrd.at(row, col) << "\t"
             << se_beta2_sqrd.at(row, col) << "\t"
             << neglog10_pvl2_sqrd.at(row, col);
      row++;
    }

    buffer << "\t" << abs_err_val2 << "\t" << rel_err_val2 << "\t" << iter2
           << "\t" << ll1 << "\t" << ll2 << "\t" << lrs << "\t" << lrs_pval
           << "\t" << num_G << "\t" << rank << "\t" << lls.at(col, 6) << "\t"
           << lls.at(col, 7) << std::endl;
  }
  outfile << buffer.str();
  outfile.close();
}

void VLAResultf::write_headers(std::string dir, std::string file_name,
                               std::string pheno_name,
                               std::vector<std::string> row_names) {
  fs::create_directory(dir);
  fs::create_directory(dir + "/" + pheno_name);

  std::stringstream ss;
  ss << dir << "/" << pheno_name << "/" << file_name << "_headers_"
     << getLocalTime() << ".tsv";
  std::string result_file = ss.str();
  std::ofstream outfile;

  outfile.open(result_file);
  if (!outfile.is_open()) {
    Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
  }

  int idx = cov_n_cols + 1;
  int row = idx;

  std::stringstream buffer;
  buffer << "vid"
         << "\t"
         << "N"
         << "\t"
         << "DF_fit1"
         << "\t"
         << "Est_fit1"
         << "\t"
         << "StdErr_fit1"
         << "\t"
         << "mlog10_Pval_fit1"
         << "\t"
         << "AbsErr_fit1"
         << "\t"
         << "RelErr_fit1"
         << "\t"
         << "iter_fit1"
         << "\t"
         << "DF_fit2"
         << "\t"
         << "EstG_fit2"
         << "\t"
         << "StdErrG_fit2"
         << "\t"
         << "mlog10_PvalG_fit2";
  while (row < (int)beta_est2.n_rows) {
    buffer << "\tEstG*X" << row - idx << "_fit2"
           << "\tStdErrG*X" << row - idx << "_fit2"
           << "\tmlog10_PvalG*X" << row - idx << "_fit2";
    row++;
  }
  buffer << "\t"
         << "EstGsq_fit2"
         << "\t"
         << "StdErrGsq_fit2"
         << "\t"
         << "mlog10_PvalGsq_fit2";
  row++;
  int sqrd_idx = row;
  while (row < (int)beta_est2_sqrd.n_rows) {
    buffer << "\tEstGsq*X" << (row - sqrd_idx) << "_fit2"
           << "\tStdErrGsq*X" << (row - sqrd_idx) << "_fit2"
           << "\tmlog10_PvalGsq*X" << (row - sqrd_idx) << "_fit2";
    row++;
  }
  buffer << "\t"
         << "AbsErr_fit2"
         << "\t"
         << "RelErr_fit2"
         << "\t"
         << "iter_fit2"
         << "\t"
         << "LogLikelihood_fit1"
         << "\t"
         << "LogLikelihood_fit2"
         << "\t"
         << "LRS"
         << "\t"
         << "LRS_Pval"
         << "\t"
         << "num_G"
         << "\t"
         << "rank"
         << "\t"
         << "caf"
         << "\t"
         << "het_count" << std::endl;
  outfile << buffer.str();
  outfile.close();
}