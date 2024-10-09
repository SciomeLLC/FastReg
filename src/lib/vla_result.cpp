#include <vla_result.h>

VLAResult::VLAResult(arma::mat &covar_matrix, arma::mat &poi_matrix,
                     arma::mat &no_interactions, arma::mat &interactions,
                     arma::mat &interactions_sqrd) {
  num_parms = no_interactions.n_cols + covar_matrix.n_cols;
  num_parms2 = interactions.n_cols + covar_matrix.n_cols;
  num_parms2_sqrd = 2 * interactions.n_cols + covar_matrix.n_cols;

  beta_est = arma::mat(num_parms, poi_matrix.n_cols, arma::fill::zeros);
  se_beta = arma::mat(num_parms, poi_matrix.n_cols, arma::fill::zeros);
  neglog10_pvl = arma::mat(num_parms, poi_matrix.n_cols, arma::fill::zeros);

  beta_est2 = arma::mat(num_parms2, poi_matrix.n_cols, arma::fill::zeros);
  se_beta2 = arma::mat(num_parms2, poi_matrix.n_cols, arma::fill::zeros);
  beta_est2_sqrd =
      arma::mat(num_parms2_sqrd, poi_matrix.n_cols, arma::fill::zeros);
  se_beta2_sqrd =
      arma::mat(num_parms2_sqrd, poi_matrix.n_cols, arma::fill::zeros);
  neglog10_pvl2 = arma::mat(num_parms2, poi_matrix.n_cols, arma::fill::zeros);
  neglog10_pvl2_sqrd =
      arma::mat(num_parms2_sqrd, poi_matrix.n_cols, arma::fill::zeros);

  beta_rel_errs = arma::colvec(poi_matrix.n_cols, arma::fill::zeros);
  beta_abs_errs = arma::colvec(poi_matrix.n_cols, arma::fill::zeros);
  beta_rel_errs2 = arma::colvec(poi_matrix.n_cols, arma::fill::zeros);
  beta_abs_errs2 = arma::colvec(poi_matrix.n_cols, arma::fill::zeros);

  iters = arma::mat(poi_matrix.n_cols, 2, arma::fill::zeros);
  lls = arma::mat(poi_matrix.n_cols, 6, arma::fill::zeros);
  W2 = arma::mat(poi_matrix.n_rows, poi_matrix.n_cols, arma::fill::ones);

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

void VLAResult::set_lls(double ll1, double ll2, double lrs, double lrs_pval,
                        int num_g, int idx, int rank) {
  lls.at(idx, 0) = ll1;
  lls.at(idx, 1) = ll2;
  lls.at(idx, 2) = lrs;
  lls.at(idx, 3) = lrs_pval;
  lls.at(idx, 4) = num_g;
  lls.at(idx, 5) = rank;
}

void VLAResult::set_betas_fit1(arma::colvec &beta, arma::colvec &se,
                               arma::colvec &pval, int idx) {
  beta_est.col(idx) = beta;
  se_beta.col(idx) = se;
  neglog10_pvl.col(idx) = pval;
}

void VLAResult::set_betas_fit2(arma::colvec &beta, arma::colvec &se,
                               arma::colvec &pval, int idx) {
  beta_est2.col(idx) = beta;
  se_beta2.col(idx) = se;
  neglog10_pvl2.col(idx) = pval;
}

void VLAResult::set_betas_fit2_sqrd(arma::colvec &beta, arma::colvec &se,
                                    arma::colvec &pval, int idx) {
  beta_est2_sqrd.col(idx) = beta;
  se_beta2_sqrd.col(idx) = se;
  neglog10_pvl2_sqrd.col(idx) = pval;
}

void VLAResult::write_to_file(std::string dir, std::string file_name) {
  int n_parms = beta_est.n_rows;
  int n_parms2 = beta_est2.n_rows;

  // create dir if it doesn't exist
  fs::create_directory(dir);
  // Rcpp::Rcout << "Dir created or already exists" << std::endl;
  if (srt_cols.size() != beta_est.n_cols) {
    Rcpp::Rcout << "Error: The size of poi_names does not match the number of "
                   "columns in the beta matrix."
                << std::endl;
    // return;
  }

  std::stringstream ss;
  ss << dir << "/" << file_name << ".tsv";
  std::string result_file = ss.str();
  std::ofstream outfile;
  // TODO: figure out column index for interactions
  //   auto tmp_idx = std::find(cov_int_names_sqrd.begin(),
  //   cov_int_names_sqrd.end(), "poi");

  //   if (tmp_idx == cov_int_names_sqrd.end()) {
  //     Rcpp::stop("Couldn't find a poi column");
  //   }
  //   int idx = tmp_idx - cov_int_names_sqrd.begin();
  int idx = beta_est.n_rows + 1;
  if (fs::exists(result_file)) {
    outfile.open(result_file, std::ios::app);
    if (!outfile.is_open()) {
      Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
    }
    // Rcpp::Rcout << "File already exists and opened for writing/appending." <<
    // std::endl;
  } else {
    outfile.open(result_file);
    if (!outfile.is_open()) {
      Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
    }

    outfile << "POI\tN\tDF_fit1\tEstimate_fit1_[poi]\tSE_fit1_[poi]\tmlog10P_"
               "fit1_[poi]\tAbs_"
               "Err_fit1\tRel_Err_fit1\titers_fit1\tDF_fit2";

    // TODO: figure out column names for interactions to write to file
    // estimate, se, pval for poi, poi^2 ineraction terms
    // for (int i = idx; i < cov_int_names_sqrd.size(); i++) {
    //   outfile << "\tEstimate_[" << cov_int_names_sqrd[i] << "]"
    //           << "\tSE_[" << cov_int_names_sqrd[i] << "]"
    //           << "\tmlog10P_[" << cov_int_names_sqrd[i] << "]";
    // }

    outfile << "\tAbs_Err_fit2\tRel_Err_fit2\titers_fit2\tLL_fit1\tLL_"
               "fit2\tLRS\tLRS_mlog10pvl\tnumG\trank"
            << std::endl;
  }

  outfile << std::fixed << std::setprecision(17);

  std::stringstream buffer;
  for (int col = 0; col < (int)beta_est.n_cols; col++) {
    std::string poi_name = file_name + std::to_string(col);
    double abs_err_val = beta_abs_errs.at(col);
    double rel_err_val = beta_rel_errs.at(col);
    double iter1 = iters.at(col, 0);
    double iter2 = iters.at(col, 1);
    double abs_err_val2 = beta_abs_errs2.at(col);
    double rel_err_val2 = beta_rel_errs2.at(col);
    double ll1 = lls.at(col, 0);
    double ll2 = lls.at(col, 1);
    double lrs = lls.at(col, 2);
    double lrs_pval = lls.at(col, 3);
    double num_G = lls.at(col, 4);
    double rank = lls.at(col, 5);

    int N = arma::as_scalar(arma::sum(W2.col(col), 0));
    int df = N - n_parms;
    n_parms2 = (num_G == 3) ? beta_est2_sqrd.n_rows : beta_est2.n_rows;
    int df2 = N - n_parms2;
    int adder2 = 1;
    buffer << poi_name << "\t" << N << "\t" << df << "\t"
           << beta_est.at(idx, col) << "\t" << se_beta.at(idx, col) << "\t"
           << neglog10_pvl.at(idx, col) << "\t" << abs_err_val << "\t"
           << rel_err_val << "\t" << iter1 << "\t" << df2;

    int row2 = idx;
    int row3 = (int)beta_est2.n_rows;
    if (num_G == 2) {
      for (; row2 < (int)beta_est2.n_rows; row2 += adder2) {
        buffer << "\t" << beta_est2.at(row2, col) << "\t"
               << se_beta2.at(row2, col) << "\t" << neglog10_pvl2.at(row2, col);
      }
    } else {
      for (; row2 < (int)beta_est2.n_rows; row2 += adder2) {
        buffer << "\t" << beta_est2_sqrd.at(row2, col) << "\t"
               << se_beta2_sqrd.at(row2, col) << "\t"
               << neglog10_pvl2_sqrd.at(row2, col);
      }
    }

    for (; row3 < (int)beta_est2_sqrd.n_rows; row3 += adder2) {
      buffer << "\t" << beta_est2_sqrd.at(row3, col) << "\t"
             << se_beta2_sqrd.at(row3, col) << "\t"
             << neglog10_pvl2_sqrd.at(row3, col);
    }

    buffer << "\t" << abs_err_val2 << "\t" << rel_err_val2 << "\t" << iter2
           << "\t" << ll1 << "\t" << ll2 << "\t" << lrs << "\t" << lrs_pval
           << "\t" << num_G << "\t" << rank << std::endl;
  }
  outfile << buffer.str();
  outfile.close();
}

void VLAResult::concatenate(std::string output_dir,
                            std::string file_name_prefix,
                            std::string file_concatenation_prefix) {
  std::set<int> stratums;
  std::map<int, std::set<std::string>> stratum_files;
  for (const auto &entry : fs::directory_iterator(output_dir)) {
    if (entry.path().extension() == ".tsv" &&
        entry.path().filename().string().find(file_name_prefix) !=
            std::string::npos) {
      // Extract stratum from filename
      std::string filename = entry.path().filename().stem().string();
      std::vector<std::string> tokens;
      std::stringstream ss(filename);
      std::string token;
      while (std::getline(ss, token, '_')) {
        tokens.push_back(token);
      }
      if (tokens.size() >= 2) {
        try {
          int stratum =
              std::stoi(tokens[tokens.size() - 2]); // Second-to-last token
          stratums.insert(stratum);
          stratum_files[stratum].insert(filename);
        } catch (const std::invalid_argument &ex) {
          Rcpp::Rcerr << "Invalid stratum found in filename: " << filename
                      << std::endl;
        }
      }
    }
  }

  // Concatenate files for each unique stratum
  for (int stratum : stratums) {
    std::string outputFile = output_dir + "/" + file_concatenation_prefix +
                             "_" + file_name_prefix + "_stratum_" +
                             std::to_string(stratum) + ".tsv";

    // Check if the output file already exists
    if (fs::exists(outputFile)) {
      Rcpp::Rcerr << "Output file already exists: " << outputFile
                  << ". It will be overwritten.\n";
    }

    std::ofstream out(outputFile);

    if (!out.is_open()) {
      Rcpp::Rcerr << "Failed to open the output file: " << outputFile
                  << std::endl;
      continue;
    }
    bool header_written = false;
    for (const auto &filename : stratum_files[stratum]) {
      fs::path file_path = output_dir + "/" + filename + ".tsv";
      std::ifstream in(file_path.string());

      if (!in.is_open()) {
        Rcpp::Rcerr << "Failed to open " << file_path << std::endl;
        continue;
      }
      // Skip header line if it has already been written
      if (header_written) {
        std::string headerLine;
        std::getline(in, headerLine);
      } else {
        std::string headerLine;
        std::getline(in, headerLine);
        out << headerLine << std::endl;
        header_written = true;
      }
      out << in.rdbuf();
      in.close();
      fs::remove(file_path);
    }
    out.close();
  }
}