#include <fr_matrix.h>
#include <fr_result.h>

FRResult::FRResult(FRMatrix &covar_matrix, FRMatrix &poi_matrix,
                   FRMatrix &no_interactions, FRMatrix &interactions)
{
  num_parms = no_interactions.data.n_cols + covar_matrix.data.n_cols;
  num_parms2 = interactions.data.n_cols + covar_matrix.data.n_cols;
  num_parms2_sqrd = 2 * interactions.data.n_cols + covar_matrix.data.n_cols;

  beta_est = arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
  se_beta = arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
  neglog10_pvl =
      arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);

  beta_est2 = arma::fmat(num_parms2, poi_matrix.data.n_cols, arma::fill::zeros);
  se_beta2 = arma::fmat(num_parms2, poi_matrix.data.n_cols, arma::fill::zeros);
  beta_est2_sqrd = arma::fmat(num_parms2_sqrd, poi_matrix.data.n_cols, arma::fill::zeros);
  se_beta2_sqrd = arma::fmat(num_parms2_sqrd, poi_matrix.data.n_cols, arma::fill::zeros);
  neglog10_pvl2 =
      arma::fmat(num_parms2, poi_matrix.data.n_cols, arma::fill::zeros);
  neglog10_pvl2_sqrd =
      arma::fmat(num_parms2_sqrd, poi_matrix.data.n_cols, arma::fill::zeros);

  beta_rel_errs = arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);
  beta_abs_errs = arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);
  beta_rel_errs2 = arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);
  beta_abs_errs2 = arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);

  iters = arma::fmat(poi_matrix.data.n_cols, 2, arma::fill::zeros);
  lls = arma::fmat(poi_matrix.data.n_cols, 5, arma::fill::zeros);
  W2 = arma::fmat(poi_matrix.data.n_rows, poi_matrix.data.n_cols,
                  arma::fill::ones);

  cov_no_int_names.resize(num_parms);
  cov_int_names.resize(num_parms2);
  // set row names
  for (auto &col_name : covar_matrix.col_names)
  {
    row_names[col_name.first] = col_name.second;
    row_names2[col_name.first] = col_name.second;
    cov_no_int_names.at(col_name.second) = col_name.first;
    cov_int_names.at(col_name.second) = col_name.first;
  }

  // set interaction names
  for (auto &col_name : no_interactions.col_names)
  {
    row_names["poi"] = covar_matrix.col_names.size();
    cov_no_int_names.at(covar_matrix.col_names.size()) = "poi";
  }

  for (auto &col_name : interactions.col_names)
  {
    row_names2[col_name.first] =
        covar_matrix.col_names.size() + col_name.second;
    cov_int_names.at(col_name.second + covar_matrix.col_names.size()) =
        col_name.first;
  }

  // set col names
  col_names = covar_matrix.row_names;
  col_names2 = covar_matrix.row_names;

  for (arma::uword v = 0; v < poi_matrix.data.n_cols; v++)
  {
    arma::uvec G_na = arma::find_nonfinite(poi_matrix.data.col(v));
    for (arma::uword i = 0; i < G_na.n_elem; i++)
    {
      W2(G_na(i), v) = 0.0;
      poi_matrix.data(G_na(i), v) = 0.0;
    }
  }

  poi_sqrd_mat = arma::square(poi_matrix.data);
  srt_cols = poi_matrix.sort_map(false);
  srt_rows = covar_matrix.sort_map(true);
  this->interactions = &interactions;
  this->no_interactions = &no_interactions;
}

void FRResult::set_lls(double ll1, double ll2, double lrs, double lrs_pval,
                       int num_g, int idx)
{
  lls.at(idx, 0) = ll1;
  lls.at(idx, 1) = ll2;
  lls.at(idx, 2) = lrs;
  lls.at(idx, 3) = lrs_pval;
  lls.at(idx, 4) = num_g;
}

void FRResult::set_betas_fit1(arma::fcolvec &beta, arma::fcolvec &se,
                              arma::fcolvec &pval, int idx)
{
  beta_est.col(idx) = beta;
  se_beta.col(idx) = se;
  neglog10_pvl.col(idx) = pval;
}

void FRResult::set_betas_fit2(arma::fcolvec &beta, arma::fcolvec &se,
                              arma::fcolvec &pval, int idx)
{
  beta_est2.col(idx) = beta;
  se_beta2.col(idx) = se;
  neglog10_pvl2.col(idx) = pval;
}

void FRResult::write_to_file(std::string dir, std::string file_name,
                             int stratum, int process_id)
{
  int n_parms = beta_est.n_rows;
  int n_parms2 = beta_est2.n_rows;

  std::vector<std::string> sorted_row_names = srt_cols;
  std::vector<std::string> sorted_row_names2 =
      srt_cols; // TODO: double check this size and names

  // create dir if it doesn't exist
  fs::create_directory(dir);
  // Rcpp::Rcout << "Dir created or already exists" << std::endl;
  if (srt_cols.size() != beta_est.n_cols)
  {
    Rcpp::Rcout << "Error: The size of poi_names does not match the number of "
                   "columns in the beta matrix."
                << std::endl;
    // return;
  }

  std::stringstream ss;
  ss << dir << "/" << file_name << "_stratum_" << stratum + 1 << "_"
     << process_id << ".tsv";
  std::string result_file = ss.str();
  std::ofstream outfile;

  auto tmp_idx = std::find(cov_int_names.begin(), cov_int_names.end(), "poi");

  if (tmp_idx == cov_int_names.end())
  {
    Rcpp::stop("Couldn't find a poi column");
  }
  int idx = tmp_idx - cov_int_names.begin();
  Rcpp::Rcout << "POIs start from idx: " << idx << std::endl;
  int num_poi = cov_int_names.size() - idx;
  // TODO: Change this to use poi_sqrd
  int num_poi_sqrd = cov_int_names.size() - idx;

  if (fs::exists(result_file))
  {
    outfile.open(result_file, std::ios::app);
    if (!outfile.is_open())
    {
      Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
    }
    // Rcpp::Rcout << "File already exists and opened for writing/appending." <<
    // std::endl;
  }
  else
  {
    outfile.open(result_file);
    if (!outfile.is_open())
    {
      Rcpp::stop("Error: Unable to open file for writing: %s", result_file);
    }

    outfile
        << "POI\tN\tDF_fit1\tEstimate_fit1_poi\tSE_fit1_poi\tmlog10P_fit1_poi\tAbs_"
           "Err_fit1\tRel_Err_fit1\titers_fit1\tDF_fit2";

    // estimate, se, pval for each poi and poi interactions
    for (int i = idx; i < cov_int_names.size(); i++)
    {
      outfile << "\tEstimate_" << cov_int_names[i]
              << "\tSE_" << cov_int_names[i] << "\tmlog10P_" << cov_int_names[i];
    }

    // estimate, se, pval for each poi^2 and poi^2 interactions
    for (int i = idx; i < cov_int_names.size(); i++)
    {
      outfile << "\tEstimate_" << cov_int_names[i] << "_sqrd"
              << "\tSE_" << cov_int_names[i] << "_sqrd"
              << "\tmlog10P_" << cov_int_names[i] << "_sqrd";
    }

    outfile << "\tAbs_Err_fit2\tRel_Err_fit2\titers_fit2\tLL_fit1\tLL_fit2\tLRS\tLRS_mlog10pvl\tnumG" << std::endl;
  }

  outfile << std::fixed << std::setprecision(17);

  std::stringstream buffer;
  for (int col = 0; col < (int)beta_est.n_cols; col++)
  {
    std::string poi_name = srt_cols[col];
    float abs_err_val = beta_abs_errs.at(col);
    float rel_err_val = beta_rel_errs.at(col);
    float iter1 = iters.at(col, 0);
    float iter2 = iters.at(col, 1);
    float abs_err_val2 = beta_abs_errs2.at(col);
    float rel_err_val2 = beta_rel_errs2.at(col);
    float ll1 = lls.at(col, 0);
    float ll2 = lls.at(col, 1);
    float lrs = lls.at(col, 2);
    float lrs_pval = lls.at(col, 3);
    float num_G = lls.at(col, 4);
    if (num_G == 0)
    {
      continue;
    }
    int N = arma::as_scalar(arma::sum(W2.col(col), 0));
    int df = N - n_parms;
    int df2 = N - n_parms2;
    int adder2 = 1;
    buffer << poi_name << "\t" << N << "\t" << df << "\t" << beta_est.at(idx, col) << "\t" << se_beta.at(idx, col) << "\t" << neglog10_pvl.at(idx, col)
           << "\t" << abs_err_val << "\t" << rel_err_val << "\t" << iter1 << "\t" << df2;

    int row2 = idx;
    for (; row2 < (int)beta_est2.n_rows; row2 += adder2)
    {
      buffer << "\t" << beta_est2.at(row2, col) << "\t" << se_beta2.at(row2, col) << "\t"
             << neglog10_pvl2.at(row2, col);
    }
    // TODO: Change this to use beta_sqrd values
    int row3 = idx;
    for (; row3 < (int)beta_est2.n_rows; row3 += adder2)
    {
      buffer << "\t" << beta_est2.at(row3, col) << "\t" << se_beta2.at(row3, col) << "\t"
             << neglog10_pvl2.at(row3, col);
    }

    // outfile << "\tAbs_Err_fit2\tRel_Err_fit2\titers_fit2\tLL_fit1\tLL_fit2\tLRS\tLRS_mlog10pvl\tnumG" << std::endl;
    buffer << "\t" << abs_err_val2 << "\t" << rel_err_val2
           << "\t" << iter2 << "\t" << ll1 << "\t" << ll2 << "\t" << lrs
           << "\t" << lrs_pval << "\t" << num_G << std::endl;
  }
  outfile << buffer.str();
  outfile.close();
}

void FRResult::concatenate(std::string output_dir,
                           std::string file_name_prefix,
                           std::string file_concatenation_prefix)
{
  std::set<int> stratums;
  std::map<int, std::set<std::string>> stratum_files;
  for (const auto &entry : fs::directory_iterator(output_dir))
  {
    if (entry.path().extension() == ".tsv" &&
        entry.path().filename().string().find(file_name_prefix) !=
            std::string::npos)
    {
      // Extract stratum from filename
      std::string filename = entry.path().filename().stem().string();
      std::vector<std::string> tokens;
      std::stringstream ss(filename);
      std::string token;
      while (std::getline(ss, token, '_'))
      {
        tokens.push_back(token);
      }
      if (tokens.size() >= 2)
      {
        try
        {
          int stratum =
              std::stoi(tokens[tokens.size() - 2]); // Second-to-last token
          stratums.insert(stratum);
          stratum_files[stratum].insert(filename);
        }
        catch (const std::invalid_argument &ex)
        {
          Rcpp::Rcerr << "Invalid stratum found in filename: " << filename
                      << std::endl;
        }
      }
    }
  }

  // Concatenate files for each unique stratum
  for (int stratum : stratums)
  {
    std::string outputFile = output_dir + "/" + file_concatenation_prefix +
                             "_" + file_name_prefix + "_stratum_" +
                             std::to_string(stratum) + ".tsv";

    // Check if the output file already exists
    if (fs::exists(outputFile))
    {
      Rcpp::Rcerr << "Output file already exists: " << outputFile
                  << ". It will be overwritten.\n";
    }

    std::ofstream out(outputFile);

    if (!out.is_open())
    {
      Rcpp::Rcerr << "Failed to open the output file: " << outputFile
                  << std::endl;
      continue;
    }
    bool header_written = false;
    for (const auto &filename : stratum_files[stratum])
    {
      fs::path file_path = output_dir + "/" + filename + ".tsv";
      std::ifstream in(file_path.string());

      if (!in.is_open())
      {
        Rcpp::Rcerr << "Failed to open " << file_path << std::endl;
        continue;
      }
      // Skip header line if it has already been written
      if (header_written)
      {
        std::string headerLine;
        std::getline(in, headerLine);
      }
      else
      {
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