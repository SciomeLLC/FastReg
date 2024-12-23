#include <fr_matrix.h>

int FRMatrix::get_row_idx(const std::string &row_name) {
  bool exists = row_names.find(row_name) != row_names.end();
  if (!exists) {
    return -1;
  }
  return row_names[row_name];
}

int FRMatrix::get_col_idx(const std::string &col_name) {
  bool exists = col_names.find(col_name) != col_names.end();
  if (!exists) {
    return -1;
  }
  return col_names[col_name];
}

bool FRMatrix::validate_cols(std::string &names) {
  char delim = ';';

  std::vector<std::string> id_cols = split(names, delim);
  for (std::string &col : id_cols) {
    if (col_names.find(col) == col_names.end()) {
      return false;
    }
  }
  return true;
}

void FRMatrix::shed_rows(std::vector<int> &idx,
                         std::vector<std::string> &new_row_names) {
  data.shed_rows(arma::conv_to<arma::uvec>::from(idx));
  row_names = std::unordered_map<std::string, int>(new_row_names.size());
  for (size_t j = 0; j < data.n_rows; j++) {
    row_names[new_row_names[j]] = j;
  }
}

FRMatrix FRMatrix::get_submat_by_cols(const std::vector<int> &row_idx,
                                      const std::vector<std::string> &names) {
  std::vector<int> col_indices;
  FRMatrix subm;

  // Iterate over the names vector and retrieve the corresponding column indices
  for (const std::string &item : names) {
    auto itr = col_names.find(item);
    if (itr != col_names.end()) {
      col_indices.push_back(itr->second);
      subm.col_names[item] = itr->second;
    }
  }

  // Extract the submatrix using the row and column indices
  subm.data = data.submat(arma::conv_to<arma::uvec>::from(row_idx),
                          arma::conv_to<arma::uvec>::from(col_indices));

  // Assign the row names of the submatrix
  for (const int idx : row_idx) {
    auto itr = std::find_if(
        row_names.begin(), row_names.end(),
        [&](const std::pair<std::string, int> &p) { return p.second == idx; });

    if (itr != row_names.end()) {
      subm.row_names[itr->first] = idx;
    }
  }

  return subm;
}

FRMatrix FRMatrix::get_submat_by_cols(
    const std::vector<int> &row_idx,
    const std::unordered_map<std::string, int> &names) {
  std::vector<int> col_indices;
  FRMatrix subm;

  // Iterate over the names vector and retrieve the corresponding column indices
  for (auto &item : names) {
    auto itr = col_names.find(item.first);
    if (itr != col_names.end()) {
      col_indices.push_back(itr->second);
      subm.col_names[item.first] = itr->second;
    }
  }

  // Extract the submatrix using the row and column indices
  subm.data = data.submat(arma::conv_to<arma::uvec>::from(row_idx),
                          arma::conv_to<arma::uvec>::from(col_indices));

  // Assign the row names of the submatrix
  for (const int idx : row_idx) {
    auto itr = std::find_if(
        row_names.begin(), row_names.end(),
        [&](const std::pair<std::string, int> &p) { return p.second == idx; });

    if (itr != row_names.end()) {
      subm.row_names[itr->first] = idx;
    }
  }

  return subm;
}

std::vector<std::string> FRMatrix::split(const std::string &str_tokens,
                                         char delim) {
  std::string str = str_tokens;
  str.erase(std::remove_if(str.begin(), str.end(),
                           [](char i) { return (i == '\r'); }),
            str.end()); // remove the carriage return if it is there
  std::vector<std::string> toks(0);
  std::stringstream stream(str);
  std::string temp;
  // int i = 1;
  // Loop over the stringstream until newline '\n' is hit
  while (!stream.eof()) {
    std::getline(stream, temp, delim);
    temp.erase(std::remove_if(temp.begin(), temp.end(), ::isspace),
               temp.end()); // remove all whitespace from the value
    toks.push_back(temp);
  }

  return toks;
}

std::vector<std::string> FRMatrix::get_col_str(const std::string &col_name) {

  bool exists = col_names_str.find(col_name) != col_names_str.end();
  if (!exists) {
    Rcpp::stop("Column name " + col_name +
               " doesn't match any non-numeric columns.");
  }
  int col_idx = col_names_str[col_name];
  std::vector<std::string> col_vals;
  for (std::vector<std::string> row : str_data) {
    col_vals.push_back(row[col_idx]);
  }

  return col_vals;
}

std::vector<std::string> FRMatrix::sort_map(bool rows) {
  std::unordered_map<std::string, int> temp;
  if (rows) {
    temp = row_names;
  } else {
    temp = col_names;
  }

  std::vector<std::string> sorted_arr(temp.size());
  for (auto item : temp) {
    sorted_arr[item.second] = item.first;
  }
  return sorted_arr;
}

void FRMatrix::write_summary(std::string dir, std::string name, int stratum,
                             int process_id) {
  fs::create_directory(dir);
  std::stringstream ss;
  ss << dir << "/" << name << "_stratum_" << stratum + 1 << "_" << process_id
     << ".tsv";
  std::string file_name = ss.str();

  arma::fmat temp = data.t();
  std::ofstream outfile;

  if (fs::exists(file_name)) {
    outfile = std::ofstream(file_name, std::ios::app);
  } else {
    outfile = std::ofstream(file_name);
    std::vector<std::string> ordered_row_names = sort_map(true);
    for (const auto &row_name : ordered_row_names) {
      outfile << '\t' << row_name;
    }
    outfile << std::endl;
  }

  outfile << std::fixed << std::setprecision(17);
  std::stringstream buffer;
  std::vector<std::string> ordered_col_names = sort_map(false);
  int row_idx = 0;
  buffer << "POI\t";
  for (const auto &col_name : ordered_col_names) {
    buffer << col_name;
    for (size_t col_idx = 0; col_idx < temp.n_cols; ++col_idx) {
      buffer << '\t' << temp(row_idx, col_idx);
    }
    buffer << std::endl;
    row_idx++;
  }
  outfile << buffer.str();
  outfile.close();
}

void FRMatrix::print() {
  std::vector<std::string> sort_cols = sort_map(false);
  std::vector<std::string> sort_rows = sort_map(true);
  for (auto &col : sort_cols) {
    Rcpp::Rcout << col << "\t";
  }
  Rcpp::Rcout << std::endl;

  for (size_t i = 0; i < data.n_rows; i++) {
    Rcpp::Rcout << sort_rows[i] << "\t";
    data.row(i).print();
  }
}

void FRMatrix::write_results(FRMatrix &beta, FRMatrix &se_beta,
                             FRMatrix &neglog10, arma::umat &W2,
                             arma::fcolvec &rel_err, arma::fcolvec &abs_err,
                             arma::fcolvec &iters,
                             std::vector<std::string> poi_names,
                             std::string dir, std::string file_name,
                             int stratum, bool exclude_covars, int process_id) {
  int n_parms = beta.data.n_rows;

  std::vector<std::string> sorted_row_names = beta.sort_map(true);
  // create dir if it doesn't exist
  fs::create_directory(dir);
  // Rcpp::Rcout << "Dir created or already exists" << std::endl;
  if (poi_names.size() != beta.data.n_cols) {
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

  if (fs::exists(result_file)) {
    outfile.open(result_file, std::ios::app);
    if (!outfile.is_open()) {
      Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file
                  << std::endl;
      return;
    }
    // Rcpp::Rcout << "File already exists and opened for writing/appending." <<
    // std::endl;
  } else {
    outfile.open(result_file);
    if (!outfile.is_open()) {
      Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file
                  << std::endl;
      return;
    }
    outfile
        << "POI\tN\tDF\tEffect\tEstimate\tSE\tmlog10P\tAbs_Err\tRel_Err\titers"
        << std::endl;
    // Rcpp::Rcout << "File created for writing." << std::endl;
  }

  outfile << std::fixed << std::setprecision(17);

  std::stringstream buffer;
  for (int col = 0; col < (int)beta.data.n_cols; col++) {
    std::string poi_name = poi_names[col];
    arma::uvec w2_col = W2.col(col);
    double abs_err_val = abs_err.at(col);
    double rel_err_val = rel_err.at(col);
    double iter = iters.at(col);

    int N = arma::as_scalar(arma::sum(w2_col, 0));
    int df = N - n_parms;
    int adder = 1;
    int row = 0;
    if (exclude_covars) {
      adder = n_parms;
      row = n_parms - 1;
    }

    for (; row < (int)beta.data.n_rows; row += adder) {
      std::string effect_name = sorted_row_names[row];
      double estimate = beta.data.at(row, col);
      double std_error = se_beta.data.at(row, col);
      double neglog10_pval = neglog10.data.at(row, col);

      buffer << poi_name << "\t" << N << "\t" << df << "\t" << effect_name
             << "\t" << estimate << "\t" << std_error << "\t" << neglog10_pval
             << "\t" << abs_err_val << "\t" << rel_err_val << "\t" << iter
             << std::endl;
    }
  }
  outfile << buffer.str();
  outfile.close();
}

void FRMatrix::write_vla_results(
    FRMatrix &beta, FRMatrix &se_beta, FRMatrix &neglog10, arma::umat &W2,
    arma::fcolvec &rel_err, arma::fcolvec &abs_err, FRMatrix &beta2,
    FRMatrix &se_beta2, FRMatrix &neglog102, arma::fcolvec &rel_err2,
    arma::fcolvec &abs_err2, arma::fmat &lls, arma::fmat &iters,
    std::vector<std::string> poi_names, std::string dir, std::string file_name,
    int stratum, bool exclude_covars, int process_id) {
  int n_parms = beta.data.n_rows;
  int n_parms2 = beta2.data.n_rows;

  std::vector<std::string> sorted_row_names = beta.sort_map(true);
  // create dir if it doesn't exist
  fs::create_directory(dir);
  // Rcpp::Rcout << "Dir created or already exists" << std::endl;
  if (poi_names.size() != beta.data.n_cols) {
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

  if (fs::exists(result_file)) {
    outfile.open(result_file, std::ios::app);
    if (!outfile.is_open()) {
      Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file
                  << std::endl;
      return;
    }
    // Rcpp::Rcout << "File already exists and opened for writing/appending." <<
    // std::endl;
  } else {
    outfile.open(result_file);
    if (!outfile.is_open()) {
      Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file
                  << std::endl;
      return;
    }
    outfile
        << "POI\tEffect\tN\tDF_fit1\tEstimate_fit1\tSE_fit1\tmlog10P_fit1\tAbs_"
           "Err_fit1\tRel_"
           "Err_fit1\titers_fit1\tDF_fit2\tEstimate_fit2\tSE_fit2\tmlog10P_"
           "fit2\tAbs_"
           "Err_fit2\tRel_"
           "Err_fit2\titers_fit2\tLL_fit1\tLL_fit2\tLRS\tLRS_mlog10pvl\tnumG"
        << std::endl;
    // Rcpp::Rcout << "File created for writing." << std::endl;
  }

  outfile << std::fixed << std::setprecision(17);

  std::stringstream buffer;
  for (int col = 0; col < (int)beta.data.n_cols; col++) {
    std::string poi_name = poi_names[col];
    arma::uvec w2_col = W2.col(col);
    float abs_err_val = abs_err.at(col);
    float rel_err_val = rel_err.at(col);
    float iter1 = iters.at(col, 0);
    float iter2 = iters.at(col, 1);
    float abs_err_val2 = abs_err2.at(col);
    float rel_err_val2 = rel_err2.at(col);
    float ll1 = lls.at(col, 0);
    float ll2 = lls.at(col, 1);
    float lrs = lls.at(col, 2);
    float lrs_pval = lls.at(col, 3);
    float num_G = lls.at(col, 4);

    int N = arma::as_scalar(arma::sum(w2_col, 0));
    int df = N - n_parms;
    int df2 = N - n_parms2;
    int adder = 1;
    int adder2 = 1;
    int row = 0;
    int row2 = 0;
    if (exclude_covars) {
      adder = n_parms;
      row = n_parms - 1;
      adder2 = n_parms2;
      row2 = n_parms2 - 1;
    }

    for (; row < (int)beta.data.n_rows; row += adder) {
      std::string effect_name = sorted_row_names[row];
      float estimate = beta.data.at(row, col);
      float std_error = se_beta.data.at(row, col);
      float neglog10_pval = neglog10.data.at(row, col);
      float estimate2 = beta2.data.at(row2, col);
      float std_error2 = se_beta2.data.at(row2, col);
      float neglog10_pval2 = neglog102.data.at(row2, col);
      // outfile
      //       <<
      //       "POI\tEffect\tN\tDF_fit1\tEstimate_fit1\tSE_fit1\tmlog10P_fit1\tAbs
      //       "
      //          "Err_fit1\tRel
      //          Err_fit1\tDF_fit2\tEstimate_fit2\tSE_fit2\tmlog10P_fit2\tAbs "
      //          "Err_fit2\tRel
      //          Err_fit2\tLL_fit1\tLL_fit2\tLRS\tLRS_mlog10pvl\titers"
      //       << std::endl;
      buffer << poi_name << "\t" << effect_name << "\t" << N << "\t" << df
             << "\t" << estimate << "\t" << std_error << "\t" << neglog10_pval
             << "\t" << abs_err_val << "\t" << rel_err_val << "\t" << iter1
             << "\t" << df2 << "\t" << estimate2 << "\t" << std_error2 << "\t"
             << neglog10_pval2 << "\t" << abs_err_val2 << "\t" << rel_err_val2
             << iter2 << "\t" << ll1 << "\t" << ll2 << "\t" << lrs << "\t"
             << lrs_pval << "\t" << num_G << std::endl;
      row2 += adder2;
    }
  }
  outfile << buffer.str();
  outfile.close();
}

void FRMatrix::write_convergence_results(FRMatrix &beta,
                                         std::vector<std::string> poi_names,
                                         std::string dir, std::string file_name,
                                         arma::fcolvec &rel_err,
                                         arma::fcolvec &abs_err, int stratum,
                                         int process_id) {
  std::vector<std::string> sorted_row_names = beta.sort_map(true);
  // create dir if it doesn't exist
  fs::create_directory(dir);
  // Rcpp::Rcout << "Dir created or already exists" << std::endl;
  if (poi_names.size() != beta.data.n_cols) {
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

  if (fs::exists(result_file)) {
    outfile.open(result_file, std::ios::app);
    if (!outfile.is_open()) {
      Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file
                  << std::endl;
      return;
    }
    // Rcpp::Rcout << "File already exists and opened for writing/appending." <<
    // std::endl;
  } else {
    outfile.open(result_file);
    if (!outfile.is_open()) {
      Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file
                  << std::endl;
      return;
    }
    outfile << "POI\tAbs Err\tRel Err" << std::endl;
    // Rcpp::Rcout << "File created for writing." << std::endl;
  }

  outfile << std::fixed << std::setprecision(17);

  std::stringstream buffer;
  for (int col = 0; col < (int)beta.data.n_cols; col++) {
    std::string poi_name = poi_names[col];

    double abs_err_val = abs_err.at(col);
    double rel_err_val = rel_err.at(col);

    buffer << poi_name << "\t" << abs_err_val << "\t" << rel_err_val
           << std::endl;
  }
  outfile << buffer.str();
  outfile.close();
}

void FRMatrix::zip_results(std::string output_dir) {
  Rcpp::Environment utils_env("package:utils");
  Rcpp::Function zip = utils_env["zip"];
  if (fs::exists(output_dir)) {
    std::string parent_path = fs::path(output_dir).parent_path().string();
    const auto time_now = std::chrono::system_clock::now();
    const auto time_secs = std::chrono::duration_cast<std::chrono::seconds>(
                               time_now.time_since_epoch())
                               .count();
    Rcpp::Rcout << parent_path << std::endl;

    std::string archive_name =
        parent_path + "/results_" + std::to_string(time_secs) + ".zip";
    zip(archive_name, output_dir);
  }
}

void FRMatrix::concatenate_results(std::string output_dir,
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

void FRMatrix::join(FRMatrix &fmat) {

  if (data.size() == 0) {
    data = fmat.data;
    col_names = fmat.col_names;
    row_names = fmat.row_names;
    return;
  }

  if (fmat.data.n_rows != data.n_rows) {
    Rcpp::Rcout << fmat.data.n_rows << " vs " << data.n_rows << std::endl;
    Rcpp::stop("Error: cannot join_horiz matrices with different row counts.");
  }
  // join the matrices
  data = arma::join_horiz(data, fmat.data);
  size_t num_cols = col_names.size();

  // Add column names and increase indices accordingly
  col_names_arr.resize(col_names_arr.size() + fmat.col_names.size());
  for (auto el : fmat.col_names) {
    col_names.emplace(el.first, el.second + num_cols);
    col_names_arr.at(el.second + num_cols) = el.first;
  }
}
