// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <covariate.h>
#include <fr_matrix.h>
#include <string>
#include <utils.h>
#include <vector>

std::vector<std::string> Covariate::split(std::string &val, char delim) {
  // std::vector<std::string> split_result;
  val.erase(std::remove_if(val.begin(), val.end(),
                           [](char i) { return (i == '\r'); }),
            val.end()); // remove the carriage return if it is there
  std::vector<std::string> tokens(0);
  std::stringstream stream(val);
  std::string temp;

  // Loop over the stringstream until the end of the string
  while (!stream.eof()) {
    std::getline(stream, temp, delim);
    temp.erase(std::remove_if(temp.begin(), temp.end(), ::isspace),
               temp.end()); // remove all whitespace from the value
    tokens.push_back(temp);
  }

  return tokens;
}

void Covariate::set_col_idx(int idx) { col_idx = idx; }

void Covariate::generate_levels(std::vector<std::string> col_vals) {
  std::vector<std::string> unique_vals = col_vals;
  auto it =
      std::unique(unique_vals.begin(), unique_vals.end()); // skip header value
  unique_vals.resize(std::distance(unique_vals.begin(), it));

  for (unsigned int i = 0; i < unique_vals.size(); i++) {
    levels.push_back(unique_vals[i]);
  }
  if (ref_level.empty()) {
    ref_level = unique_vals.at(0);
  }
}

FRMatrix Covariate::create_numeric_matrix(std::vector<std::string> col_vals) {
  FRMatrix candidate_mat;
  candidate_mat.data = arma::fmat(col_vals.size(), 1, arma::fill::zeros);
  
  // fill col names
  candidate_mat.col_names[name] = 0;
  col_names_arr.push_back(name);
  candidate_mat.col_names_arr.push_back(name);

  for (size_t row = 0; row < col_vals.size(); row++) {
    try {
      candidate_mat.data(row, 0) = std::stof(col_vals[row]);
    } catch (const std::invalid_argument &e) {
      if (isWhitespace(col_vals[row]) || col_vals[row].empty()) {
        candidate_mat.data(row, 0) = NAN;
      } else {
        Rcpp::stop("Expected numeric value in column: %s but found invalid "
                   "string '%s'",
                   name, col_vals[row]);
      }
    }
  }
  
  return candidate_mat;
}

FRMatrix
Covariate::create_categorical_matrix(std::vector<std::string> col_vals) {
  FRMatrix candidate_mat;
  // move reference level to the start of the levels list
  levels.erase(std::find(levels.begin(), levels.end(), ref_level));
  levels.insert(levels.begin(), ref_level);
  int count = 0;

  for (const std::string &lev : levels) {
    arma::fmat res(col_vals.size(), 1);
    for (size_t i = 0; i < col_vals.size(); i++) {
      if (col_vals[i].empty() || isWhitespace(col_vals[i])) {
        res(i, 0) = NAN;
      } else {
        res(i, 0) = (col_vals[i] == lev) ? 1 : 0;
      }
    }
    // generate column names for each level
    std::string col_name;
    if (count == 0) {
      col_name = name + " (" + lev + ")";
    } else {
      col_name = name + " (" + lev + " vs. " + ref_level + ")";
    }

    candidate_mat.data.insert_cols(count, res);
    candidate_mat.col_names_arr.push_back(col_name);
    candidate_mat.col_names[col_name] = count;
    col_names_arr.push_back(col_name);
    count++;
  }

  return candidate_mat;
}

void Covariate::standardize_matrix(FRMatrix &frmat) {
  if (standardize) {
    for (arma::uword i = 0; i < frmat.data.n_cols; ++i) {
      arma::fvec col = frmat.data.col(i);
      col = (col - arma::mean(col)) / arma::stddev(col);
      frmat.data.col(i) = col;
    }
  }
}

void Covariate::create_matrix(std::vector<std::vector<std::string>> values,
                              std::vector<std::string> row_names) {
  // Get covariate column values
  std::vector<std::string> col_vals;
  for (std::vector<std::string> line : values) {
    col_vals.push_back(line.at(col_idx));
  }
  col_vals.erase(col_vals.begin()); // remove th header value

  // Create the matrix for the covariate
  if (cov_type == "numeric") {
    frmat = create_numeric_matrix(col_vals);
  } else {
    // Generate unique covariate levels if categorical and user-specified levels
    // are empty
    if (levels.empty()) {
      generate_levels(col_vals);
    }
    frmat = create_categorical_matrix(col_vals);
  }
  
  // Standardize if required
  standardize_matrix(frmat);

  // set row names
  frmat.row_names_arr.resize(row_names.size());
  for (size_t i = 0; i < row_names.size(); i++) {
    frmat.row_names[row_names[i]] = i;
    frmat.row_names_arr.at(i) = row_names[i];
  }
}

FRMatrix Covariate::filter_colinear(FRMatrix &design_mat,
                                    double colinearity_rsq) {
  int nS = frmat.data.n_rows;
  arma::fmat Z;
  arma::fmat full_Z = design_mat.data;
  std::vector<int> cols_to_keep;

  // Test each column for colinearity as they are added to full_Z matrix
  for (arma::uword i = 0; i < frmat.data.n_cols; i++) {
    std::string colname = col_names_arr.at(i);
    arma::fcolvec col_vec = arma::fcolvec(frmat.data.col(i));
    Z = full_Z;

    if (i == 0 && full_Z.n_cols == 0) {
      full_Z.insert_cols(full_Z.n_cols, col_vec);
      cols_to_keep.push_back(i);
      continue;
    }

    // Skip reference level if categorical
    if (cov_type != "numeric" && i == 0) {
      continue;
    }

    arma::fmat y = arma::reshape(col_vec, nS, 1);
    // Add current column to Z matrix to test for colinearity
    Z = arma::join_horiz(Z, col_vec);

    arma::uvec valid_rows = arma::uvec(nS, arma::fill::zeros);
    arma::uvec w = arma::uvec(nS, arma::fill::ones);
    for (int j = 0; j < nS; j++) {
      if (y.row(j).has_nan() || Z.row(j).has_nan()) {
        w[j] = 0;
      }
    }

    Z = Z.rows(arma::find(w > 0));
    y = y.rows(arma::find(w > 0));

    float ssa = arma::accu(arma::square(y));
    arma::fmat zty = Z.t() * y;
    arma::fmat g = arma::pinv(Z.t() * Z);
    float sse = arma::accu(arma::square(y - Z * g * zty));
    float rsquared = 1 - sse / ssa;

    if (rsquared <= colinearity_rsq) {
      full_Z.insert_cols(full_Z.n_cols, col_vec);
      cols_to_keep.push_back(i);
    } else {
      Rcpp::Rcout << "SSE: " << sse << " rsquared: " << rsquared
                  << " ssa: " << ssa << std::endl;
      Rcpp::Rcout
          << "Covariate column " << colname
          << " was not added to design matrix due to potential of colinearity."
          << std::endl;
    }
  }

  // Create the resulting matrix
  int col_count = 0;
  FRMatrix res_mat;
  for (int idx : cols_to_keep) {
    res_mat.data.insert_cols(col_count, frmat.data.col(idx));
    res_mat.col_names[col_names_arr.at(idx)] = col_count;
    res_mat.col_name_str_arr.push_back(col_names_arr.at(idx));
    col_count++;
  }
  res_mat.row_names = frmat.row_names;
  res_mat.row_names_arr = frmat.row_names_arr;
  return res_mat;
}

void Covariate::add_to_matrix(FRMatrix &df, FRMatrix &X_mat,
                              double colinearity_rsq) {
  int cov_idx = df.col_names[name];
  int nS = df.data.n_rows;
  int num_cols = X_mat.data.n_cols;
  std::vector<std::string> temp_col_names;
  if (cov_idx == -1) {
    Rcpp::stop("Covariate %s not found in covariate data file", name);
  }

  // Get levels and reference level if not specified
  if (levels.empty() && cov_type != "numeric") {
    std::vector<std::string> col_vals = df.get_col_str(name);
    auto it = std::unique(col_vals.begin(), col_vals.end());
    col_vals.resize(std::distance(col_vals.begin(), it));
    for (unsigned int i = 0; i < col_vals.size(); i++) {
      levels.push_back(col_vals[i]);
    }
    if (ref_level.empty()) {
      ref_level = col_vals.at(0);
    }
  }

  FRMatrix candidate_cols;
  candidate_cols.row_names = df.row_names;

  if (cov_type == "numeric") {
    candidate_cols.data = arma::fmat(df.row_names.size(), 1, arma::fill::zeros);
    candidate_cols.data.insert_cols(0, df.data.col(cov_idx));
    candidate_cols.col_names[name] = 0;
    temp_col_names.push_back(name);
  } else {
    if (std::find(levels.begin(), levels.end(), ref_level) == levels.end()) {
      Rcpp::stop("covariate.ref.level must be one of unique value of "
                 "covariate.levels.");
    }
    levels.erase(std::find(levels.begin(), levels.end(), ref_level));
    if (X_mat.data.n_cols == 0) {
      levels.insert(levels.begin(), ref_level);
    }
    // size_t num_levels = levels.size();

    int count = 0;
    std::vector<std::string> col_vals = df.get_col_str(name);
    for (const std::string &lev : levels) {
      arma::fmat res(df.data.n_rows, 1);
      for (arma::uword i = 0; i < df.data.n_rows; i++) {
        if (col_vals[i].empty() || isWhitespace(col_vals[i])) {
          res(i, 0) = NAN;
        } else {
          res(i, 0) = (col_vals[i] == lev) ? 1 : 0;
        }
      }
      candidate_cols.data.insert_cols(count, res);
      std::string col_name;
      // std::string lev_str = (lev.empty() ? "" : lev);
      if (X_mat.data.n_cols == 0) {
        col_name = name + " (" + lev + ")";
      } else {
        col_name = name + " (" + lev + " vs. " + ref_level + ")";
      }

      candidate_cols.col_names[col_name] = count;
      temp_col_names.push_back(col_name);
      count++;
    }
    // candidate_cols.data.print();
  }
  if (standardize) {
    for (arma::uword i = 0; i < candidate_cols.data.n_cols; ++i) {
      arma::fvec col = candidate_cols.data.col(i);
      col = (col - arma::mean(col)) / arma::stddev(col);
      candidate_cols.data.col(i) = col;
    }
  }

  // int nl = candidate_cols.data.n_cols;
  int temp_col_count = 0;
  std::unordered_map<std::string, int> retained_col_map;
  std::vector<std::string> retained_cols;
  for (auto &can_col_name : temp_col_names) {

    if (num_cols == 0 && candidate_cols.col_names[can_col_name] == 0) {
      retained_col_map[can_col_name] =
          candidate_cols.col_names[can_col_name] + num_cols;
      retained_cols.push_back(can_col_name);
      temp_col_count++;
      // retain[col_idx] = true;
      continue;
    }
    arma::fmat y = arma::reshape(
        candidate_cols.data.col(candidate_cols.col_names[can_col_name]), nS, 1);

    arma::fmat Z = X_mat.data;
    for (auto retained_col : retained_col_map) {
      Z = arma::join_horiz(Z,
                           arma::fmat(candidate_cols.data.col(
                               candidate_cols.col_names[retained_col.first])));
    }
    arma::uvec w = arma::uvec(nS, arma::fill::ones);
    for (int i = 0; i < nS; i++) {
      if (y.row(i).has_nan()) {
        w[i] = 0;
        continue;
      }
      if (Z.row(i).has_nan()) {
        w[i] = 0;
        continue;
      }
    }

    Z = Z.rows(arma::find(w > 0));
    y = arma::fmat(y.rows(arma::find(w > 0)));

    float ssa = arma::accu(arma::square(y));
    arma::fmat zty = Z.t() * y;
    arma::fmat g = arma::pinv(Z.t() * Z);

    float sse = arma::accu(arma::square(y - Z * g * zty));
    float rsquared = 1 - sse / ssa;

    if (rsquared <= colinearity_rsq) {
      retained_col_map[can_col_name] = temp_col_count + num_cols;
      retained_cols.push_back(can_col_name);
      temp_col_count++;
    } else {
      Rcpp::Rcout
          << "Candidate column " << can_col_name
          << " was not added to design matrix due to potential of colinearity"
          << std::endl;
    }
  }

  for (auto retained_col : retained_cols) {
    X_mat.data = arma::join_horiz(
        X_mat.data,
        candidate_cols.data.col(candidate_cols.col_names[retained_col]));
  }
  X_mat.col_names.insert(retained_col_map.begin(), retained_col_map.end());
}

void Covariate::print() {
  Rcpp::Rcout << "covariate: " << name << " type: " << cov_type;
  if (levels.size() > 0 && !levels[0].empty()) {
    Rcpp::Rcout << " levels: ";
    for (size_t i = 0; i < levels.size(); i++) {
      Rcpp::Rcout << levels.at(i) << ",";
    }
    Rcpp::Rcout << " ref_level: " << ref_level;
  }
  Rcpp::Rcout << ", standardize: " << standardize << std::endl;
}