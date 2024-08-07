// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <covariate.h>
#include <fr_matrix.h>
#include <string>
#include <vector>

std::vector<std::string> Covariate::split(std::string &val, std::string delim) {
  std::vector<std::string> split_result;

  size_t pos = 0;
  std::string token;

  while ((pos = val.find(delim)) != std::string::npos) {
    token = val.substr(0, pos);
    split_result.push_back(token);
    val.erase(0, pos + delim.length());
  }

  split_result.push_back(val);

  return split_result;
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
  if (levels.empty() && type != "numeric") {
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
  candidate_cols.data = arma::fmat(df.row_names.size(), 1, arma::fill::zeros);
  candidate_cols.row_names = df.row_names;

  if (type == "numeric") {
    // Rcpp::Rcout << "numeric data" << std::endl;
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

    int count = 0;
    std::vector<std::string> col_vals = df.get_col_str(name);
    for (const std::string &lev : levels) {
      arma::fmat res(df.data.n_rows, 1);
      for (arma::uword i = 0; i < df.data.n_rows; i++) {
        res(i, 0) = (col_vals[i] == lev) ? 1 : 0;
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

      // Rcpp::Rcout << "standardized " << name << " " <<
      // candidate_cols.col_names[name] << std::endl;
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
      // Rcpp::Rcout << "rsquared: " << rsquared << std::endl;
      // Rcpp::Rcout << "sse: " << sse << std::endl;
      // Rcpp::Rcout << "ssa: " << ssa << std::endl;
      // zty.brief_print();
      Rcpp::Rcout
          << "Candidate column " << can_col_name
          << " was not added to design matrix due to potential of colinearity"
          << std::endl;
    }
    // Rcpp::Rcout << "join horiz Z covar" << std::endl;
    //  Z = arma::join_horiz(Z, candidate_cols.data.col(can_col_name.second));
  }

  for (auto retained_col : retained_cols) {
    X_mat.data = arma::join_horiz(
        X_mat.data,
        candidate_cols.data.col(candidate_cols.col_names[retained_col]));
  }
  X_mat.col_names.insert(retained_col_map.begin(), retained_col_map.end());
}
