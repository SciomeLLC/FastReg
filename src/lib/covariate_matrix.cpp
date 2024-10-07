#include <covariate_matrix.h>

FRMatrix CovariateMatrix::create_design_matrix() {
  design_matrix.row_names_arr.resize(reader.row_names.size());
  for (size_t i = 0; i < reader.row_names.size(); i++) {
    design_matrix.row_names[reader.row_names[i]] = i;
    design_matrix.row_names_arr.at(i) = reader.row_names[i];
  }

  if (has_intercept) {
    design_matrix.data = arma::fmat(reader.row_names.size(), 1, fill::ones);
    design_matrix.col_names["Intercept"] = 0;
    design_matrix.col_names_arr.push_back("Intercept");
  }

  for (Covariate cov : covariates) {
    cov.create_matrix(reader.values, reader.row_names);
    FRMatrix cov_mat = cov.filter_colinear(design_matrix, colinearity_rsq);
    if (cov_mat.data.n_rows == 0) {
      continue;
    }
    design_matrix.join(cov_mat);
  }
  return design_matrix;
}

int CovariateMatrix::find_col_idx(std::string col_name) {
  for (size_t i = 0; i < headers.size(); i++) {
    if (col_name == headers.at(i)) {
      return i;
    }
  }
  Rcpp::Rcout << "WARNING: column name: " << col_name << " not found in file "
              << file_name << std::endl;
  return -1;
}

void CovariateMatrix::set_cov_col_idx() {
  for (size_t i = 0; i < covariates.size(); i++) {
    std::string cov_name = covariates.at(i).name;
    auto idx = std::find(headers.begin(), headers.end(), cov_name);
    if (idx != headers.end()) {
      int index = std::distance(headers.begin(), idx);
      covariates.at(i).set_col_idx(index);
    } else {
      Rcpp::stop("ERROR: Covariate column name '%s' not found in %s", cov_name,
                 file_name);
    }
  }
}

std::vector<std::string> CovariateMatrix::get_col_names(std::vector<int> idx) {
  std::vector<std::string> col_names(idx.size());
  for (size_t i = 0; i < idx.size(); i++) {
    col_names.push_back(headers.at(idx[i]));
  }
  return col_names;
}
