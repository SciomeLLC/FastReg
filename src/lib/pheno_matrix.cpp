#include <pheno_matrix.h>

// Creates the phenotype matrix (FRMatrix) and assigns it to the class member
// "mat"
FRMatrix PhenoMatrix::create_matrix() {
  // Check for id column and get the index
  auto idx = std::find(headers.begin(), headers.end(), id_col);
  if (idx == headers.end()) {
    Rcpp::stop("Could not find column name: %s in %s", id_col, file_name);
  }
  id_col_idx = idx - headers.begin();

  // Check for phenotype column and get the index
  idx = std::find(headers.begin(), headers.end(), phenotype);
  if (idx == headers.end()) {
    Rcpp::stop("Could not find column name: %s in %s", phenotype, file_name);
  }
  int pheno_idx = idx - headers.begin();

  // Set column name and index
  mat.col_names[phenotype] = 0;
  mat.col_name_str_arr.push_back(phenotype);
  mat.col_names_arr.push_back(phenotype);

  // Initialize the matrix
  mat.data = arma::fmat(reader.row_names.size(), 1, arma::fill::zeros);

  // Set row names and values
  int i = 0;
  for (auto &row : reader.values) {
    // skip the header
    if (i == 0) {
      i++;
      continue;
    }
    // Check for NANs and invalid values
    try {
      mat.data(i - 1, 0) = std::stof(row[pheno_idx]);
    } catch (const std::invalid_argument &e) {
      if (isWhitespace(row[pheno_idx]) || row[pheno_idx].empty()) {
        mat.data(i - 1, 0) = NAN;
      } else {
        Rcpp::stop("Expected numeric value in column: %s but found '%s'",
                   phenotype, row[pheno_idx]);
      }
    }
    i++;
  }
  // set the row names
  mat.row_names_arr = reader.row_names;
  for (int j = 0; j < reader.row_names.size(); j++) {
    mat.row_names[reader.row_names[j]] = j;
  }

  return mat;
}

FRMatrix PhenoMatrix::get_squared() {
  FRMatrix sqrd = mat;
  sqrd.col_names[phenotype + "^2"] = 1;
  sqrd.col_names_arr.push_back(phenotype + "^2");
  sqrd.col_name_str_arr.push_back(phenotype + "^2");
  arma::fcolvec sqr = arma::square(mat.data.col(0));
  sqrd.data.insert_cols(1, sqr);

  return sqrd;
}