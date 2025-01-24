#include <pheno_matrix.h>
#include <utils.h>
// Creates the phenotype matrix (FRMatrix) and assigns it to the class member
// "mat"
FRMatrix PhenoMatrix::create_matrix() {
  // Check for id column and get the index
  auto idx = std::find(headers.begin(), headers.end(), id_col);
  if (idx == headers.end()) {
    Rcpp::stop("Could not find column name: %s in %s", id_col, file_name);
  }
  id_col_idx = idx - headers.begin();

  // Initialize the matrix
  mat.data = arma::fmat(reader.row_names.size(), phenotypes.size(), arma::fill::zeros);
  // Check for phenotype column and get the index
  for (auto it = phenotypes.begin(); it != phenotypes.end(); ++it){
    int index = std::distance(phenotypes.begin(), it);
    idx = std::find(headers.begin(), headers.end(), phenotypes[index]);
    if (idx == headers.end()) {
      Rcpp::stop("Could not find column name: %s in %s", phenotype, file_name);
    }
    int pheno_idx = idx - headers.begin();

    // Set column name and index
    std::string phen = phenotypes[index];
    mat.col_names[phen] = index;
    mat.col_name_str_arr.push_back(phen);
    mat.col_names_arr.push_back(phen);
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
        mat.data(i - 1, index) = std::stof(row[pheno_idx]);
      } catch (const std::invalid_argument &e) {
        if (isWhitespace(row[pheno_idx]) || row[pheno_idx].empty()) {
          mat.data(i - 1, index) = NAN;
        } else {
          Rcpp::stop("Expected numeric value in column: %s but found '%s'",
                    phen, row[pheno_idx]);
        }
      }
      i++;
    }
  }
  // set the row names
  mat.row_names_arr = reader.row_names;
  for (size_t j = 0; j < reader.row_names.size(); j++) {
    mat.row_names[reader.row_names[j]] = j;
  }

  return mat;
}