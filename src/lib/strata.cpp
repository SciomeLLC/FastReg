
#include <strata.h>

void Strata::stratify(const std::vector<std::string> &split_by,
                      FRMatrix &covar_df,
                      const std::vector<std::string> &common_individuals) {
  int nsplit_vars = split_by.size();
  if (!split_by[0].empty()) {
    Rcpp::Rcout << "Stratifying - split_by value exists" << std::endl;
    arma::uvec common_individuals_idx(common_individuals.size());
    int ct = 0;
    for (const auto &ind : common_individuals) {
      if (ind.empty()) {
        Rcpp::Rcout << "Stratify found empty id string: " << ind << std::endl;
        continue;
      }
      int row_idx = covar_df.get_row_idx(ind);
      if (row_idx == -1) {
        continue;
      }
      common_individuals_idx[ct] = row_idx;
      ct++;
    }

    ct = 0;
    arma::uvec split_by_idx(split_by.size());

    for (const auto &col : split_by) {
      if (col.empty()) {
        Rcpp::Rcout << "Stratify found empty col string: " << col << std::endl;
        continue;
      }
      int col_idx = covar_df.get_col_idx(col);
      if (col_idx == -1) {
        continue;
      }
      split_by_idx[ct] = col_idx;
      ct++;
    }

    arma::fmat unique_combinations = arma::unique(
        covar_df.data.submat(common_individuals_idx, split_by_idx));
    nstrata = unique_combinations.n_rows;
    ids.resize(nstrata);

    for (int stratum = 0; stratum < nstrata; ++stratum) {
      std::string stratum_key = "";
      std::vector<std::string> stratum_indices;
      for (int col = 0; col < nsplit_vars; ++col) {
        stratum_key += "_" + split_by[col] + "=" +
                       std::to_string((int)unique_combinations(stratum, col));
      }

      for (const auto &row_name : common_individuals) {
        int row_idx = covar_df.get_row_idx(row_name);
        bool match = true;

        for (int col = 0; col < nsplit_vars && match; ++col) {
          int col_idx = covar_df.get_col_idx(split_by[col]);
          if (col_idx >= 0 && row_idx >= 0 &&
              covar_df.data(row_idx, col_idx) !=
                  unique_combinations(stratum, col)) {
            match = false;
          }
        }
        if (match) {
          stratum_indices.push_back(row_name);
        }
      }

      ids[stratum] = stratum_key;
      index_list[stratum_key] = stratum_indices;
    }
  } else {
    nstrata = 1;
    ids.push_back("_");
    index_list["_"] = common_individuals;
  }
}
