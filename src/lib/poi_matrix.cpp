#include <poi_matrix.h>

FRMatrix POIMatrix::get_chunk(const std::vector<std::string> &rows,
                              const std::vector<std::string> &cols) {
  return reader->read_chunk(rows, cols);
}

void POIMatrix::check_dominant_recessive(FRMatrix &poi_matrix) {
  if (effect_type == "dominant") {
    poi_matrix.data.elem(find(poi_matrix.data == 2)).fill(1);
  } else if (effect_type == "recessive") {
    poi_matrix.data.elem(find(poi_matrix.data == 1)).fill(0);
    poi_matrix.data.elem(find(poi_matrix.data == 2)).fill(1);
  }
}

void POIMatrix::filter_rows(FRMatrix &chunk, const std::vector<std::string> r_names) {
  std::vector<std::string> drop_rows = set_diff(individuals, r_names);
  int num_dropped = individuals.size() - r_names.size();

  arma::uvec drop_row_idx(drop_rows.size());
  for (size_t i = 0; i < drop_rows.size(); i++) {
    drop_row_idx[i] = chunk.row_names[drop_rows[i]];
  }

  std::unordered_map<std::string, int> new_row_names(r_names.size());
  for (auto &ind : r_names) {
    new_row_names[ind] = chunk.row_names[ind] - num_dropped;
  }
  chunk.data.shed_rows(drop_row_idx);
  chunk.row_names = new_row_names;
}

FRMatrix POIMatrix::filter_genotype(FRMatrix &poi_matrix) {
  FRMatrix filtered;
  if (poi_type != "genotype") {
    return filtered;
  }

  check_threshold(filtered, poi_matrix);
  arma::uvec filtered_col = arma::find(filtered.data.row(5) == 0);

  if (filtered.data.n_cols == 0 ||
      filtered_col.n_elem == poi_matrix.data.n_cols) {
    Rcpp::stop("no POI passed filtering");
  }

  std::vector<std::string> poi_col_names = filtered.sort_map(false);
  int cols_erased = 0;

  for (unsigned int i = 0; i < poi_col_names.size(); i++) {
    if ((unsigned)cols_erased < filtered_col.n_elem &&
        filtered_col[cols_erased] == i) {
      poi_matrix.col_names.erase(poi_col_names[i]);
      cols_erased++;
    } else {
      poi_matrix.col_names[poi_col_names[i]] =
          poi_matrix.col_names[poi_col_names[i]] - cols_erased;
    }
  }
  poi_matrix.data.shed_cols(filtered_col);
  check_dominant_recessive(poi_matrix);
  return filtered;
}

void POIMatrix::check_threshold(FRMatrix &res, FRMatrix &poi_matrix) {

  res.col_names = poi_matrix.col_names;
  res.row_names = {{"Freq (a)", 0},  {"Freq (b)", 1},   {"MAF", 2},
                   {"HWE Chisq", 3}, {"HWE Pvalue", 4}, {"keep", 5}};
  res.data = arma::fmat(res.row_names.size(), poi_matrix.data.n_cols,
                        arma::fill::ones);

  umat isnan_mat =
      umat(poi_matrix.data.n_rows, poi_matrix.data.n_cols, arma::fill::zeros);
  for (uword col = 0; col < poi_matrix.data.n_cols; ++col) {
    uvec nonfinite_indices = arma::find_nonfinite(poi_matrix.data.col(col));
    // isnan_mat.col(col).
    for (uword i = 0; i < nonfinite_indices.n_elem; ++i) {
      isnan_mat(nonfinite_indices(i), col) = 1;
    }
  }

  FRMatrix G_copy = poi_matrix;
  G_copy.data.replace(datum::nan, 0);

  frowvec nS =
      conv_to<frowvec>::from(poi_matrix.data.n_rows - arma::sum(isnan_mat, 0));
  // Rcpp::Rcout << nS << " G n_cols: " << G_copy.data.n_cols << std::endl;
  frowvec nS_2 = 2 * nS;
  frowvec a_freq = arma::sum(G_copy.data, 0) / (nS_2);
  frowvec b_freq = 1 - a_freq;

  frowvec maf_freq = min(a_freq, b_freq);
  // Rcpp::Rcout << " maf_freq: " << maf_freq << std::endl;
  frowvec aa_of = conv_to<frowvec>::from(arma::sum(poi_matrix.data == 0, 0));
  frowvec ab_of = conv_to<frowvec>::from(arma::sum(G_copy.data == 1, 0));
  frowvec bb_of = conv_to<frowvec>::from(arma::sum(G_copy.data == 2, 0));

  frowvec p = (2 * aa_of + 1 * ab_of) / (nS_2);
  frowvec aa_ef = square(p) % nS;
  frowvec ab_ef = 2 * p % (1 - p) % nS;
  frowvec bb_ef = square(1 - p) % nS;

  frowvec HWE_chisq = (square(aa_of - aa_ef) / aa_ef) +
                      (square(ab_of - ab_ef) / ab_ef) +
                      (square(bb_of - bb_ef) / bb_ef);
  arma::rowvec wtf =
      1 - (Rcpp::pchisq(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(HWE_chisq)), 1,
                        true, false));
  frowvec HWE_pval = arma::conv_to<frowvec>::from(wtf);
  frowvec keep = conv_to<frowvec>::from(
      conv_to<frowvec>::from((maf_freq >= maf_threshold)) &&
      conv_to<frowvec>::from((HWE_pval >= hwe_threshold)));

  res.data.row(0) = a_freq;
  res.data.row(1) = b_freq;
  res.data.row(2) = maf_freq;
  res.data.row(3) = HWE_chisq;
  res.data.row(4) = HWE_pval;
  res.data.row(5) = keep;
}