

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <fr_matrix.h>
#include <covariate.h>
#include <omp.h>
#pragma once

using namespace arma;

// void trim(std::string &s)
// {
//     s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
//     s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
// }

void transform_poi(FRMatrix &G, std::string effect_type = "additive") {
    if (effect_type == "dominant") {
        G.data.elem(find(G.data == 2)).fill(1);
    } else if (effect_type == "recessive") {
        G.data.elem(find(G.data == 1)).fill(0);
        G.data.elem(find(G.data == 2)).fill(1);
    }
}

FRMatrix filter_poi(FRMatrix &G, double maf_threshold = 0.01, double hwe_threshold = 0.05) {

    FRMatrix res;
    res.col_names = G.col_names;
    res.row_names = {
        {"Freq (a)", 0}, {"Freq (b)", 1}, {"MAF", 2}, {"HWE Chisq", 3}, {"HWE Pvalue", 4}, {"keep", 5}};
    res.data = arma::mat(res.row_names.size(), G.data.n_cols, arma::fill::ones);
    umat isnan_mat = umat(G.data.n_rows, G.data.n_cols, arma::fill::zeros);
    for (uword col = 0; col < G.data.n_cols; ++col) {
      uvec nonfinite_indices = arma::find_nonfinite(G.data.col(col));
      // isnan_mat.col(col).
      for (uword i = 0; i < nonfinite_indices.n_elem; ++i) {
        isnan_mat(nonfinite_indices(i), col) = 1;
      }
    }
    FRMatrix G_copy = G;
    G_copy.data.replace(datum::nan, 0);
    
    rowvec nS = conv_to<rowvec>::from(G.data.n_rows - arma::sum(isnan_mat, 0));
    rowvec nS_2 = 2 * nS;
    rowvec a_freq = arma::sum(G_copy.data, 0) / (nS_2);
    rowvec b_freq = 1 - a_freq;

    rowvec maf_freq = min(a_freq, b_freq);
    rowvec aa_of = conv_to<rowvec>::from(arma::sum(G.data == 0, 0));
    rowvec ab_of = conv_to<rowvec>::from(arma::sum(G_copy.data == 1, 0));
    rowvec bb_of = conv_to<rowvec>::from(arma::sum(G_copy.data == 2, 0));
    
    rowvec p = (2*aa_of + 1*ab_of) / (nS_2);
    rowvec aa_ef = square(p) % nS;
    rowvec ab_ef = 2 * p % (1 - p) % nS;
    rowvec bb_ef = square(1 - p) % nS;

    rowvec HWE_chisq = (square(aa_of - aa_ef) / aa_ef) + (square(ab_of - ab_ef) / ab_ef) + (square(bb_of - bb_ef) / bb_ef);
    rowvec HWE_pval = 1 - (Rcpp::pchisq(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(HWE_chisq)), 1, true, false));
    rowvec keep = conv_to<rowvec>::from((maf_freq >= maf_threshold) % (HWE_pval >= hwe_threshold));
    res.data.row(0) = a_freq;
    res.data.row(1) = b_freq;
    res.data.row(2) = maf_freq;
    res.data.row(3) = HWE_chisq;
    res.data.row(4) = HWE_pval;
    res.data.row(5) = keep;
    
    return res;
}

// void filter_matrices(FRMatrix& X, FRMatrix& Y, FRMatrix& Z, FRMatrix& G) {
//     std::vector<int> X_index, Y_index, G_index, Z_index;
//     std::vector<std::string> common_subjects, intersection_xy, intersection_xyz;
//
//     std::sort(X.row_names.begin(), X.row_names.end());
//     std::sort(Y.row_names.begin(), Y.row_names.end());
//     std::sort(G.row_names.begin(), G.row_names.end());
//     std::sort(Z.row_names.begin(), Z.row_names.end());
//
//     std::set_intersection(X.row_names.begin(), X.row_names.end(),
//                             Y.row_names.begin(), Y.row_names.end(),
//                             std::back_inserter(intersection_xy));
//
//     std::set_intersection(intersection_xy.begin(), intersection_xy.end(),
//                             G.row_names.begin(), G.row_names.end(),
//                             std::back_inserter(intersection_xyz));
//
//     std::set_intersection(intersection_xyz.begin(), intersection_xyz.end(),
//                             Z.row_names.begin(), Z.row_names.end(),
//                             std::back_inserter(common_subjects));
//
//     for (const std::string& subject : common_subjects) {
//         X_index.push_back(X.row_names[subject]);
//         Y_index.push_back(Y.row_names[subject]);
//         G_index.push_back(G.row_names[subject]);
//         Z_index.push_back(Z.row_names[subject]);
//     }
//
//     X = X.get_submat_by_cols(X_index, X.col_names);
//     Y = Y.get_submat_by_cols(Y_index, Y.col_names);
//     G = G.get_submat_by_cols(G_index, G.col_names);
//     Z = Z.get_submat_by_cols(Z_index, Z.col_names);
// }

void create_Z_matrix(FRMatrix& df, const std::vector<std::string>& poi_covar_interactions, FRMatrix& Z)
{
    std::vector<std::string> interaction_cols;
    if(poi_covar_interactions.size() == 1 && poi_covar_interactions[0].empty()) {
        return;
    }

    for (const auto& interaction : poi_covar_interactions) {
        for (const auto& col_name : df.col_names) {
            if (col_name.first.find(interaction) != std::string::npos) {
                interaction_cols.push_back(col_name.first);
            }
        }
    }

    if(!interaction_cols.empty()) {
        Rcpp::Rcout << "Found poi covar interactions" << std::endl;
        std::vector<int> col_idx;
        int count = 0;
        for (const auto& col_name : interaction_cols) {
            int idx = df.get_col_idx(col_name);
            if (idx != -1) {
                col_idx.push_back(idx);
                Z.col_names["poi*" + col_name] = count;
                count++;
            }
        }

        arma::mat interaction_mat = conv_to<mat>::from(df.data.cols(conv_to<uvec>::from(col_idx)));
        Z.data = arma::join_rows(Z.data, interaction_mat);
    }
}

template <typename T>
T convert_to(vec input);

template <>
std::vector<std::string> convert_to(vec input) {
  std::vector<std::string> output(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = std::to_string(input[i]);
  }
  return output;
}

template <>
uvec convert_to(vec input) {
  uvec output(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = static_cast<unsigned int>(input[i]);
  }
  return output;
}

template <typename T>
std::vector<T> fr_unique(const rowvec &vec) {
  std::vector<T> result(vec);
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}


FRMatrix create_design_matrix(
    FRMatrix &df,
    std::vector<Covariate> covariates,
    bool no_intercept=false,
    double colinearity_rsq=0.99
) {
    FRMatrix X;
    X.row_names = df.row_names;

    if (!no_intercept) {
        X.data = mat(df.data.n_rows, 1, fill::ones);
        X.col_names["Intercept"] = 0;
    } else {
        X.data = mat(df.data.n_rows, 0, fill::zeros);
    }
    if (covariates.empty()) return X;

    for(auto &cv : covariates) {
      cv.add_to_matrix(df, X, colinearity_rsq);
    }
    Rcpp::Rcout << "Created design matrix" << std::endl;
    return X;
}

    // if (!std::all_of(
    //     covariates.begin(),
    //     covariates.end(),
    //     [&df](const std::string &cov) { return std::find(df.col_names.begin(), df.col_names.end(), cov) != df.col_names.end(); }
    //     )
    // ) {
    //     stop("invalid covariates argument");
    // }

    // if (!std::all_of(
    //         covariate_type.begin(),
    //         covariate_type.end(),
    //         [](const std::string& kv) { return kv == "numeric" || kv == "categorical" || kv == "count"; }
    //     )
    // ) {
    //     stop("invalid covariate.type");
    // }

    // std::vector<std::string> ccov;
    // for (const auto &cv_type : covariate_type) {
    //     if (cv_type == "categorical") {
    //         size_t index = &cv_type - &covariate_type[0];
    //         ccov.push_back(covariates[index]);
    //     }
    // }

    // if (covariate_levels.empty()) {
    //     int count = 0;
    //     for(const std::string& cv : covariates) {
    //         int cv_idx = std::find(df.col_names.begin(), df.col_names.end(), cv) - df.col_names.begin();
    //         covariate_levels[count] = conv_to<std::string>::from(unique(df.data.col(cv_idx)));
    //         count++;
    //     }
    // }
    //
    // if(covariate_ref_level.empty()) {
    //     for(int i = 0; i < ccov.size(); i++) {
    //         covariate_ref_level[i] = covariate_levels[i];
    //     }
    // }

    // std::map<std::string, bool> covariate_standardize_map;
    // for (const auto &cv : covariates) {
    //     covariate_standardize_map[cv] = false;
    // }

    // for (const auto &cv: covariate_standardize) {
    //     if (covariate_standardize_map.count(cv) > 0) {
    //         covariate_standardize_map[cv] = true;
    //     }
    // }

    // for(const auto &cv : ccov) {
    //     int cv_idx = std::find(df.col_names.begin(), df.col_names.end(), cv) - df.col_names.begin();
    //     uvec not_na_idx = find_finite(df.data.col(cv_idx));
    //     auto c_idx = std::find(covariate_levels.begin(), covariate_levels.end(), cv);
    //     if(!all_finite(intersect(conv_to<uvec>::from(df.data.col(cv_idx)), conv_to<uvec>::from(covariate_levels[c_idx])))) {
    //         stop("invalid covariate.levels");
    //     }
    //     if(std::find(covariate_levels[cv].begin(), covariate_levels[cv].end(), covariate_ref_level[cv]) == covariate_levels[cv].end()) {
    //         stop("invalid covariate.ref.level");
    //     }
    // }
//
// void add_covar(
//     FRMatrix &df,
//     std::string cv,
//     std::string cv_type,
//     bool standardize,
//     FRMatrix X,
//     std::vector<std::string> levels,
//     std::string ref_level,
//     double colinearity_rsq
// ) {
//   unsigned int nS = df.data.n_rows;
//   unsigned int col_idx = std::distance(df.col_names.begin(), std::find(df.col_names.begin(), df.col_names.end(), cv));
//
//   std::map<std::string, vec> candidate_cols;
//
//   if (cv_type == "numeric") {
//     candidate_cols[cv] = df.data.col(col_idx);
//   } else {
//     if (levels.empty()) {
//
//       levels = conv_to<std::vector<std::string>>::from(unique(df.data.col(col_idx)));
//     }
//     if (ref_level.empty()) {
//       ref_level = levels[0];
//     }
//     if (std::find(levels.begin(), levels.end(), ref_level) == levels.end()) {
//       stop("ref.level must be one of unique values of var");
//     }
//
//     levels.erase(std::remove(levels.begin(), levels.end(), ref_level), levels.end());
//     if (X.data.n_cols == 0) {
//       levels.insert(levels.begin(), ref_level);
//     }
//
//     for (const auto &lev : levels) {
//       std::string col_name = cv + "(" + lev + (X.data.n_cols == 0 ? "" : " vs. " + ref_level) + ")";
//       candidate_cols[col_name] = conv_to<vec>::from(df.data.col(col_idx) == std::stod(lev));
//     }
//
//     if (standardize) {
//       stop("standardization of non-numeric variable not permitted");
//     }
//   }
//
//   if (standardize) {
//     for (auto &col : candidate_cols) {
//       double mean_val = mean(col.second);
//       double sd_val = stddev(col.second);
//       col.second = (col.second - mean_val) / sd_val;
//     }
//   }
//
//   std::vector<bool> retain(candidate_cols.size(), false);
//   unsigned int index = 0;
//
//   for (const auto &cc : candidate_cols) {
//     if (X.data.n_cols == 0 && index == 0) {
//       retain[index] = true;
//       index++;
//       continue;
//     }
//
//     vec y = cc.second;
//     mat Z = join_horiz(X.data, mat(nS, 1, fill::zeros));
//     uvec not_na_idx = find_finite(y);
//
//     Z = Z.rows(not_na_idx);
//     y = y.rows(not_na_idx);
//
//     double ssa = sum(square(y));
//     double sse = sum(square(y - Z * pinv(Z.t() * Z) * (Z.t() * y)));
//     double rsquared = 1 - sse / ssa;
//
//     if (rsquared <= colinearity_rsq) {
//       retain[index] = true;
//     } else {
//       Rcout << cc.first << " was not added to the design matrix due to the potential of colinearity\n";
//     }
//
//     index++;
//   }
//
//   index = 0;
//   for (const auto &cc : candidate_cols) {
//     if (retain[index]) {
//       X.data = join_horiz(X.data, cc.second);
//     }
//     index++;
//   }
// }


