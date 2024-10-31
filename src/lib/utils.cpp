#include <utils.h>

/**
 * @brief Function to check if an R interrupt has occurred.
 *
 * This function is called internally to check for user interrupts in R.
 */
void chkIntFn(void *dummy) { R_CheckUserInterrupt(); };

/**
 * @brief Checks if the user has triggered an interrupt.
 *
 * This function checks for user interrupts and throws an exception if one is
 * detected, stopping the current process.
 */
void checkInterrupt()
{
  if (R_ToplevelExec(chkIntFn, NULL) == FALSE)
  {
    Rcpp::stop("Received user interrupt. Stopping FastReg...");
  }
};

/**
 * @brief Checks if a directory exists.
 *
 * @param path The path to the directory.
 * @return true if the directory exists, false otherwise.
 */
bool dir_exists(const std::string &path)
{
  fs::path directory(path);
  return fs::is_directory(directory);
};
/**
 * @brief Deletes a directory and its contents.
 *
 * @param path The path to the directory to be deleted.
 */
void delete_dir(const std::string &path)
{
  fs::path directory(path);

  if (fs::exists(directory) && fs::is_directory(directory))
  {
    fs::remove_all(directory);
    Rcpp::Rcout << "Directory deleted: " << path << std::endl;
  }
  else
  {
    Rcpp::Rcout << "Directory does not exist: " << path << std::endl;
  }
};
/**
 * @brief Checks if a value is present in a container.
 *
 * @tparam T The type of the value.
 * @tparam U The type of the container.
 * @param value The value to search for.
 * @param container The container in which to search.
 * @return true if the value is found in the container, false otherwise.
 */

template <typename T, typename U>
bool isin(const T &value, const U &container)
{
  return std::find(std::begin(container), std::end(container), value) !=
         std::end(container);
};
/**
 * @brief Finds the intersection of two sets of row names.
 *
 * @param a The first vector of row names.
 * @param b The second vector of row names.
 * @return A vector containing the intersection of the two sets of row names.
 */
std::vector<std::string> intersect_row_names(const std::vector<std::string> &a,
                                             const std::vector<std::string> &b)
{
  std::vector<std::string> result;
  result.reserve(std::min(a.size(), b.size()));
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::back_inserter(result));
  return result;
};
/**
 * @brief Computes the set difference between two sets of row names.
 *
 * @param a The first vector of row names.
 * @param b The second vector of row names.
 * @return A vector containing the elements of `a` that are not in `b`.
 */
std::vector<std::string> set_diff(const std::vector<std::string> &a,
                                  const std::vector<std::string> &b)
{
  std::vector<std::string> result;
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::back_inserter(result));
  return result;
};
/**
 * @brief Transforms the Predictors of Interest (POI) matrix based on the
 * specified effect type.
 *
 * This function modifies the POI matrix by applying transformations for
 * dominant and recessive effect types.
 *
 * @param G The POI matrix to transform.
 * @param effect_type The type of effect to apply ("additive", "dominant", or
 * "recessive").
 */
void transform_poi(FRMatrix &G, std::string effect_type)
{
  if (effect_type == "dominant")
  {
    G.data.elem(find(G.data == 2)).fill(1);
  }
  else if (effect_type == "recessive")
  {
    G.data.elem(find(G.data == 1)).fill(0);
    G.data.elem(find(G.data == 2)).fill(1);
  }
};

/**
 * @brief Filters the POI matrix based on Minor Allele Frequency (MAF) and
 * Hardy-Weinberg Equilibrium (HWE) thresholds.
 *
 * This function returns a matrix with filtered POI data and associated
 * statistics such as MAF and HWE p-values.
 *
 * @param G The POI matrix to filter.
 * @param maf_threshold The MAF threshold.
 * @param hwe_threshold The HWE threshold.
 * @return A filtered FRMatrix containing the POI data and statistics.
 */
FRMatrix filter_poi(FRMatrix &G, double maf_threshold = 0.01,
                    double hwe_threshold = 0.05)
{

  FRMatrix res;
  res.col_names = G.col_names;
  res.row_names = {{"Freq (a)", 0}, {"Freq (b)", 1}, {"MAF", 2}, {"HWE Chisq", 3}, {"HWE Pvalue", 4}, {"keep", 5}};
  res.data = arma::fmat(res.row_names.size(), G.data.n_cols, arma::fill::ones);

  umat isnan_mat = umat(G.data.n_rows, G.data.n_cols, arma::fill::zeros);
  for (uword col = 0; col < G.data.n_cols; ++col)
  {
    uvec nonfinite_indices = arma::find_nonfinite(G.data.col(col));
    // isnan_mat.col(col).
    for (uword i = 0; i < nonfinite_indices.n_elem; ++i)
    {
      isnan_mat(nonfinite_indices(i), col) = 1;
    }
  }

  FRMatrix G_copy = G;
  G_copy.data.replace(datum::nan, 0);

  frowvec nS = conv_to<frowvec>::from(G.data.n_rows - arma::sum(isnan_mat, 0));
  // Rcpp::Rcout << nS << " G n_cols: " << G_copy.data.n_cols << std::endl;
  frowvec nS_2 = 2 * nS;
  frowvec a_freq = arma::sum(G_copy.data, 0) / (nS_2);
  frowvec b_freq = 1 - a_freq;

  frowvec maf_freq = min(a_freq, b_freq);
  // Rcpp::Rcout << " maf_freq: " << maf_freq << std::endl;
  frowvec aa_of = conv_to<frowvec>::from(arma::sum(G.data == 0, 0));
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

  return res;
};
/**
 * @brief Creates a matrix of interactions between POIs and covariates.
 *
 * @param df The FRMatrix containing the data.
 * @param poi_covar_interactions A vector of strings specifying the
 * POI-covariate interactions.
 * @param Z The matrix where the interactions will be stored.
 */
void create_Z_matrix(FRMatrix &df,
                     const std::vector<std::string> &poi_covar_interactions,
                     FRMatrix &Z)
{
  std::vector<std::string> interaction_cols;
  if (poi_covar_interactions.size() == 1 && poi_covar_interactions[0].empty())
  {
    return;
  }
  for (const auto &interaction : poi_covar_interactions)
  {
    for (const auto &col_name : df.col_names)
    {
      if (col_name.first.find(interaction) != std::string::npos)
      {
        interaction_cols.push_back(col_name.first);
      }
    }
  }

  if (!interaction_cols.empty())
  {
    Rcpp::Rcout << "Found poi covar interactions" << std::endl;
    std::vector<int> col_idx;
    int count = 0;
    for (const auto &col_name : interaction_cols)
    {
      int idx = df.get_col_idx(col_name);
      if (idx != -1)
      {
        col_idx.push_back(idx);
        Z.col_names["poi*" + col_name] = count;
        count++;
      }
    }

    arma::fmat interaction_mat =
        conv_to<fmat>::from(df.data.cols(conv_to<uvec>::from(col_idx)));
    Z.data = arma::join_rows(Z.data, interaction_mat);
  }
};
/**
 * @brief Returns the unique elements of a floating point vector.
 *
 * @tparam T The type of the elements.
 * @param vec The vector to process.
 * @return A vector containing the unique elements of the input vector.
 */
template <typename T>
std::vector<T> fr_unique(const frowvec &vec)
{
  std::vector<T> result(vec);
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
};
/**
 * @brief Creates a design matrix from covariates and filters based on
 * collinearity.
 *
 * This function builds the design matrix from covariates and optionally adds an
 * intercept term.
 *
 * @param df The FRMatrix containing the data.
 * @param covariates The vector of covariates to include in the design matrix.
 * @param no_intercept Whether to exclude the intercept term.
 * @param colinearity_rsq The collinearity threshold for filtering covariates.
 * @return An FRMatrix representing the design matrix.
 */

FRMatrix create_design_matrix(FRMatrix &df, std::vector<Covariate> covariates,
                              bool no_intercept, double colinearity_rsq)
{
  FRMatrix X;
  X.row_names = df.row_names;

  if (!no_intercept)
  {
    X.data = fmat(df.data.n_rows, 1, fill::ones);
    X.col_names["Intercept"] = 0;
  }
  else
  {
    X.data = fmat(df.data.n_rows, 0, fill::zeros);
  }
  if (covariates.empty())
    return X;

  for (auto &cv : covariates)
  {
    cv.add_to_matrix(df, X, colinearity_rsq);
  }
  return X;
};

/**
 * @brief Checks if a string contains only whitespace characters.
 *
 * @param str The string to check.
 * @return true if the string contains only whitespace, false otherwise.
 */
bool isWhitespace(const std::string &str)
{
  return std::all_of(str.begin(), str.end(),
                     [](char c)
                     { return std::isspace(c); });
};