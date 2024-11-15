#define ARMA_WARN_LEVEL 0

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "BEDMatrix.h"
#include <blas_library_manager.h>
#include <regression.h>
#include <vla_result.h>

#if !defined(__APPLE__) && !defined(__MACH__)
#include <omp.h>
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/wait.h>
#endif

using namespace Rcpp;
using namespace arma;

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#if __cplusplus >= 201703L && __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#endif

#include "BEDMatrix.h"
#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include <RcppEigen.h>

void delete_dir(const std::string &path) {
  fs::path directory(path);

  if (fs::exists(directory) && fs::is_directory(directory)) {
    fs::remove_all(directory);
    Rcpp::Rcout << "Directory deleted: " << path << std::endl;
  } else {
    Rcpp::Rcout << "Directory does not exist: " << path << std::endl;
  }
}

double recode_genotype2(int genotype) {
  switch (genotype) {
  case 0:
    return 2.0; // two copies of A1
  case 3:
    return 0.0; // zero copies of A1
  case 2:
    return 1.0; // one copy of A1
  default:
    return arma::datum::nan; // invalid/missing
  }
  return arma::datum::nan; // invalid/missing
}

float recode_genotype2f(int genotype) {
  switch (genotype) {
  case 0:
    return 2.0; // two copies of A1
  case 3:
    return 0.0; // zero copies of A1
  case 2:
    return 1.0; // one copy of A1
  default:
    return arma::datum::nan; // invalid/missing
  }
  return arma::datum::nan; // invalid/missing
}

arma::mat scanBEDMatrix(SEXP xptr, arma::ivec i, arma::ivec j) {
  struct BEDMatrix *state =
      static_cast<struct BEDMatrix *>(R_ExternalPtrAddr(xptr));
  if (state == NULL) {
    Rf_error("BEDMatrix instance has been unmapped.");
  }
  int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
  int ni = i.n_rows;
  int nj = j.n_rows;
  arma::mat out(ni, nj);
  for (int cj = 0; cj < nj; cj++) {
    int jj = j[cj];
    for (int ci = 0; ci < ni; ci++) {
      int ii = i[ci];
      out(ci, cj) = recode_genotype2(extract_genotype_cartesian(
          state->data, ii - 1, jj - 1, num_bytes_per_variant));
    }
  }
  return out;
}

arma::fmat scanBEDMatrixf(SEXP xptr, arma::ivec i, arma::ivec j) {
  struct BEDMatrix *state =
      static_cast<struct BEDMatrix *>(R_ExternalPtrAddr(xptr));
  if (state == NULL) {
    Rf_error("BEDMatrix instance has been unmapped.");
  }
  int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
  int ni = i.n_rows;
  int nj = j.n_rows;
  arma::fmat out(ni, nj);

  for (int cj = 0; cj < nj; cj++) {
    int jj = j[cj];
    for (int ci = 0; ci < ni; ci++) {
      int ii = i[ci];
      out(ci, cj) = recode_genotype2f(extract_genotype_cartesian(
          state->data, ii - 1, jj - 1, num_bytes_per_variant));
    }
  }
  return out;
}
/// @brief Standardize (center/scale) X and remove features with low
/// explainability; it is essentially PCA dimension reduction
/// @param x Matrix with zero'ed out columns for X/Y missing
/// @param dev_to_acc Amount of explained variance (summed Eigenvalues) to keep
/// @param N Number of non-zero rows (number of obs)
/// @return A clean X matrix for further calculations
arma::mat standardize_and_derank(arma::mat x, double dev_to_acc = 0.95,
                                 double N = 1.0) {
  // scale (subtract mean and divide by SD) in place
  x.each_col([&N](arma::vec &y) {
    y -= (arma::accu(y) / N);
    y /= sqrt(arma::accu(arma::square(y)) / (N - 1));
  });
  // get principal components
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::princomp(coeff, score, latent, x);
  // coeff = rotation in R, latent is the eigenvalues (sdev^2)
  //  arma::uvec indices = sort_index(-1.0 * latent);
  //  arma::vec sorted = sort(latent);
  arma::uvec ans = arma::find(cumsum(latent) / accu(latent) > dev_to_acc);
  return x * coeff.cols(0, ans[0]);
}

/// @brief Standardize (center/scale) X and remove features with low
/// explainability; it is essentially PCA dimension reduction WITH PRINTING
/// @param x Matrix with zero'ed out columns for X/Y missing
/// @param dir Directory
/// @param namee Name of Y
/// @param suffix Suffix from above just in case
/// @param dev_to_acc Amount of explained variance (summed Eigenvalues) to keep
/// @param N Number of non-zero rows (number of obs)
/// @return A clean X matrix for further calculations
//@export
// [[Rcpp::export]]
arma::mat standardize_and_derank_print(arma::mat x, const std::string &dir,
                                       const std::string &namee,
                                       const std::string &suffix,
                                       double dev_to_acc = 0.95,
                                       double N = 1.0) {
  // scale (subtract mean and divide by SD) in place

  arma::mat ans(2, x.n_cols);
  ans.row(0) = (arma::sum(x, 0.0) / N);
  x.each_row() -= ans.row(0);
  ans.row(1) = sqrt(arma::sum(arma::square(x), 0.0) / (N - 1));
  // ensure it doesn't go zero and divide by zero which introduces NAs
  ans.row(1).clamp(1e-8, 1e99);
  x.each_row() /= ans.row(1);
  // get principal components
  arma::mat coeff;  // x.n_cols, x.n_cols
  arma::mat score;  // x.n_rows, x.n_cols
  arma::vec latent; // x.n_cols
  arma::princomp(coeff, score, latent, x);

  arma::vec cumulative_var = cumsum(latent) / accu(latent);
  arma::uvec temp = arma::find(cumulative_var > dev_to_acc);

  std::stringstream ss;
  // ss << dir << "/" << namee << "_" << suffix << "_metadata.tsv";
  fs::create_directory(dir);
  fs::create_directory(dir + "/" + namee);
  ss << dir << "/" << namee;
  ss << "/PCA_metadata.tsv";
  std::string result_file = ss.str();
  std::ofstream outfile;
  bool need_writing = !(fs::exists(result_file));
  //if we need to write it, do it
  std::stringstream buffer;
  if(need_writing){
    outfile.open(result_file);
    outfile << std::fixed << std::setprecision(8);
    // std::stringstream buffer;
    // add column headers
    buffer << "Means"
          << "\t"
          << "SDs";
  }
  // if not all PCs, use only subset
  if (temp.n_elem > 0) {
    if(need_writing){
      arma::mat output = arma::join_rows(ans.t(), coeff.cols(0, temp[0]));
      for (uword ii = 0; ii < (temp[0] + 1); ii++) {
        buffer << "\t"
              << "Loading" << (ii + 1);
      }
      buffer << std::endl;
      for (uword j = 0; j < output.n_rows; j++) {
        for (uword k = 0; k < (output.n_cols - 1); k++) {
          buffer << output(j, k) << "\t";
        }
        buffer << output(j, output.n_cols - 1);
        buffer << std::endl;
      }
      buffer << 0 << "\t" << 0;
      for (uword ii = 0; ii < (temp[0] + 1); ii++) {
        buffer << "\t" << cumulative_var[ii];
      }
      buffer << std::endl;
      outfile << buffer.str();
      outfile.close();
    }
    return x * coeff.cols(0, temp[0]);
  } else { // needed to prevent crashes
    if(need_writing){
      arma::mat output = arma::join_rows(ans.t(), coeff);
      for (uword ii = 0; ii < coeff.n_cols; ii++) {
        buffer << "\t"
              << "Loading" << int(ii);
      }
      buffer << std::endl;
      for (uword j = 0; j < output.n_rows; j++) {
        for (uword k = 0; k < (output.n_cols - 1); k++) {
          buffer << output(j, k) << "\t";
        }
        buffer << output(j, output.n_cols - 1);
        buffer << std::endl;
      }
      buffer << 0 << "\t" << 0;
      for (uword ii = 0; ii < coeff.n_cols; ii++) {
        buffer << "\t" << cumulative_var[ii];
      }
      buffer << std::endl;
      outfile << buffer.str();
      outfile.close();
    }
    return x * coeff;
  }
}

// [[Rcpp::export]]
void FastVLA_logisticf(const arma::mat &Y, SEXP Gptr, const arma::ivec &v_index,
                       const arma::ivec &i_index, const arma::mat &X,
                       const int &chunk_size, const std::string &dir,
                       const std::vector<std::string> cnames,
                       const std::vector<std::string> vnames,
                       const std::string suffix, const float epss = 1e-7,
                       const double &mafthresh = 0.005,
                       const double &pca_var_explained = 0.95,
                       const int max_iter = 6, const int max_threads = 1,
                       const int max_blas_threads = 1, const bool do_pca = true,
                       const bool add_intercept = true) {

  // delete_dir(dir);
  BLASLibraryManager blas_mgr;
  blas_mgr.detect_lib();
  int cur_blas_threads = blas_mgr.get_num_threads();
  if (cur_blas_threads != max_blas_threads) {
    blas_mgr.set_num_threads(max_blas_threads);
  }
  int total_variants = v_index.n_elem;
  int num_chunks = total_variants / chunk_size;
  if (total_variants % chunk_size != 0) {
    num_chunks += 1;
  }

  // get boolean inclusion vectors
  arma::fmat xy_bools(Y.n_rows, Y.n_cols);
  uvec badX = find_nan(arma::sum(X, 1.0));
  // put in each column booleans
  for (uword i = 0; i < Y.n_cols; i++) {
    // get X NANs for removal
    arma::fcolvec x_bools(Y.n_rows, fill::ones);
    x_bools.elem(badX).fill(0.0);
    // add Y NANs for removal
    arma::uvec badY = find_nan(Y.col(i));
    x_bools.elem(badY).fill(0.0);
    xy_bools.col(i) = x_bools;
  }

#if !defined(__APPLE__) && !defined(__MACH__)
  omp_set_dynamic(0);
  omp_set_num_threads(max_threads);
#endif
  // cnames = prefix per pheno - make into folder for each pheno
  // vnames = 1 vname per poi - aka rownames

  for (int i = 0; i < Y.n_cols; i++) {
    arma::mat cov_d = X;
    arma::ivec _index = i_index;
    arma::fmat pheno = arma::conv_to<arma::fcolvec>::from(Y.col(i));
    arma::uvec nan_idx = arma::find(xy_bools.col(i) == 0.0);
    _index.shed_rows(nan_idx);
    pheno.shed_rows(nan_idx);
    cov_d.shed_rows(nan_idx);

    // arma::fmat cov = arma::conv_to<arma::fmat>::from(
    //     standardize_and_derank(cov_d, pca_var_explained,
    //     double(cov_d.n_rows)));

    arma::fmat cov;
    if (do_pca) {
      auto start = std::chrono::high_resolution_clock::now();
      cov = arma::conv_to<arma::fmat>::from(standardize_and_derank_print(
          cov_d, dir, cnames[i], suffix, pca_var_explained,
          double(cov_d.n_rows)));

      auto end = std::chrono::high_resolution_clock::now();
      double timing =
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                        start)
              .count();
      Rcpp::Rcout << "PCA ananlysis timing: " << timing << "ms" << std::endl;
    } else {
      cov = arma::conv_to<arma::fmat>::from(cov_d);
    }
    arma::fmat no_interactions = arma::fmat(cov.n_rows, 1, arma::fill::ones);
    arma::fmat interactions = arma::join_rows(no_interactions, cov);
    arma::fmat interactions_sqrd = arma::join_rows(no_interactions, cov);
    // add column of ones if required (same as interactions matrix now)
    if (add_intercept) {
      cov = interactions;
    }

    std::time_t t = std::time(nullptr);
    std::tm *lt = std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(lt, "%d-%m-%Y %H-%M-%S");
    std::string local_time = oss.str();
    for (int chunk = 0; chunk < num_chunks; ++chunk) {
      int start_idx = chunk * chunk_size;
      int end_idx = std::min(start_idx + chunk_size - 1, total_variants - 1);

      arma::ivec current_variant_indices = v_index.subvec(start_idx, end_idx);

      // Extract corresponding variant names
      std::vector<std::string> variant_names_chunk(
          vnames.begin() + start_idx, vnames.begin() + end_idx + 1);
      arma::fmat G = scanBEDMatrixf(Gptr, _index, current_variant_indices);
      VLAResultf res(cov, G, no_interactions, interactions, interactions_sqrd,
                     local_time);

      LogisticRegression regression;
      regression.run_vla(cov, pheno, G, res, max_iter, true, mafthresh, epss);
      res.write_to_file(dir, suffix, cnames[i], variant_names_chunk);
      if (chunk == 0) {
        res.write_headers(dir, suffix, cnames[i], variant_names_chunk);
      }
    }
  }
  blas_mgr.set_num_threads(cur_blas_threads);
}

// [[Rcpp::export]]
void FastVLA_logistic(const arma::mat &Y, SEXP Gptr, const arma::ivec &v_index,
                      const arma::ivec &i_index, const arma::mat &X,
                      const int &chunk_size, const std::string &dir,
                      const std::vector<std::string> cnames,
                      const std::vector<std::string> vnames,
                      const std::string suffix, const double epss = 1e-7,
                      const double &mafthresh = 0.005,
                      const double &pca_var_explained = 0.950,
                      const int max_iter = 6, const int max_threads = 1,
                      const int max_blas_threads = 1, const bool do_pca = true,
                      const bool add_intercept = true) {

  // delete_dir(dir);
  BLASLibraryManager blas_mgr;
  blas_mgr.detect_lib();
  int cur_blas_threads = blas_mgr.get_num_threads();
  if (cur_blas_threads != max_blas_threads) {
    blas_mgr.set_num_threads(max_blas_threads);
  }

  int total_variants = v_index.n_elem;
  int num_chunks = total_variants / chunk_size;
  if (total_variants % chunk_size != 0) {
    num_chunks += 1;
  }

  // get boolean inclusion vectors
  arma::mat xy_bools(Y.n_rows, Y.n_cols);
  uvec badX = find_nan(arma::sum(X, 1.0));
  // put in each column booleans
  for (uword i = 0; i < Y.n_cols; i++) {
    // get X NANs for removal
    arma::vec x_bools(Y.n_rows, fill::ones);
    x_bools.elem(badX).fill(0.0);
    // add Y NANs for removal
    arma::uvec badY = find_nan(Y.col(i));
    x_bools.elem(badY).fill(0.0);
    xy_bools.col(i) = x_bools;
  }

#if !defined(__APPLE__) && !defined(__MACH__)
  omp_set_dynamic(0);
  omp_set_num_threads(max_threads);
#endif
  // cnames = prefix per pheno - make into folder for each pheno
  // vnames = 1 vname per poi - aka rownames

  for (int i = 0; i < Y.n_cols; i++) {
    arma::mat cov = X;
    arma::ivec _index = i_index;
    arma::mat pheno = Y.col(i);
    arma::uvec nan_idx = arma::find(xy_bools.col(i) == 0.0);
    _index.shed_rows(nan_idx);
    pheno.shed_rows(nan_idx);
    cov.shed_rows(nan_idx);

    // cov = standardize_and_derank(cov, pca_var_explained, double(cov.n_rows));
    if (do_pca) {
      auto start = std::chrono::high_resolution_clock::now();
      cov = standardize_and_derank_print(cov, dir, cnames[i], suffix,
                                         pca_var_explained, double(cov.n_rows));

      auto end = std::chrono::high_resolution_clock::now();
      double timing =
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                        start)
              .count();
      Rcpp::Rcout << "PCA ananlysis timing: " << timing << "ms" << std::endl;
    }
    arma::mat no_interactions = arma::mat(cov.n_rows, 1, arma::fill::ones);
    arma::mat interactions = arma::join_rows(no_interactions, cov);
    arma::mat interactions_sqrd = arma::join_rows(no_interactions, cov);
    // add column of ones if required (same as interactions matrix now)
    if (add_intercept) {
      cov = interactions;
    }

    std::time_t t = std::time(nullptr);
    std::tm *lt = std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(lt, "%d-%m-%Y %H-%M-%S");
    std::string local_time = oss.str();
    for (int chunk = 0; chunk < num_chunks; ++chunk) {
      int start_idx = chunk * chunk_size;
      int end_idx = std::min(start_idx + chunk_size - 1, total_variants - 1);

      arma::ivec current_variant_indices = v_index.subvec(start_idx, end_idx);

      // Extract corresponding variant names
      std::vector<std::string> variant_names_chunk(
          vnames.begin() + start_idx, vnames.begin() + end_idx + 1);
      arma::mat G = scanBEDMatrix(Gptr, _index, current_variant_indices);
      VLAResult res(cov, G, no_interactions, interactions, interactions_sqrd,
                    local_time);
      LogisticRegression regression;
      regression.run_vla(cov, pheno, G, res, max_iter, true, mafthresh, epss);
      res.write_to_file(dir, suffix, cnames[i], variant_names_chunk);
      if (chunk == 0) {
        res.write_headers(dir, suffix, cnames[i], variant_names_chunk);
      }
    }
  }
  blas_mgr.set_num_threads(cur_blas_threads);
}

arma::mat manual_svd(arma::mat x, double &rankk, double epsilonn = 0.0001) {
  arma::mat U(x.n_cols, x.n_cols);
  arma::mat V(x.n_cols, x.n_cols);
  arma::vec s(x.n_cols);
  arma::svd(U, s, V, x.t() * x);
  arma::vec ones(x.n_cols, arma::fill::ones);
  // make sure epsilon >= epsilon in armadillo
  epsilonn *= s.max();
  // epsilonn = std::max(epsilonn, datum::eps);
  // ones.elem(arma::find(s < (epsilonn * max(s)))).fill(0.0); //remove bad
  // eigenvals
  ones.elem(arma::find(s < epsilonn)).fill(0.0); // remove bad eigenvals

  s = ones / s;
  // s.clean(epsilonn * max(s));
  // ones.elem(s == 0).fill(0.0);
  rankk = arma::accu(ones);
  return V * arma::diagmat(s) * U.t();
}

arma::mat manual_svd2(arma::mat x, double &rankk, double epsilonn = 0.0001) {
  arma::mat U;
  arma::mat V;
  arma::vec s;
  arma::svd_econ(U, s, V, x.t() * x, "left");
  // s = ones / s;
  // s.clean(epsilonn * max(s));
  // ones.elem(s == 0).fill(0.0);
  arma::uvec temp = find(s < epsilonn * s.max());
  s.elem(temp).fill(arma::datum::inf);
  rankk = (double)(s.n_elem - temp.n_elem);
  return U * arma::diagmat(1.0 / s) * U.t();
}

/// @brief Outputs 1-D estimate of beta (for G) in various regressions
/// @param value Row number
/// @param offset Starting place for outputting this regression
/// @param df Degrees of freedom
/// @param estimat Estimate of beta
/// @param MSE Mean Squared Error
/// @param inv_var Inverse of variance (for SE of beta)
/// @param output Output matrix
/// @param logg10 -log10 constant
/// @param logg2 log2 constant
void output_1d(const uword &value, uword offset, double &df, double &estimat,
               double &MSE, double inv_var, mat &output, const double &logg10,
               const double &logg2) {
  output(value, offset) = df;
  output(value, offset + 1) = estimat;
  double se = sqrt(MSE / inv_var);
  output(value, offset + 2) = se;
  output(value, offset + 3) =
      std::max(0.0, ((R::pt(abs(estimat / se), df, 0, 1) + logg2) / logg10));
  // output(value, offset + 3) =
  //     (R::pt(abs(estimat / se), df, 0, 1) + logg2) / logg10;
  return;
}

/// @brief Outputs 2-D estimates of betas (for G, G^2 by construction) in
/// regressions
/// @param value Row number
/// @param offset Starting place for outputting this regression
/// @param df Degrees of freedom
/// @param estimat Estimate of betas
/// @param MSE Mean Squared Error
/// @param var_mat Variance matrix to rescale for SE of Betas
/// @param output Output matrix
/// @param logg10 -log10 constant
/// @param logg2 log2 constant
void output_2d(const uword &value, uword offset, double &df, vec &estimat,
               double &MSE, mat &var_mat, mat &output, const double &logg10,
               const double &logg2) {
  output(value, offset) = df;
  output(value, offset + 1) = estimat[0];
  double se = sqrt(MSE * var_mat(0, 0));
  output(value, offset + 2) = se;
  output(value, offset + 3) =
      std::max(0.0, (R::pt(abs(estimat[0] / se), df, 0, 1) + logg2) / logg10);
  se = sqrt(MSE * var_mat(1, 1));
  output(value, offset + 4) = estimat[1];
  output(value, offset + 5) = se;
  output(value, offset + 6) =
      std::max(0.0, (R::pt(abs(estimat[1] / se), df, 0, 1) + logg2) / logg10);
  return;
}

/// @brief Regresses X out of Y,G,H
/// @param yadj Y before X removal
/// @param gadj G before G removal
/// @param hadj H before G removal
/// @param temp_matmul Memory placeholder for Y,G,H
/// @param x covariates
/// @param Xadj Inverse of X^T X (from manual_svd2 output)
/// @param w Inclusion mask
void regress_out_X(vec &yadj, vec &gadj, vec &hadj, mat &temp_matmul, mat &x,
                   mat &Xadj, vec &w) {
  temp_matmul = join_rows(yadj, gadj, hadj);
  temp_matmul -= x * (Xadj * temp_matmul);
  // make sure zeros for irrelevant rows
  temp_matmul.each_col() %= w;
  yadj = temp_matmul.col(0);
  gadj = temp_matmul.col(1);
  hadj = temp_matmul.col(2);
  return;
}

/// @brief Given a file name, prints output to file
/// @param file_name File name
/// @param output Output matrix to print to file
void printer(std::string &file_name, mat &output, const bool &append,
             const std::vector<std::string> &vnames) {
  std::ofstream outfile;
  if (append) {
    outfile.open(file_name, std::ios_base::app);
  } else {
    outfile.open(file_name);
  }
  outfile << std::fixed << std::setprecision(8);
  //     //output.save(outfile);
  std::stringstream buffer;
  for (uword j = 0; j < output.n_rows; j++) {
    buffer << vnames[j];
    for (uword k = 0; k < (output.n_cols); k++) {
      buffer << "\t" << output(j, k);
    }
    buffer << std::endl;
    // buffer << output(j, output.n_cols - 1);
    // buffer << std::endl;
  }
  outfile << buffer.str();
  outfile.close();
}

/// @brief Gets the correct column of Y, X, G, and w (the mask)
/// @param i Y column
/// @param value G column
/// @param y Y placeholder to fill
/// @param g G col placeholder to fill
/// @param x X mat placeholder to fill
/// @param w Mask placeholder to fill
/// @param Y Y matrix
/// @param G G matrix
/// @param X X matrix
/// @param mask Mask of Y/X
void get_ygxw(uword &i, const uword &value, vec &y, vec &g, mat &x, vec &w,
              const mat &Y, const mat &G, const mat &X, vec &mask) {
  y = Y.col(i);
  x = X;
  g = G.col(value);

  w = mask;
  w(find_nan(g)).fill(0.0);
  // w = ws.col(value);

  uvec bad_inds = find(w == 0);

  g.elem(bad_inds).fill(0.0);
  y.elem(bad_inds).fill(0.0);
  x.rows(bad_inds).fill(0.0);
  return;
}

/// @brief Workhorse function of VLA regressions
/// @param Y Phenotypes/outcomes matrix
/// @param G Matrix of POIs (subset by the calling function)
/// @param X Matrix of covariates
/// @param dir Directory for output
/// @param Ybad Matrix of indicators if X/Y is NaN (0) or should be included (1)
/// @param cnames Subfolder names of Y columns
/// @param vnames Rownames for output (POIs)
/// @param suffix Name to append for output in correspond dir/cnames folder
/// @param epss Epsilon threshold for SVD/inversion
/// @param mafthresh Threshold for Mean allele frequency
/// @param append Boolean -- create a new file or append to old?
/// @param pca_var_explained -- Double describing how much explained variance in
/// X needed (for rank reduction)
/// @return 1 if successful
int FastVLA_cpp_internal(const arma::mat &Y, const arma::mat &G,
                         const arma::mat &X_orig, const std::string &dir,
                         const arma::mat &Ybad,
                         const std::vector<std::string> cnames,
                         const std::vector<std::string> vnames,
                         const std::string suffix, const double epss = 1e-6,
                         const double &mafthresh = 0.01,
                         const bool &append = false,
                         const double pca_var_explained = 0.95) {

  // precreate objects
  double F;
  mat x = X_orig;       // copy of X
  arma::mat X = X_orig; // another copy of X for Deepak's rank reduction

  vec Wtemp2(Y.n_rows); // for X and Y bads
  vec w(Y.n_rows);

  // to store necessary additional vectors
  vec g_unscaled(Y.n_rows);
  vec yadj(Y.n_rows);
  vec xadj(Y.n_rows);
  vec gadj(Y.n_rows);
  vec hadj(Y.n_rows);
  mat Xadj(X.n_cols, Y.n_rows);
  vec yvQTL(Y.n_rows);

  mat temp_matmul(Y.n_rows, 3);
  mat Xg(Y.n_rows, 2);
  mat Vg(2, 2);

  vec est(2);

  vec z(Y.n_rows);

  double p, gmean, df, SST, SSE, MSE;

  double rankk;

  // for 1d
  double ans, estt;

  // constants for conversion from log pval to -log10 pval
  const double logg10 = -1.0 * log(10.0);
  const double logg2 = log(2.0);

  // output holder
  mat output(G.n_cols, 33);
  // condition checks placeholders
  double condition, caf;
  bool ones, twos;
  bool temp;
  double one_counter;
  vec ps(G.n_cols);
  vec whichcondition(G.n_cols);

  // create output directory
  fs::create_directory(dir);

  // loop over Ys
  for (uword i = 0; i < Y.n_cols; i++) {

    X = X_orig;

    // fill output with NANs
    output.fill(datum::nan);

    // create file name for this Y
    std::stringstream ss;
    ss << dir << "/" << cnames[i];
    std::string result_file = ss.str();

    // make sure result directory is there then add the file inside it
    fs::create_directory(result_file);

    ss << "/" << suffix << ".tsv";
    result_file = ss.str();

    // y adj= Y.col(i);
    // get bad X,Y joint for now
    Wtemp2 = Ybad.col(i);

    // fix x
    // first zero out rows with NAs

    arma::uvec nansxy = find(Wtemp2 == 0);
    X.rows(nansxy).fill(0.0);
    // then center/standardize/drop cols
    // X = standardize_and_derank(X, pca_var_explained,
    //                            double(X.n_rows - nansxy.n_elem));
    X = standardize_and_derank_print(X, dir, cnames[i], suffix,
                                     pca_var_explained,
                                     double(X.n_rows - nansxy.n_elem));

    // loop over POIs/genes and check if good gene, bad gene, or minimal gene
    // (only 2 unique values)
    for (uword j = 0; j < G.n_cols; j++) {
      gadj = G.col(j);
      // condition = 10.0;
      w = Wtemp2;
      w(find_nan(gadj)).fill(0.0);

      uvec bad_inds = find(w == 0);
      gadj.elem(bad_inds).fill(0.0);

      // put it into w storage
      //  ws.col(j) = w;
      p = arma::accu(w);
      gmean = dot(gadj.t(), w) / p;
      output(j, 0) = p;
      caf = gmean / 2.0;
      // changed as of 10/10 to be CAF instead of MAF
      //    output(j,1) = std::min(caf, 1.0 - caf);
      output(j, 1) = caf;
      //   whichcondition[j] = check_condition(gadj, caf, mafthresh);
      condition = -1.0 * ((caf < mafthresh) |
                          (caf > (1.0 - mafthresh))); // 1 is all NANs
      ps[j] = p;
      if ((caf < mafthresh) | (caf > (1.0 - mafthresh))) {
        whichcondition[j] = 2.0;
      } else {
        // check if ones and twos
        // ones = false;
        // twos = false;
        // one_counter = 0.0;
        // // 10/10 add a counter for number of ones;
        // for (const auto &value : gadj) {
        //   temp = (value == 1.0);
        //   ones = ones | temp;
        //   one_counter += 1.0 * (temp);
        //   twos = twos | (value == 2.0);
        // }
        ones = arma::any(gadj == 1.0);
        twos = arma::any(gadj == 2.0);
        one_counter = arma::accu(gadj == 1.0);
        whichcondition[j] = 2.0 - ones - twos;
        output(j, 2) = one_counter;
      }
    }

    // find conditions that need further calculations
    arma::uvec oneg = find(whichcondition == 1.0);
    arma::uvec good = find(whichcondition == 0.0);

    // do 1D G
#pragma omp parallel for private(x, w, g_unscaled, Xadj, temp_matmul, yadj,    \
                                 gadj, hadj, Xg, Vg, SST, SSE, MSE, df, z, F,  \
                                 est, ans, estt)
    for (const auto &value : oneg) {

      get_ygxw(i, value, yadj, gadj, x, w, Y, G, X, Wtemp2);
      try {
        Xadj = manual_svd2(x, rankk, epss);
        output(value, 3) = rankk;
        output(value, 4) = 1.0;

        // create un-adjusted copy for vQTL
        g_unscaled = gadj;

        // as of 10/3 no more scaling of G
        //  g -= gmeans[value];
        // get G^2
        hadj = arma::square(gadj) % w;

        Xadj = Xadj * x.t();

        // sGWAS regress out X
        regress_out_X(yadj, gadj, hadj, temp_matmul, x, Xadj, w);

        // common between sGWAS, fGWAS
        df = ps[value] - rankk - 1.0;
        SST = arma::accu(arma::square(yadj));

        // sGWAS: regress out G
        ans = dot(gadj.t(), gadj);
        estt = dot(gadj.t(), yadj) / ans;
        z = arma::square((yadj - gadj * estt) % w);
        SSE = arma::accu(z);
        MSE = SSE / df;

        output_1d(value, 5, df, estt, MSE, ans, output, logg10, logg2);

        // vQTL: regression on squared residuals (z) from sGWAS
        Xg = arma::join_rows(w, g_unscaled);
        Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);
        // decrease by 2 more for vQTL to account for vQTL betas [intercept + G]
        df -= 2.0;
        est = Vg * (Xg.t() * z);
        MSE = arma::accu(arma::square((z - Xg * est) % w)) / df;
        // se = sqrt(MSE * Vg(1,1));
        output_1d(value, 9, df, est[1], MSE, 1.0 / Vg(1, 1), output, logg10,
                  logg2);

        // add back the 2 DF removed for vQTL
        df += 2.0;

        // Regression 1 of VLA
        // No g^2 so no need for matrix
        // estimate betas and subtract out for SSE
        g_unscaled = gadj;
        ans = dot(gadj.t(), gadj);
        estt = dot(gadj.t(), yadj) / ans;
        gadj *= estt;
        z = arma::square((yadj - gadj) % w);
        SSE = arma::accu(z);
        MSE = SSE / df;

        output_1d(value, 13, df, estt, MSE, ans, output, logg10, logg2);

        // add in F statistic (importance of adding G to the regression)
        F = ((SST - SSE) / 1.0) / MSE;
        output(value, 20) = F;
        output(value, 21) = (R::pf(F, 1, df, 0, 1)) / logg10;

        // STEP 2: Regress squared residuals from above on X + G
        // regress out X in VLA
        hadj = (z - x * (Xadj * z)) % w;
        SST = arma::accu(arma::square(hadj));

        // update DF! Remove added X + POI estimate dfs
        df -= (rankk + 1.0);

        // regress on G and get betas
        estt = dot(g_unscaled.t(), hadj) / ans;
        g_unscaled *= estt;
        z = arma::square((hadj - g_unscaled) % w);
        SSE = arma::accu(z);
        MSE = SSE / df;

        output_1d(value, 22, df, estt, MSE, ans, output, logg10, logg2);

        // add in F-stat (importance of G in VLA variance regression) and
        // SST,SSE
        F = ((SST - SSE) / 1.0) / MSE;
        output(value, 29) = SST;
        output(value, 30) = SSE;
        output(value, 31) = F;
        output(value, 32) = (R::pf(F, 1, df, 0, 1)) / logg10;

      } catch (...) { // catch-all for issues
        output(value, 2) = 100;
      }
    }

    // do 2D G (has 0, 1, and 2 values)
#pragma omp parallel for private(x, w, g_unscaled, Xadj, temp_matmul, yadj,    \
                                 gadj, hadj, Xg, Vg, SST, SSE, MSE, df, z, F,  \
                                 est, ans, estt)
    for (const auto &value : good) {
      get_ygxw(i, value, yadj, gadj, x, w, Y, G, X, Wtemp2);
      // try to compute everything with a catch-all for failure
      try {
        Xadj = manual_svd2(x, rankk, epss);
        output(value, 3) = rankk;
        output(value, 4) = 2.0;

        // create a copy pre-regressing out X
        g_unscaled = gadj;

        // as of 10/3 no more scaling of G
        //  g -= gmeans[value];
        hadj = arma::square(gadj) % w;

        Xadj = Xadj * x.t();

        regress_out_X(yadj, gadj, hadj, temp_matmul, x, Xadj, w);

        df = ps[value] - rankk - 1.0;
        SST = arma::accu(arma::square(yadj));

        // sGWAS (standard)
        ans = dot(gadj.t(), gadj);
        estt = dot(gadj.t(), yadj) / ans;
        z = arma::square((yadj - gadj * estt) % w);
        SSE = arma::accu(z);
        MSE = SSE / df;

        // output relevant requested statistics of the regression
        output_1d(value, 5, df, estt, MSE, ans, output, logg10, logg2);

        // vQTL is 8-11 (standard for comparison)
        Xg = arma::join_rows(w, g_unscaled);
        Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);
        // decrease by 2 more for vQTL to account for vQTL betas
        df -= 2.0;
        est = Vg * (Xg.t() * z);
        MSE = arma::accu(arma::square((z - Xg * est) % w)) / df;

        output_1d(value, 9, df, est[1], MSE, 1.0 / Vg(1, 1), output, logg10,
                  logg2);

        // NEW VLA STEP 1: regress on G and G^2

        Xg = join_rows(gadj, hadj);
        Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);

        // add back 1 df from vQTL original
        df += 1.0;
        // get beta estimates
        est = Vg * (Xg.t() * yadj);
        z = arma::square((yadj - Xg * est) % w);
        SSE = arma::accu(z);
        MSE = SSE / df;
        // print output
        output_2d(value, 13, df, est, MSE, Vg, output, logg10, logg2);
        // add in STEP 1 F value
        F = ((SST - SSE) / 2.0) / MSE;
        output(value, 20) = F;
        output(value, 21) = (R::pf(F, 2, df, 0, 1)) / logg10;

        // STEP 2: Regress X out of residuals^2 for VLA
        hadj = (z - x * (Xadj * z)) % w;

        // update DF NEW TO REMOVE X AND G AND G^2 AGAIN
        df -= (rankk + 2.0);

        // get SST with only X regressed
        SST = arma::accu(arma::square(hadj));
        est = Vg * (Xg.t() * hadj);
        z = arma::square((hadj - Xg * est) % w);

        // get SSE with G and G^2 regressed as well
        SSE = arma::accu(z);
        MSE = SSE / df;

        // output relevant statistics for G and G^2
        output_2d(value, 22, df, est, MSE, Vg, output, logg10, logg2);

        // add in F for Level 2 VLA regression (and SST,SSE)
        F = ((SST - SSE) / 2.0) / MSE;
        output(value, 29) = SST;
        output(value, 30) = SSE;
        output(value, 31) = F;
        output(value, 32) = (R::pf(F, 2, df, 0, 1)) / logg10;
      } catch (...) { // catch-all for issues
        output(value, 3) = 100;
      }
    }
    // output the necessary results to file from the matrix
    printer(result_file, output, append, vnames);
  }

  return 1; // success!
}

//' This function is the exported wrapper that can be called from R. It has the
// inputs, ' and chunks G (the POIs) based on chunk_size to be able to stay
// within RAM requirements. ' @useDynLib FastVLA ' @importFrom Rcpp sourceCpp '
//@export
// [[Rcpp::export]]
int FastVLA_chunked_sota(
    const arma::mat &Y, SEXP Gptr, const arma::ivec &v_index,
    const arma::ivec &i_index, const arma::mat &X, const int &chunk_size,
    const std::string &dir, const std::vector<std::string> cnames,
    const std::vector<std::string> vnames, const std::string suffix,
    const double epss = 1e-6, const double &mafthresh = 0.005,
    const double pca_var_explained = 0.95, const int max_threads = 2,
    const int max_blas_threads = 1) {
  // Clean up previous run
  // delete_dir(dir);
  BLASLibraryManager blas_mgr;
  blas_mgr.detect_lib();
  int cur_blas_threads = blas_mgr.get_num_threads();
  if (cur_blas_threads != max_blas_threads) {
    blas_mgr.set_num_threads(max_blas_threads);
  }

  omp_set_dynamic(0);
  omp_set_num_threads(max_threads);
  int num_chunks = std::floor(v_index.n_elem / chunk_size);

  // get boolean inclusion vectors
  arma::mat xy_bools(Y.n_rows, Y.n_cols);
  uvec badX = find_nan(arma::sum(X, 1.0));
  // put in each column booleans
  for (uword i = 0; i < Y.n_cols; i++) {
    // get X NANs for removal
    arma::vec x_bools(Y.n_rows, fill::ones);
    x_bools.elem(badX).fill(0.0);
    // add Y NANs for removal
    arma::uvec badY = find_nan(Y.col(i));
    x_bools.elem(badY).fill(0.0);
    xy_bools.col(i) = x_bools;
  }
  // std::vector<std::string>::const_iterator first = vnames.begin();
  // std::vector<std::string>::const_iterator last = vnames.begin() +
  // chunk_size;
  std::vector<std::string> vnames_subset(vnames.begin(),
                                         vnames.begin() + chunk_size);
  // do first chunk
  arma::mat G = scanBEDMatrix(Gptr, i_index, v_index.rows(0, chunk_size - 1));
  int a =
      FastVLA_cpp_internal(Y, G, X, dir, xy_bools, cnames, vnames_subset,
                           suffix, epss, mafthresh, false, pca_var_explained);
  // do rest of chunks, appending output
  // #pragma omp parallel for
  for (int i = 1; i < num_chunks; i++) {
    G = scanBEDMatrix(Gptr, i_index,
                      v_index.rows(chunk_size * i, chunk_size * (i + 1) - 1));
    // first = vnames.begin() + i*chunk_size;
    // last = vnames.begin() + (i+1)*chunk_size - 1;
    vnames_subset.assign(vnames.begin() + i * chunk_size,
                         vnames.begin() + (i + 1) * chunk_size);
    a = std::min(a, FastVLA_cpp_internal(Y, G, X, dir, xy_bools, cnames,
                                         vnames_subset, suffix, epss, mafthresh,
                                         true, pca_var_explained));
  }
  // do last chunk if needed, appending output
  if ((int)v_index.n_elem > chunk_size * num_chunks) {
    arma::mat G = scanBEDMatrix(
        Gptr, i_index,
        v_index.rows(chunk_size * num_chunks, v_index.n_elem - 1));
    // first = vnames.begin() + num_chunks*chunk_size;
    // last = vnames.end();
    std::vector<std::string> vnames_subset2(
        vnames.begin() + num_chunks * chunk_size, vnames.end());
    a = std::min(a, FastVLA_cpp_internal(Y, G, X, dir, xy_bools, cnames,
                                         vnames_subset2, suffix, epss,
                                         mafthresh, true, pca_var_explained));
  }
  blas_mgr.set_num_threads(cur_blas_threads);
  return a;
}

/// @brief Gets the correct column of Y, X, G, and w (the mask)
/// @param i Y column
/// @param value G column
/// @param y Y placeholder to fill
/// @param g G col placeholder to fill
/// @param x X mat placeholder to fill
/// @param w Mask placeholder to fill
/// @param Y Y matrix
/// @param G G matrix
/// @param X X matrix
/// @param mask Mask of Y/X
void get_ygxw2(const uword &value, vec &y, vec &g, mat &x, vec &w, const mat &Y,
               const mat &G, const mat &X) {
  y = Y.col(0);
  x = X;
  g = G.col(value);

  w = g;
  w.fill(1.0);
  w(find_nan(g)).fill(0.0);
  // w = ws.col(value);

  uvec bad_inds = find(w == 0);

  g.elem(bad_inds).fill(0.0);
  y.elem(bad_inds).fill(0.0);
  x.rows(bad_inds).fill(0.0);
  return;
}

/// @brief Workhorse function of VLA regressions with a Single Y with no added
/// SVD
/// @param Y Phenotypes/outcomes matrix
/// @param G Matrix of POIs (subset by the calling function)
/// @param X Matrix of covariates
/// @param dir Directory for output
/// @param cnames Subfolder names of Y columns
/// @param vnames Rownames for output (POIs)
/// @param suffix Name to append for output in correspond dir/cnames folder
/// @param epss Epsilon threshold for SVD/inversion
/// @param mafthresh Threshold for Mean allele frequency
/// @param append Boolean -- create a new file or append to old?
/// @return 1 if successful
int FastVLA_cpp_internal3(const arma::mat &Y, const arma::mat &G,
                          const arma::mat &X, const std::string &dir,
                          const std::vector<std::string> cnames,
                          const std::vector<std::string> vnames,
                          const std::string suffix, const double epss = 1e-6,
                          const double &mafthresh = 0.01,
                          const bool &append = false) {

  // precreate objects
  double F;
  mat x = X; // copy of X

  vec w(Y.n_rows);

  // to store necessary additional vectors
  vec g_unscaled(Y.n_rows);
  vec yadj(Y.n_rows);
  vec xadj(Y.n_rows);
  vec gadj(Y.n_rows);
  vec hadj(Y.n_rows);
  mat Xadj(X.n_cols, Y.n_rows);
  vec yvQTL(Y.n_rows);

  mat temp_matmul(Y.n_rows, 3);
  mat Xg(Y.n_rows, 2);
  mat Vg(2, 2);

  vec est(2);

  vec z(Y.n_rows);

  double p, gmean, df, SST, SSE, MSE;

  double rankk;

  // for 1d
  double ans, estt;

  // constants for conversion from log pval to -log10 pval
  const double logg10 = -1.0 * log(10.0);
  const double logg2 = log(2.0);

  // output holder
  mat output(G.n_cols, 33);
  // condition checks placeholders
  double condition, caf;
  bool ones, twos;
  bool temp;
  double one_counter;
  vec ps(G.n_cols);
  vec whichcondition(G.n_cols);

  // create output directory
  fs::create_directory(dir);

  // X = X_orig;

  // fill output with NANs
  output.fill(datum::nan);

  // create file name for this Y
  std::stringstream ss;
  ss << dir << "/" << cnames[0];
  std::string result_file = ss.str();

  // make sure result directory is there then add the file inside it
  fs::create_directory(result_file);

  ss << "/" << suffix << ".tsv";
  result_file = ss.str();

  // loop over POIs/genes and check if good gene, bad gene, or minimal gene
  // (only 2 unique values)
  for (uword j = 0; j < G.n_cols; j++) {
    gadj = G.col(j);
    // condition = 10.0;
    w.fill(1.0);
    w(find_nan(gadj)).fill(0.0);

    uvec bad_inds = find(w == 0);
    gadj.elem(bad_inds).fill(0.0);

    // put it into w storage
    //  ws.col(j) = w;
    p = arma::accu(w);
    gmean = dot(gadj.t(), w) / p;
    output(j, 0) = p;
    caf = gmean / 2.0;
    // changed as of 10/10 to be CAF instead of MAF
    //    output(j,1) = std::min(caf, 1.0 - caf);
    output(j, 1) = caf;
    //   whichcondition[j] = check_condition(gadj, caf, mafthresh);
    condition =
        -1.0 * ((caf < mafthresh) | (caf > (1.0 - mafthresh))); // 1 is all NANs
    ps[j] = p;
    if ((caf < mafthresh) | (caf > (1.0 - mafthresh))) {
      whichcondition[j] = 2.0;
    } else {
      // check if ones and twos
      ones = false;
      twos = false;
      one_counter = 0.0;
      // 10/10 add a counter for number of ones;
      for (const auto &value : gadj) {
        temp = (value == 1.0);
        ones = ones | temp;
        one_counter += 1.0 * (temp);
        twos = twos | (value == 2.0);
      }
      whichcondition[j] = 2.0 - ones - twos;
      output(j, 2) = one_counter;
    }
  }

  // find conditions that need further calculations
  arma::uvec oneg = find(whichcondition == 1.0);
  arma::uvec good = find(whichcondition == 0.0);

  // do 1D G
#pragma omp parallel for private(x, w, g_unscaled, Xadj, temp_matmul, yadj,    \
                                 gadj, hadj, Xg, Vg, SST, SSE, MSE, df, z, F,  \
                                 est, ans, estt)
  for (const auto &value : oneg) {

    get_ygxw2(value, yadj, gadj, x, w, Y, G, X);
    try {
      Xadj = x;
      // rankk = arma::conv_to<double>::from(Xadj.n_cols);
      rankk = (double)Xadj.n_cols;
      // Xadj = manual_svd2(x, rankk, epss);
      output(value, 3) = rankk;
      output(value, 4) = 1.0;

      // create un-adjusted copy for vQTL
      g_unscaled = gadj;

      hadj = arma::square(gadj) % w;

      Xadj = Xadj * x.t();

      // sGWAS regress out X
      regress_out_X(yadj, gadj, hadj, temp_matmul, x, Xadj, w);

      // common between sGWAS, fGWAS
      df = ps[value] - rankk - 1.0;
      SST = arma::accu(arma::square(yadj));

      // sGWAS: regress out G
      ans = dot(gadj.t(), gadj);
      estt = dot(gadj.t(), yadj) / ans;
      z = arma::square((yadj - gadj * estt) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      output_1d(value, 5, df, estt, MSE, ans, output, logg10, logg2);

      // vQTL: regression on squared residuals (z) from sGWAS
      Xg = arma::join_rows(w, g_unscaled);
      Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);
      // decrease by 2 more for vQTL to account for vQTL betas [intercept + G]
      df -= 2.0;
      est = Vg * (Xg.t() * z);
      MSE = arma::accu(arma::square((z - Xg * est) % w)) / df;
      // se = sqrt(MSE * Vg(1,1));
      output_1d(value, 9, df, est[1], MSE, 1.0 / Vg(1, 1), output, logg10,
                logg2);

      // add back the 2 DF removed for vQTL
      df += 2.0;

      // Regression 1 of VLA
      // No g^2 so no need for matrix
      // estimate betas and subtract out for SSE
      g_unscaled = gadj;
      ans = dot(gadj.t(), gadj);
      estt = dot(gadj.t(), yadj) / ans;
      gadj *= estt;
      z = arma::square((yadj - gadj) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      output_1d(value, 13, df, estt, MSE, ans, output, logg10, logg2);

      // add in F statistic (importance of adding G to the regression)
      F = ((SST - SSE) / 1.0) / MSE;
      output(value, 20) = F;
      output(value, 21) = (R::pf(F, 1, df, 0, 1)) / logg10;

      // STEP 2: Regress squared residuals from above on X + G
      // regress out X in VLA
      hadj = (z - x * (Xadj * z)) % w;
      SST = arma::accu(arma::square(hadj));

      // update DF! Remove added X + POI estimate dfs
      df -= (rankk + 1.0);

      // regress on G and get betas
      estt = dot(g_unscaled.t(), hadj) / ans;
      g_unscaled *= estt;
      z = arma::square((hadj - g_unscaled) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      output_1d(value, 22, df, estt, MSE, ans, output, logg10, logg2);

      // add in F-stat (importance of G in VLA variance regression) and
      // SST,SSE
      F = ((SST - SSE) / 1.0) / MSE;
      output(value, 29) = SST;
      output(value, 30) = SSE;
      output(value, 31) = F;
      output(value, 32) = (R::pf(F, 1, df, 0, 1)) / logg10;

    } catch (...) { // catch-all for issues
      output(value, 2) = 100;
    }
  }

  // do 2D G (has 0, 1, and 2 values)
#pragma omp parallel for private(x, w, g_unscaled, Xadj, temp_matmul, yadj,    \
                                 gadj, hadj, Xg, Vg, SST, SSE, MSE, df, z, F,  \
                                 est, ans, estt)
  for (const auto &value : good) {
    get_ygxw2(value, yadj, gadj, x, w, Y, G, X);
    // try to compute everything with a catch-all for failure
    try {
      Xadj = x;
      // rankk = arma::conv_to<double>::from(Xadj.n_cols);
      rankk = (double)Xadj.n_cols;
      // Xadj = manual_svd2(x, rankk, epss);
      output(value, 3) = rankk;
      output(value, 4) = 2.0;

      // create a copy pre-regressing out X
      g_unscaled = gadj;

      // as of 10/3 no more scaling of G
      //  g -= gmeans[value];
      hadj = arma::square(gadj) % w;

      Xadj = Xadj * x.t();

      regress_out_X(yadj, gadj, hadj, temp_matmul, x, Xadj, w);

      df = ps[value] - rankk - 1.0;
      SST = arma::accu(arma::square(yadj));

      // sGWAS (standard)
      ans = dot(gadj.t(), gadj);
      estt = dot(gadj.t(), yadj) / ans;
      z = arma::square((yadj - gadj * estt) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      // output relevant requested statistics of the regression
      output_1d(value, 5, df, estt, MSE, ans, output, logg10, logg2);

      // vQTL is 8-11 (standard for comparison)
      Xg = arma::join_rows(w, g_unscaled);
      Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);
      // decrease by 2 more for vQTL to account for vQTL betas
      df -= 2.0;
      est = Vg * (Xg.t() * z);
      MSE = arma::accu(arma::square((z - Xg * est) % w)) / df;

      output_1d(value, 9, df, est[1], MSE, 1.0 / Vg(1, 1), output, logg10,
                logg2);

      // NEW VLA STEP 1: regress on G and G^2

      Xg = join_rows(gadj, hadj);
      Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);

      // add back 1 df from vQTL original
      df += 1.0;
      // get beta estimates
      est = Vg * (Xg.t() * yadj);
      z = arma::square((yadj - Xg * est) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;
      // print output
      output_2d(value, 13, df, est, MSE, Vg, output, logg10, logg2);
      // add in STEP 1 F value
      F = ((SST - SSE) / 2.0) / MSE;
      output(value, 20) = F;
      output(value, 21) = (R::pf(F, 2, df, 0, 1)) / logg10;

      // STEP 2: Regress X out of residuals^2 for VLA
      hadj = (z - x * (Xadj * z)) % w;

      // update DF NEW TO REMOVE X AND G AND G^2 AGAIN
      df -= (rankk + 2.0);

      // get SST with only X regressed
      SST = arma::accu(arma::square(hadj));
      est = Vg * (Xg.t() * hadj);
      z = arma::square((hadj - Xg * est) % w);

      // get SSE with G and G^2 regressed as well
      SSE = arma::accu(z);
      MSE = SSE / df;

      // output relevant statistics for G and G^2
      output_2d(value, 22, df, est, MSE, Vg, output, logg10, logg2);

      // add in F for Level 2 VLA regression (and SST,SSE)
      F = ((SST - SSE) / 2.0) / MSE;
      output(value, 29) = SST;
      output(value, 30) = SSE;
      output(value, 31) = F;
      output(value, 32) = (R::pf(F, 2, df, 0, 1)) / logg10;
    } catch (...) { // catch-all for issues
      output(value, 3) = 100;
    }
  }
  // output the necessary results to file from the matrix
  printer(result_file, output, append, vnames);

  return 1; // success!
}

/// @brief Workhorse function of VLA regressions with a Single Y
/// @param Y Phenotypes/outcomes matrix
/// @param G Matrix of POIs (subset by the calling function)
/// @param X Matrix of covariates
/// @param dir Directory for output
/// @param cnames Subfolder names of Y columns
/// @param vnames Rownames for output (POIs)
/// @param suffix Name to append for output in correspond dir/cnames folder
/// @param epss Epsilon threshold for SVD/inversion
/// @param mafthresh Threshold for Mean allele frequency
/// @param append Boolean -- create a new file or append to old?
/// @return 1 if successful
int FastVLA_cpp_internal2(const arma::mat &Y, const arma::mat &G,
                          const arma::mat &X, const std::string &dir,
                          const std::vector<std::string> cnames,
                          const std::vector<std::string> vnames,
                          const std::string suffix, const double epss = 1e-6,
                          const double &mafthresh = 0.01,
                          const bool &append = false) {

  // precreate objects
  double F;
  mat x = X; // copy of X

  vec w(Y.n_rows);

  // to store necessary additional vectors
  vec g_unscaled(Y.n_rows);
  vec yadj(Y.n_rows);
  vec xadj(Y.n_rows);
  vec gadj(Y.n_rows);
  vec hadj(Y.n_rows);
  mat Xadj(X.n_cols, Y.n_rows);
  vec yvQTL(Y.n_rows);

  mat temp_matmul(Y.n_rows, 3);
  mat Xg(Y.n_rows, 2);
  mat Vg(2, 2);

  vec est(2);

  vec z(Y.n_rows);

  double p, gmean, df, SST, SSE, MSE;

  double rankk;

  // for 1d
  double ans, estt;

  // constants for conversion from log pval to -log10 pval
  const double logg10 = -1.0 * log(10.0);
  const double logg2 = log(2.0);

  // output holder
  mat output(G.n_cols, 33);
  // condition checks placeholders
  double condition, caf;
  bool ones, twos;
  bool temp;
  double one_counter;
  vec ps(G.n_cols);
  vec whichcondition(G.n_cols);

  // create output directory
  fs::create_directory(dir);

  // X = X_orig;

  // fill output with NANs
  output.fill(datum::nan);

  // create file name for this Y
  std::stringstream ss;
  ss << dir << "/" << cnames[0];
  std::string result_file = ss.str();

  // make sure result directory is there then add the file inside it
  fs::create_directory(result_file);

  ss << "/" << suffix << ".tsv";
  result_file = ss.str();

  // loop over POIs/genes and check if good gene, bad gene, or minimal gene
  // (only 2 unique values)
  for (uword j = 0; j < G.n_cols; j++) {
    gadj = G.col(j);
    // condition = 10.0;
    w.fill(1.0);
    w(find_nan(gadj)).fill(0.0);

    uvec bad_inds = find(w == 0);
    gadj.elem(bad_inds).fill(0.0);

    // put it into w storage
    //  ws.col(j) = w;
    p = arma::accu(w);
    gmean = dot(gadj.t(), w) / p;
    output(j, 0) = p;
    caf = gmean / 2.0;
    // changed as of 10/10 to be CAF instead of MAF
    //    output(j,1) = std::min(caf, 1.0 - caf);
    output(j, 1) = caf;
    //   whichcondition[j] = check_condition(gadj, caf, mafthresh);
    condition =
        -1.0 * ((caf < mafthresh) | (caf > (1.0 - mafthresh))); // 1 is all NANs
    ps[j] = p;
    if ((caf < mafthresh) | (caf > (1.0 - mafthresh))) {
      whichcondition[j] = 2.0;
    } else {
      // check if ones and twos
      ones = false;
      twos = false;
      one_counter = 0.0;
      // 10/10 add a counter for number of ones;
      for (const auto &value : gadj) {
        temp = (value == 1.0);
        ones = ones | temp;
        one_counter += 1.0 * (temp);
        twos = twos | (value == 2.0);
      }
      whichcondition[j] = 2.0 - ones - twos;
      output(j, 2) = one_counter;
    }
  }

  // find conditions that need further calculations
  arma::uvec oneg = find(whichcondition == 1.0);
  arma::uvec good = find(whichcondition == 0.0);

  // do 1D G
#pragma omp parallel for private(x, w, g_unscaled, Xadj, temp_matmul, yadj,    \
                                 gadj, hadj, Xg, Vg, SST, SSE, MSE, df, z, F,  \
                                 est, ans, estt)
  for (const auto &value : oneg) {

    get_ygxw2(value, yadj, gadj, x, w, Y, G, X);
    try {
      Xadj = manual_svd2(x, rankk, epss);
      output(value, 3) = rankk;
      output(value, 4) = 1.0;

      // create un-adjusted copy for vQTL
      g_unscaled = gadj;

      hadj = arma::square(gadj) % w;

      Xadj = Xadj * x.t();

      // sGWAS regress out X
      regress_out_X(yadj, gadj, hadj, temp_matmul, x, Xadj, w);

      // common between sGWAS, fGWAS
      df = ps[value] - rankk - 1.0;
      SST = arma::accu(arma::square(yadj));

      // sGWAS: regress out G
      ans = dot(gadj.t(), gadj);
      estt = dot(gadj.t(), yadj) / ans;
      z = arma::square((yadj - gadj * estt) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      output_1d(value, 5, df, estt, MSE, ans, output, logg10, logg2);

      // vQTL: regression on squared residuals (z) from sGWAS
      Xg = arma::join_rows(w, g_unscaled);
      Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);
      // decrease by 2 more for vQTL to account for vQTL betas [intercept + G]
      df -= 2.0;
      est = Vg * (Xg.t() * z);
      MSE = arma::accu(arma::square((z - Xg * est) % w)) / df;
      // se = sqrt(MSE * Vg(1,1));
      output_1d(value, 9, df, est[1], MSE, 1.0 / Vg(1, 1), output, logg10,
                logg2);

      // add back the 2 DF removed for vQTL
      df += 2.0;

      // Regression 1 of VLA
      // No g^2 so no need for matrix
      // estimate betas and subtract out for SSE
      g_unscaled = gadj;
      ans = dot(gadj.t(), gadj);
      estt = dot(gadj.t(), yadj) / ans;
      gadj *= estt;
      z = arma::square((yadj - gadj) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      output_1d(value, 13, df, estt, MSE, ans, output, logg10, logg2);

      // add in F statistic (importance of adding G to the regression)
      F = ((SST - SSE) / 1.0) / MSE;
      output(value, 20) = F;
      output(value, 21) = (R::pf(F, 1, df, 0, 1)) / logg10;

      // STEP 2: Regress squared residuals from above on X + G
      // regress out X in VLA
      hadj = (z - x * (Xadj * z)) % w;
      SST = arma::accu(arma::square(hadj));

      // update DF! Remove added X + POI estimate dfs
      df -= (rankk + 1.0);

      // regress on G and get betas
      estt = dot(g_unscaled.t(), hadj) / ans;
      g_unscaled *= estt;
      z = arma::square((hadj - g_unscaled) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      output_1d(value, 22, df, estt, MSE, ans, output, logg10, logg2);

      // add in F-stat (importance of G in VLA variance regression) and
      // SST,SSE
      F = ((SST - SSE) / 1.0) / MSE;
      output(value, 29) = SST;
      output(value, 30) = SSE;
      output(value, 31) = F;
      output(value, 32) = (R::pf(F, 1, df, 0, 1)) / logg10;

    } catch (...) { // catch-all for issues
      output(value, 2) = 100;
    }
  }

  // do 2D G (has 0, 1, and 2 values)
#pragma omp parallel for private(x, w, g_unscaled, Xadj, temp_matmul, yadj,    \
                                 gadj, hadj, Xg, Vg, SST, SSE, MSE, df, z, F,  \
                                 est, ans, estt)
  for (const auto &value : good) {
    get_ygxw2(value, yadj, gadj, x, w, Y, G, X);
    // try to compute everything with a catch-all for failure
    try {
      Xadj = manual_svd2(x, rankk, epss);
      output(value, 3) = rankk;
      output(value, 4) = 2.0;

      // create a copy pre-regressing out X
      g_unscaled = gadj;

      // as of 10/3 no more scaling of G
      //  g -= gmeans[value];
      hadj = arma::square(gadj) % w;

      Xadj = Xadj * x.t();

      regress_out_X(yadj, gadj, hadj, temp_matmul, x, Xadj, w);

      df = ps[value] - rankk - 1.0;
      SST = arma::accu(arma::square(yadj));

      // sGWAS (standard)
      ans = dot(gadj.t(), gadj);
      estt = dot(gadj.t(), yadj) / ans;
      z = arma::square((yadj - gadj * estt) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;

      // output relevant requested statistics of the regression
      output_1d(value, 5, df, estt, MSE, ans, output, logg10, logg2);

      // vQTL is 8-11 (standard for comparison)
      Xg = arma::join_rows(w, g_unscaled);
      Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);
      // decrease by 2 more for vQTL to account for vQTL betas
      df -= 2.0;
      est = Vg * (Xg.t() * z);
      MSE = arma::accu(arma::square((z - Xg * est) % w)) / df;

      output_1d(value, 9, df, est[1], MSE, 1.0 / Vg(1, 1), output, logg10,
                logg2);

      // NEW VLA STEP 1: regress on G and G^2

      Xg = join_rows(gadj, hadj);
      Vg = arma::inv_sympd(Xg.t() * Xg, inv_opts::allow_approx);

      // add back 1 df from vQTL original
      df += 1.0;
      // get beta estimates
      est = Vg * (Xg.t() * yadj);
      z = arma::square((yadj - Xg * est) % w);
      SSE = arma::accu(z);
      MSE = SSE / df;
      // print output
      output_2d(value, 13, df, est, MSE, Vg, output, logg10, logg2);
      // add in STEP 1 F value
      F = ((SST - SSE) / 2.0) / MSE;
      output(value, 20) = F;
      output(value, 21) = (R::pf(F, 2, df, 0, 1)) / logg10;

      // STEP 2: Regress X out of residuals^2 for VLA
      hadj = (z - x * (Xadj * z)) % w;

      // update DF NEW TO REMOVE X AND G AND G^2 AGAIN
      df -= (rankk + 2.0);

      // get SST with only X regressed
      SST = arma::accu(arma::square(hadj));
      est = Vg * (Xg.t() * hadj);
      z = arma::square((hadj - Xg * est) % w);

      // get SSE with G and G^2 regressed as well
      SSE = arma::accu(z);
      MSE = SSE / df;

      // output relevant statistics for G and G^2
      output_2d(value, 22, df, est, MSE, Vg, output, logg10, logg2);

      // add in F for Level 2 VLA regression (and SST,SSE)
      F = ((SST - SSE) / 2.0) / MSE;
      output(value, 29) = SST;
      output(value, 30) = SSE;
      output(value, 31) = F;
      output(value, 32) = (R::pf(F, 2, df, 0, 1)) / logg10;
    } catch (...) { // catch-all for issues
      output(value, 3) = 100;
    }
  }
  // output the necessary results to file from the matrix
  printer(result_file, output, append, vnames);

  return 1; // success!
}

// needed to change logic for single Y
// This is the preferred function to match with Shail's logic for logistic
// regression. Note: Larger chunk_size is better if RAM permits. ' This function
// is the exported wrapper that can be called from R. It has the
//  inputs, ' and chunks G (the POIs) based on chunk_size to be able to stay
//  within RAM requirements. ' @useDynLib FastVLA ' @importFrom Rcpp sourceCpp '
//@param do_pca: Do PCA dimension reduction?
//@param add_intercept: Add an intercept to front of X?
//@export
//  [[Rcpp::export]]
int FastVLA_single_Y(const arma::vec &Y, SEXP Gptr, const arma::ivec &v_index,
                     const arma::ivec &i_index, const arma::mat &X,
                     const int &chunk_size, const std::string &dir,
                     const std::vector<std::string> cnames,
                     const std::vector<std::string> vnames,
                     const std::string suffix, const double epss = 1e-6,
                     const double &mafthresh = 0.005,
                     const double pca_var_explained = 0.95,
                     const bool do_pca = true,
                     const bool add_intercept = true) {
  // omp_set_num_threads(2);
  int num_chunks = std::floor(v_index.n_elem / chunk_size);

  // get boolean inclusion vector
  uvec badX = find_nan(arma::sum(X, 1.0));
  // put in each column booleans
  // get X NANs for removal
  arma::vec x_bools(Y.n_rows, fill::ones);
  x_bools.elem(badX).fill(0.0);
  // add Y NANs for removal
  arma::uvec badY = find_nan(Y);
  x_bools.elem(badY).fill(0.0);

  arma::mat cov = X;
  arma::ivec _index = i_index;
  arma::mat pheno(Y.n_rows, 1);
  pheno.col(0) = Y;
  arma::uvec nan_idx = arma::find(x_bools == 0.0);
  _index.shed_rows(nan_idx);
  pheno.shed_rows(nan_idx);
  cov.shed_rows(nan_idx);

  auto start = std::chrono::high_resolution_clock::now();
  if (do_pca) {
    auto start = std::chrono::high_resolution_clock::now();
    cov = standardize_and_derank_print(cov, dir, cnames[0], suffix,
                                       pca_var_explained, double(cov.n_rows));

    auto end = std::chrono::high_resolution_clock::now();
    double timing =
        (double)std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start)
            .count();
    Rcpp::Rcout << "PCA ananlysis timing: " << timing << "ms" << std::endl;
  }
  auto end_t = std::chrono::high_resolution_clock::now();

  double timing = (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                      end_t - start)
                      .count();

  Rcpp::Rcout << "PCA ananlysis timing: " << timing << "ms" << std::endl;
  if (add_intercept) {
    arma::mat no_interactions = arma::mat(cov.n_rows, 1, arma::fill::ones);
    cov = arma::join_rows(no_interactions, cov);
  }

  // std::vector<std::string>::const_iterator first = vnames.begin();
  // std::vector<std::string>::const_iterator last = vnames.begin() +
  // chunk_size;
  std::vector<std::string> vnames_subset(vnames.begin(),
                                         vnames.begin() + chunk_size);
  // do first chunk
  arma::mat G = scanBEDMatrix(Gptr, _index, v_index.rows(0, chunk_size - 1));
  int a = FastVLA_cpp_internal2(pheno, G, cov, dir, cnames, vnames_subset,
                                suffix, epss, mafthresh, false);
  // do rest of chunks, appending output
  // #pragma omp parallel for
  for (int i = 1; i < num_chunks; i++) {
    G = scanBEDMatrix(Gptr, _index,
                      v_index.rows(chunk_size * i, chunk_size * (i + 1) - 1));
    // first = vnames.begin() + i*chunk_size;
    // last = vnames.begin() + (i+1)*chunk_size - 1;
    vnames_subset.assign(vnames.begin() + i * chunk_size,
                         vnames.begin() + (i + 1) * chunk_size);
    a = std::min(a, FastVLA_cpp_internal2(pheno, G, cov, dir, cnames,
                                          vnames_subset, suffix, epss,
                                          mafthresh, true));
  }
  // do last chunk if needed, appending output
  if ((int)v_index.n_elem > chunk_size * num_chunks) {
    arma::mat G = scanBEDMatrix(
        Gptr, _index,
        v_index.rows(chunk_size * num_chunks, v_index.n_elem - 1));
    // first = vnames.begin() + num_chunks*chunk_size;
    // last = vnames.end();
    std::vector<std::string> vnames_subset2(
        vnames.begin() + num_chunks * chunk_size, vnames.end());
    a = std::min(a, FastVLA_cpp_internal2(pheno, G, cov, dir, cnames,
                                          vnames_subset2, suffix, epss,
                                          mafthresh, true));
  }
  return a;
}

// needed to change logic for single Y
// This is the preferred function to match with Shail's logic for logistic
// regression. Note: Larger chunk_size is better if RAM permits. ' This function
// is the exported wrapper that can be called from R. It has the
//  inputs, ' and chunks G (the POIs) based on chunk_size to be able to stay
//  within RAM requirements. ' @useDynLib FastVLA ' @importFrom Rcpp sourceCpp '
//@param do_pca: Do PCA dimension reduction?
//@export
//  [[Rcpp::export]]
int FastVLA_single_Y_fast(
    const arma::vec &Y, SEXP Gptr, const arma::ivec &v_index,
    const arma::ivec &i_index, const arma::mat &X, const int &chunk_size,
    const std::string &dir, const std::vector<std::string> cnames,
    const std::vector<std::string> vnames, const std::string suffix,
    const double epss = 1e-6, const double &mafthresh = 0.005,
    const double pca_var_explained = 0.95, const bool do_pca = true,
    const bool add_intercept = true) {
  omp_set_num_threads(2);
  int num_chunks = std::floor(v_index.n_elem / chunk_size);

  // get boolean inclusion vector
  uvec badX = find_nan(arma::sum(X, 1.0));
  // put in each column booleans
  // get X NANs for removal
  arma::vec x_bools(Y.n_rows, fill::ones);
  x_bools.elem(badX).fill(0.0);
  // add Y NANs for removal
  arma::uvec badY = find_nan(Y);
  x_bools.elem(badY).fill(0.0);

  arma::mat cov = X;
  arma::ivec _index = i_index;
  arma::mat pheno(Y.n_rows, 1);
  pheno.col(0) = Y;
  arma::uvec nan_idx = arma::find(x_bools == 0.0);
  _index.shed_rows(nan_idx);
  pheno.shed_rows(nan_idx);
  cov.shed_rows(nan_idx);

  auto start = std::chrono::high_resolution_clock::now();
  if (do_pca) {
    auto start = std::chrono::high_resolution_clock::now();
    cov = standardize_and_derank_print(cov, dir, cnames[0], suffix,
                                       pca_var_explained, double(cov.n_rows));

    auto end = std::chrono::high_resolution_clock::now();
    double timing =
        (double)std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start)
            .count();
    Rcpp::Rcout << "PCA ananlysis timing: " << timing << "ms" << std::endl;
  }
  if (add_intercept) {
    arma::mat no_interactions = arma::mat(cov.n_rows, 1, arma::fill::ones);
    cov = arma::join_rows(no_interactions, cov);
  }

  // std::vector<std::string>::const_iterator first = vnames.begin();
  // std::vector<std::string>::const_iterator last = vnames.begin() +
  // chunk_size;
  std::vector<std::string> vnames_subset(vnames.begin(),
                                         vnames.begin() + chunk_size);
  // do first chunk
  arma::mat G = scanBEDMatrix(Gptr, _index, v_index.rows(0, chunk_size - 1));
  int a = FastVLA_cpp_internal3(pheno, G, cov, dir, cnames, vnames_subset,
                                suffix, epss, mafthresh, false);
  // do rest of chunks, appending output
  // #pragma omp parallel for
  for (int i = 1; i < num_chunks; i++) {
    G = scanBEDMatrix(Gptr, _index,
                      v_index.rows(chunk_size * i, chunk_size * (i + 1) - 1));
    // first = vnames.begin() + i*chunk_size;
    // last = vnames.begin() + (i+1)*chunk_size - 1;
    vnames_subset.assign(vnames.begin() + i * chunk_size,
                         vnames.begin() + (i + 1) * chunk_size);
    a = std::min(a, FastVLA_cpp_internal3(pheno, G, cov, dir, cnames,
                                          vnames_subset, suffix, epss,
                                          mafthresh, true));
  }
  // do last chunk if needed, appending output
  if ((int)v_index.n_elem > chunk_size * num_chunks) {
    arma::mat G = scanBEDMatrix(
        Gptr, _index,
        v_index.rows(chunk_size * num_chunks, v_index.n_elem - 1));
    // first = vnames.begin() + num_chunks*chunk_size;
    // last = vnames.end();
    std::vector<std::string> vnames_subset2(
        vnames.begin() + num_chunks * chunk_size, vnames.end());
    a = std::min(a, FastVLA_cpp_internal3(pheno, G, cov, dir, cnames,
                                          vnames_subset2, suffix, epss,
                                          mafthresh, true));
  }
  return a;
}
