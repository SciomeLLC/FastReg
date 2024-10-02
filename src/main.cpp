#define ARMA_WARN_LEVEL 0

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <chunker.h>
#include <config.h>
#include <covariate.h>
#include <covariate_matrix.h>
#include <ctime>
#include <fr_matrix.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <pheno_matrix.h>
#include <regression.h>
#include <sstream>
#include <stdio.h>
#include <strata.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "BEDMatrix.h"
#include <fr_result.h>
#include <poi_matrix.h>
#include <reader.h>
#include <utils.h>

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

#include <RcppEigen.h>

struct ProcResult {
  int timing_results[4] = {0, 0, 0, 0};
  double process_nonconvergence_status = 0.0;
  double process_total_filtered_pois = 0.0;
  std::mutex mtx;

  void accumulate(ProcResult &proc_res) {
    std::lock_guard<std::mutex> lock(mtx);
    for (int i = 0; i < 4; ++i) {
      timing_results[i] += proc_res.timing_results[i];
    }
    process_nonconvergence_status += proc_res.process_nonconvergence_status;
    process_total_filtered_pois += proc_res.process_total_filtered_pois;
  }

  void print_convergence_percentage(double nonconvergence_status,
                                    double filtered_pois) {
    double noncovergence_percent =
        (nonconvergence_status / filtered_pois) * 100;
    if (noncovergence_percent > 0.0) {
      Rcpp::Rcout << nonconvergence_status << " out of " << filtered_pois
                  << " (" << std::setprecision(2) << std::fixed
                  << noncovergence_percent
                  << "%) POIs did not meet relative and absolute convergence "
                     "threshold."
                  << std::endl;
    }
  }
  void print_nonconvergence_summary() {
    double noncovergence_percent =
        (process_nonconvergence_status / process_total_filtered_pois) * 100;
    Rcpp::Rcout << process_nonconvergence_status << " out of "
                << process_total_filtered_pois << " (" << std::setprecision(2)
                << std::fixed << noncovergence_percent
                << "%) POIs did not meet relative and absolute convergence "
                   "threshold."
                << std::endl;
  }
  void print_timing_summary(int process_id) {
    auto end =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    Rcpp::Rcout << "Timing Summary for process: " << process_id + 1
                << std::endl;
    Rcpp::Rcout << "Reading HDF5: " << timing_results[0] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Writing results: " << timing_results[1] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Memory allocation: " << timing_results[2] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Regression: " << timing_results[3] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Completed process " << process_id + 1
                << " at: " << std::ctime(&end) << std::endl;
  }

  void print_totals_summary(double concatenation_time, double compression_time,
                            std::string regression_type, size_t individuals,
                            int num_pois, int num_threads) {
    Rcpp::Rcout << "-----------------------------------------" << std::endl;
    Rcpp::Rcout << "Timing Summary: " << std::endl;
    Rcpp::Rcout << "Reading HDF5: " << timing_results[0] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Writing results: " << timing_results[1] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Memory allocation: " << timing_results[2] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Regression: " << timing_results[3] / 1000.0 << "s"
                << std::endl;
    Rcpp::Rcout << "Results concatenation: " << concatenation_time / 1000.0
                << "s" << std::endl;

    auto end =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    Rcpp::Rcout << "Completed " << regression_type << " regression -"
                << std::endl;
    Rcpp::Rcout << "\t\tnum individuals: " << individuals << std::endl;
    Rcpp::Rcout << "\t\t~num pois: " << num_pois << std::endl;
    Rcpp::Rcout << "\t\twith openmp thread(s): " << num_threads << std::endl;
    Rcpp::Rcout << "at: " << std::ctime(&end) << std::endl;
    Rcpp::Rcout << "-----------------------------------------" << std::endl;
  }
};

void process_chunk(int process_id, Config &config, FRMatrix &pheno_df,
                   FRMatrix &covar_df, std::string poi_file_path,
                   int chunk_size, int num_threads, bool use_blas,
                   ProcResult &proc_res) {
  // Load POI file
  POI poi(poi_file_path);
  poi.open(true);
  poi.get_values_dataset_id();
  poi.get_data_type();
  poi.get_names();
  poi.get_individuals();

  // Find common individuals
  std::vector<std::string> poi_names = poi.names;
  std::vector<std::string> common_ind =
      intersect_row_names(pheno_df.sort_map(true), covar_df.sort_map(true));
  std::vector<std::string> intersected_ind =
      intersect_row_names(common_ind, poi.individuals);
  if (intersected_ind.empty()) {
    stop("No overlapping individuals found in POI, pheno, and covar files");
  }

  // if (!config.POI_subset_file.empty()) {

  //     FRMatrix poi_subset(config.POI_subset_file,
  //     config.POI_subset_file_delim, config.POI_subset_rowname_col);
  //     std::vector<std::string> poi_names =
  //     intersect_row_names(poi_subset.str_data[0], poi.names);
  // }

  int num_poi = poi_names.size();
  if (num_poi == 0) {
    stop("No overlapping individuals found in POI, pheno, covar files");
  }
  // Rcpp::Rcout << "Stratifying data" << std::endl;
  // Stratify data
  Strata stratums;
  stratums.stratify(config.split_by, covar_df, intersected_ind);
  double memory_allocation_time = 0.0;
  double file_writing_time = 0.0;
  double poi_reading_time = 0.0;
  double regression_time = 0.0;
  for (int stratum = 0; stratum < stratums.nstrata; ++stratum) {
    std::string outfile_suffix = stratums.ids[stratum];
    if (!config.split_by[0].empty()) {
      Rcpp::Rcout << "Processing stratum: " << outfile_suffix.substr(1)
                  << std::endl;
    }
    std::vector<std::string> ind_set = stratums.index_list[outfile_suffix];
    int ct = 0;
    std::vector<int> ind_set_idx(ind_set.size());
    for (std::string ind : ind_set) {
      ind_set_idx[ct] = pheno_df.get_row_idx(ind);
      ct++;
    }

    // Rcpp::Rcout << "Init matrices" << std::endl;
    // initialize matrices
    FRMatrix pheno_matrix = pheno_df; // check col names

    // n individuals x # covariates
    FRMatrix covar_matrix = covar_df;
    FRMatrix covar_poi_interaction_matrix;
    covar_poi_interaction_matrix.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    covar_poi_interaction_matrix.row_names = covar_matrix.row_names;
    covar_poi_interaction_matrix.col_names = {{"poi", 0}};
    // n individuals x 1 or 1 + num interacting poi covars
    create_Z_matrix(covar_matrix, config.POI_covar_interactions,
                    covar_poi_interaction_matrix, !config.no_intercept);
    std::vector<int> nan_idx;
    std::vector<std::string> ind_set_filtered;

    // Rcpp::Rcout << "Identify missing" << std::endl;
    // identify missing values for covar, pheno matrix
    for (size_t i = 0; i < covar_matrix.data.n_rows; i++) {
      arma::uvec covar_nan_idx = arma::find_nonfinite(covar_matrix.data.row(i));
      arma::uvec pheno_nan_idx = arma::find_nonfinite(pheno_matrix.data.row(i));
      if (covar_nan_idx.size() > 0 || pheno_nan_idx.size() > 0) {
        nan_idx.push_back(i);
      } else {
        ind_set_filtered.push_back(ind_set[i]);
      }
    }

    // Rcpp::Rcout << "Removing missing" << std::endl;
    // remove from covar, pheno
    covar_matrix.data.shed_rows(arma::conv_to<arma::uvec>::from(nan_idx));
    pheno_matrix.data.shed_rows(arma::conv_to<arma::uvec>::from(nan_idx));
    covar_poi_interaction_matrix.data.shed_rows(
        arma::conv_to<arma::uvec>::from(nan_idx));
    covar_matrix.row_names =
        std::unordered_map<std::string, int>(ind_set_filtered.size());
    pheno_matrix.row_names =
        std::unordered_map<std::string, int>(ind_set_filtered.size());
    for (size_t j = 0; j < covar_matrix.data.n_rows; j++) {
      covar_matrix.row_names[ind_set_filtered[j]] = j;
      pheno_matrix.row_names[ind_set_filtered[j]] = j;
      covar_poi_interaction_matrix.row_names[ind_set_filtered[j]] = j;
    }
    std::vector<std::string> strat_individuals(ind_set_filtered.size());
    std::transform(
        ind_set_filtered.begin(), ind_set_filtered.end(),
        strat_individuals.begin(), [&intersected_ind](const std::string &elem) {
          return intersected_ind[std::distance(
              intersected_ind.begin(),
              std::find(intersected_ind.begin(), intersected_ind.end(), elem))];
        });
    double nonconvergence_status = 0.0;
    double filtered_pois = 0.0;

    int num_parallel_poi_blocks =
        (int)std::ceil((double)num_poi / (double)chunk_size);

    FRMatrix poi_matrix;
    auto start_time = std::chrono::high_resolution_clock::now();
    // allocate memory space for H5 file to read into
    for (int block = 0; block < num_parallel_poi_blocks; block++) {
      int start_chunk = block * chunk_size;
      int end_chunk = start_chunk + chunk_size;
      if (end_chunk >= num_poi) {
        end_chunk = num_poi;
      }

      // Rcpp::Rcout << "Reading poi chunk" << std::endl;
      std::vector<std::string> poi_names_chunk(poi_names.begin() + start_chunk,
                                               poi_names.begin() + end_chunk);
      // Rcpp::Rcout << "start_chunk: " << start_chunk << "\nend_chunk: " <<
      // end_chunk << "\nblock: " << block << "/" << num_parallel_poi_blocks <<
      // std::endl;
      FRMatrix poi_matrix = poi.read_chunk(poi.individuals, poi_names_chunk);
      // Rcpp::Rcout << "loaded chunk" << std::endl;
      std::vector<std::string> srt_cols_2 = poi_matrix.sort_map(false);
      std::vector<std::string> drop_rows =
          set_diff(poi.individuals, strat_individuals);

      int num_dropped = poi.individuals.size() - strat_individuals.size();
      arma::uvec drop_row_idx(drop_rows.size());
      for (size_t i = 0; i < drop_rows.size(); i++) {
        drop_row_idx[i] = poi_matrix.row_names[drop_rows[i]];
      }

      std::unordered_map<std::string, int> new_row_names(
          strat_individuals.size());
      for (auto &ind : strat_individuals) {
        new_row_names[ind] = poi_matrix.row_names[ind] - num_dropped;
      }
      poi_matrix.data.shed_rows(drop_row_idx);
      poi_matrix.row_names = new_row_names;

      srt_cols_2 = poi_matrix.sort_map(false);
      auto end_time = std::chrono::high_resolution_clock::now();
      poi_reading_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();

      if (config.POI_type == "genotype") {
        FRMatrix filtered =
            filter_poi(poi_matrix, config.maf_threshold, config.hwe_threshold);
        arma::uvec filtered_col = arma::find(filtered.data.row(5) == 0);

        if (filtered.data.n_cols == 0 ||
            filtered_col.n_elem == poi_matrix.data.n_cols) {
          Rcpp::Rcout << "no POI passed filtering" << std::endl;
          return;
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
        srt_cols_2 = poi_matrix.sort_map(false);
        transform_poi(poi_matrix, config.POI_effect_type);
        start_time = std::chrono::high_resolution_clock::now();
        std::string summary_name = "POI_Summary";
        filtered.write_summary(config.output_dir, summary_name, stratum,
                               process_id);
        end_time = std::chrono::high_resolution_clock::now();
        file_writing_time +=
            (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time)
                .count();
      }

      filtered_pois += poi_matrix.data.n_cols;
      // Rcpp::Rcout << "filtered pois" << std::endl;
      start_time = std::chrono::high_resolution_clock::now();
      FRMatrix beta_est;
      FRMatrix se_beta;
      int num_parms =
          covar_poi_interaction_matrix.data.n_cols + covar_matrix.data.n_cols;
      beta_est.data =
          arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
      arma::fcolvec beta_rel_errs =
          arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);
      arma::fcolvec beta_abs_errs =
          arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);
      arma::fcolvec iters =
          arma::fcolvec(poi_matrix.data.n_cols, arma::fill::zeros);

      se_beta.data =
          arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);

      FRMatrix neglog10_pvl;
      neglog10_pvl.data =
          arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);

      for (auto &col_name : covar_matrix.col_names) {
        beta_est.row_names[col_name.first] = col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] =
            beta_est.row_names[col_name.first];
      }
      for (auto &col_name : covar_poi_interaction_matrix.col_names) {
        beta_est.row_names[col_name.first] =
            covar_matrix.col_names.size() + col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] =
            beta_est.row_names[col_name.first];
      }

      beta_est.col_names = covar_matrix.row_names;
      se_beta.col_names = beta_est.col_names;
      neglog10_pvl.col_names = beta_est.col_names;
      std::vector<std::string> srt_cols = poi_matrix.sort_map(false);

      // Rcpp::Rcout << "Creating W2" << std::endl;
      arma::umat W2 = arma::umat(poi_matrix.data.n_rows, poi_matrix.data.n_cols,
                                 arma::fill::ones);
      for (arma::uword v = 0; v < poi_matrix.data.n_cols; v++) {
        arma::uvec G_na = arma::find_nonfinite(poi_matrix.data.col(v));
        for (arma::uword i = 0; i < G_na.n_elem; i++) {
          W2(G_na(i), v) = 0;
          poi_matrix.data(G_na(i), v) = 0;
        }
      }

      end_time = std::chrono::high_resolution_clock::now();
      memory_allocation_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();

#if !defined(__APPLE__) && !defined(__MACH__)
      omp_set_num_threads(num_threads);
#endif
      start_time = std::chrono::high_resolution_clock::now();
      std::unique_ptr<RegressionBase> regression;
      if (config.regression_type == "logistic") {
        regression.reset(new LogisticRegression());
      } else {
        regression.reset(new LinearRegression());
      }

      regression->run(
          covar_matrix, pheno_matrix, poi_matrix, covar_poi_interaction_matrix,
          W2, beta_est, se_beta, neglog10_pvl, beta_rel_errs, beta_abs_errs,
          iters, config.max_iter, config.p_value_type == "t.dist", use_blas);
      end_time = std::chrono::high_resolution_clock::now();
      regression_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();
      start_time = std::chrono::high_resolution_clock::now();
      FRMatrix::write_results(beta_est, se_beta, neglog10_pvl, W2,
                              beta_rel_errs, beta_abs_errs, iters, srt_cols,
                              config.output_dir, "Results", stratum,
                              config.output_exclude_covar, process_id + 1);
      end_time = std::chrono::high_resolution_clock::now();
      file_writing_time +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                start_time)
              .count();
      poi_matrix.col_names.clear();

      arma::fcolvec convergence = arma::conv_to<fcolvec>::from(
          (beta_rel_errs > config.rel_conv_tolerance) &&
          (beta_abs_errs > config.abs_conv_tolerance));
      nonconvergence_status = arma::sum(convergence);
      // proc_res.print_convergence_percentage(nonconvergence_status,
      // filtered_pois);
      proc_res.process_nonconvergence_status += nonconvergence_status;
      proc_res.process_total_filtered_pois += filtered_pois;
    }
  }

  poi.close_all();

  auto end =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  proc_res.timing_results[0] = poi_reading_time;
  proc_res.timing_results[1] = file_writing_time;
  proc_res.timing_results[2] = memory_allocation_time;
  proc_res.timing_results[3] = regression_time;
  proc_res.print_timing_summary(process_id);
}

bool is_BLAS_faster() {
  bool use_blas = true;
  arma::fmat tempM(250, 10000);
  tempM.randu(250, 10000);
  Eigen::MatrixXf mddata =
      Eigen::Map<Eigen::MatrixXf>(tempM.memptr(), tempM.n_rows, tempM.n_cols);
  ///////////////
  auto start = std::chrono::system_clock::now();
  tempM = tempM.t() * tempM;
  auto end = std::chrono::system_clock::now();
  ///////////////

  ///////////////
  auto start_a = std::chrono::system_clock::now();
  mddata = mddata.transpose() * mddata;
  auto end_a = std::chrono::system_clock::now();
  ///////////////
  const auto ms_int_a =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  const auto ms_int_b =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end_a - start_a);
  double BLAS = std::chrono::duration<double>(ms_int_a).count();
  double EIGEN = std::chrono::duration<double>(ms_int_b).count();

  if (BLAS / EIGEN < 1) {
    use_blas = true;
    Rcpp::Rcout << "Using loaded BLAS ilbrary." << std::endl;
  } else {
    use_blas = false;
    Rcpp::Rcout
        << "Using RcppEigen as it appears faster than loaded BLAS library. "
        << std::endl;
  }
  return use_blas;
}
// [[Rcpp::export]]
void FastRegCpp(
    const std::string phenotype, const std::string regression_type,
    const std::string pvalue_dist, bool output_exclude_covar,
    double maf_threshold, double hwe_threshold, bool no_intercept,
    double colinearity_rsq, int poi_block_size, int max_iter,
    double rel_conv_tolerance, double abs_conv_tolderance,
    int max_openmp_threads, const std::string pheno_file,
    const std::string pheno_rowname_cols, const std::string pheno_file_delim,
    const std::string covar_file, const std::string covar_rowname_cols,
    const std::string covar_file_delim, const std::string poi_file_dir,
    const std::string poi_file_delim, const std::string poi_file_format,
    const std::string poi_type, const std::string poi_effect_type,
    const Rcpp::StringVector covariates,
    const Rcpp::StringVector covariate_type,
    const Rcpp::LogicalVector covariate_standardize,
    const Rcpp::StringVector covariate_levels,
    const Rcpp::StringVector covariate_ref_level,
    const Rcpp::StringVector POI_covar_interactions_str,
    const Rcpp::StringVector split_by_str, const std::string output_dir,
    bool compress_results, int max_workers) {
  Config config(
      phenotype, regression_type, pvalue_dist, output_exclude_covar,
      maf_threshold, hwe_threshold, no_intercept, colinearity_rsq,
      poi_block_size, max_iter, rel_conv_tolerance, abs_conv_tolderance,
      max_openmp_threads, pheno_file, pheno_rowname_cols, pheno_file_delim,
      covar_file, covar_rowname_cols, covar_file_delim, poi_file_dir,
      poi_file_delim, poi_file_format, poi_type, poi_effect_type, covariates,
      covariate_type, covariate_standardize, covariate_levels,
      covariate_ref_level, POI_covar_interactions_str, split_by_str, output_dir,
      compress_results, max_workers);

  config.print();
  // Check for BLAS speed
  bool use_blas = is_BLAS_faster();
  // Clean up previous run
  if (dir_exists(config.output_dir)) {
    delete_dir(config.output_dir);
  }

  // FRMatrix pheno_df(config.pheno_file, config.pheno_file_delim,
  //                   config.pheno_rowname_cols, config.phenotype);
  PhenoMatrix pheno_matrix(config.pheno_file, config.pheno_file_delim,
                           config.pheno_rowname_cols, config.phenotype);
  FRMatrix pheno_df = pheno_matrix.create_matrix();

  CovariateMatrix cov_mat = CovariateMatrix(
      config.covar_file, config.covar_file_delim, config.covar_rowname_cols,
      config.covs, config.colinearity_rsq, config.no_intercept);

  FRMatrix covar_df = cov_mat.create_design_matrix();

  // Load the first POI file to calculate chunks

  std::string poi_file_path = config.poi_files[0];
  // Rcpp::Rcout << "POI file path: " << poi_file_path << std::endl;
  POI poi(poi_file_path);
  poi.open(true);
  poi.get_values_dataset_id();
  poi.get_names();
  poi.get_individuals();

  // Find common individuals
  std::vector<std::string> poi_names = poi.names;
  std::vector<std::string> common_ind =
      intersect_row_names(pheno_df.sort_map(true), covar_df.sort_map(true));
  std::vector<std::string> intersected_ind =
      intersect_row_names(common_ind, poi.individuals);

  Rcout << intersected_ind.size()
        << " common unique subjects in pheno.file, "
           "covar.file, and POI.file"
        << std::endl;
  if (intersected_ind.empty()) {
    stop("No overlapping individuals found in POI, pheno, and covar files");
  }

  int num_poi = poi_names.size();
  if (num_poi == 0) {
    stop("No overlapping individuals found in POI, pheno, covar files");
  }

  // Stratify data
  int num_poi_files = config.poi_files.size();
  Strata stratums;
  stratums.stratify(config.split_by, covar_df, intersected_ind);
  Chunker chunker =
      Chunker(num_poi, intersected_ind.size(), config.max_openmp_threads,
              config.poi_block_size, num_poi_files, config.max_workers);

  // clean up memory
  poi.close_all();

  // setup parallel processing
  // total_num_chunks
  int parallel_chunk_size = chunker.get_chunk_size();
  int num_threads = chunker.get_openmp_threads();
  const int timing_results_size = 4;
  double concatenation_time, compression_time = 0.0;
  ProcResult total_proc_res;
  int total_timing_results[timing_results_size] = {0, 0, 0, 0};
#if !defined(__APPLE__) && !defined(__MACH__)
  omp_set_num_threads(num_threads);
#endif
#ifdef _WIN32
  for (int i = 0; i < num_poi_files; i++) {
    ProcResult proc_res;
    int timing_results[] = {0, 0, 0, 0};
    process_chunk(i, config, pheno_df, covar_df, config.poi_files[i],
                  parallel_chunk_size, num_threads, use_blas, proc_res);

    total_proc_res.accumulate(proc_res);
  }
#else

  int num_processes_total = chunker.get_total_workers();
  int max_processes = chunker.get_num_workers();
  std::vector<int> pipe_file_descriptors(num_processes_total * 2);
  std::vector<pid_t> process_ids(num_processes_total);
  std::vector<bool> has_completed(num_processes_total, false);

  int num_processes_started = 0;
  int num_processes_completed = 0;

  while (num_processes_completed < num_processes_total) {
    while ((num_processes_started - num_processes_completed) < max_processes) {
      checkInterrupt();
      if (num_processes_started == num_processes_total) {
        break;
      }
      int i = num_processes_started;
      if (pipe(&pipe_file_descriptors[i * 2]) == -1) {
        perror("pipe");
        return;
      }
      process_ids[i] = fork();
      if (process_ids[i] == -1) {
        perror("fork");
        return;
      }
      std::string poi_file_path = config.poi_files[i];
      if (process_ids[i] == 0) {             // child process
        close(pipe_file_descriptors[i * 2]); // close read pipe

        ProcResult proc_res;
        process_chunk(i, config, pheno_df, covar_df, poi_file_path,
                      parallel_chunk_size, num_threads, use_blas, proc_res);
        ssize_t res = write(pipe_file_descriptors[i * 2 + 1], &proc_res,
                            sizeof(proc_res));
        close(pipe_file_descriptors[i * 2 + 1]);
        _exit(EXIT_SUCCESS);
        return;
      } else {
        close(pipe_file_descriptors[i * 2 + 1]);
        num_processes_started++;
      }
    }
    // Check for finished processes
    for (int i = 0; i < num_processes_started; i++) {
      checkInterrupt();
      if (process_ids[i] != 0) { // parent process
        if (!has_completed[i] && waitpid(process_ids[i], NULL, WNOHANG) > 0) {
          has_completed[i] = true;

          ProcResult proc_res;
          ssize_t res =
              read(pipe_file_descriptors[i * 2], &proc_res, sizeof(proc_res));
          close(pipe_file_descriptors[i * 2]);
          total_proc_res.accumulate(proc_res);
          num_processes_completed++;
        }
      }
    }
  }
#endif
  auto start_time = std::chrono::high_resolution_clock::now();
  FRMatrix::concatenate_results(config.output_dir, "Results", "Full");
  FRMatrix::concatenate_results(config.output_dir, "Convergence", "Full");
  if (config.POI_type == "genotype") {
    FRMatrix::concatenate_results(config.output_dir, "POI_Summary", "Full");
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  concatenation_time =
      (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                    start_time)
          .count();

  if (config.compress_results) {
    start_time = std::chrono::high_resolution_clock::now();
    FRMatrix::zip_results(config.output_dir);
    end_time = std::chrono::high_resolution_clock::now();
    compression_time +=
        (double)std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time)
            .count();
    Rcpp::Rcout << "Results compression: " << compression_time / 1000.0 << "s"
                << std::endl;
  }
  total_proc_res.print_nonconvergence_summary();
  total_proc_res.print_totals_summary(
      concatenation_time, compression_time, config.regression_type,
      intersected_ind.size(), num_poi * num_poi_files, num_threads);
  // Rcpp::Rcout << "-----------------------------------------" << std::endl;
  // Rcpp::Rcout << "Convergence Summary: " << std::endl;
  // Rcpp::Rcout << "Reading HDF5: " << total_timing_results[0] / 1000.0 << "s";
  // Rcpp::Rcout << "-----------------------------------------" << std::endl;
}

Rcpp::DataFrame arma_2_df(const arma::fmat &mat,
                          std::vector<std::string> row_names,
                          std::vector<std::string> col_names) {
  int n_rows = mat.n_rows;
  int n_cols = mat.n_cols;

  if (col_names.size() != static_cast<size_t>(n_cols)) {
    Rcpp::stop("Number of column names does not match number of columns in "
               "matrix. n_cols: %s, col_names: %s",
               n_cols, col_names.size());
  }

  if (row_names.size() != static_cast<size_t>(n_rows)) {
    Rcpp::stop("Number of row names does not match number of rows in matrix.");
  }

  for (size_t i = 0; i < col_names.size(); ++i) {
    if (col_names[i].empty()) {
      Rcpp::stop("Column name at position %d is empty.", i);
    }
  }

  for (size_t i = 0; i < row_names.size(); ++i) {
    if (row_names[i].empty()) {
      Rcpp::stop("Row name at position %d is empty.", i);
    }
  }
  Rcpp::List df_cols;
  for (int i = 0; i < n_cols; i++) {
    df_cols[col_names[i]] =
        Rcpp::NumericVector(mat.colptr(i), mat.colptr(i) + n_rows);
  }

  Rcpp::DataFrame df(df_cols);
  df.attr("row.names") = Rcpp::wrap(row_names);

  return df;
}

void process_chunk_vla(int process_id, Config &config, FRMatrix &pheno_df,
                       FRMatrix &covar_df, std::string poi_file_path,
                       int chunk_size, int num_threads, bool use_blas,
                       ProcResult &proc_res) {
  // Load POI file
  BEDReader bed_reader(poi_file_path);

  // Find common individuals
  std::vector<std::string> poi_names = bed_reader.get_names();
  std::vector<std::string> poi_individuals = bed_reader.get_individuals();
  std::vector<std::string> common_ind =
      intersect_row_names(pheno_df.sort_map(true), covar_df.sort_map(true));
  std::vector<std::string> intersected_ind =
      intersect_row_names(common_ind, poi_individuals);
  if (intersected_ind.empty()) {
    stop("No overlapping individuals found in POI, pheno, and covar files");
  }

  int num_poi = poi_names.size();
  if (num_poi == 0) {
    stop("No overlapping individuals found in POI, pheno, covar files");
  }

  // Rcpp::Rcout << "Stratifying data" << std::endl;
  // Stratify data
  Strata stratums;
  stratums.stratify(config.split_by, covar_df, intersected_ind);
  double memory_allocation_time = 0.0;
  double file_writing_time = 0.0;
  double poi_reading_time = 0.0;
  double regression_time = 0.0;
  for (int stratum = 0; stratum < stratums.nstrata; ++stratum) {
    std::string outfile_suffix = stratums.ids[stratum];
    if (!config.split_by[0].empty()) {
      Rcpp::Rcout << "Processing stratum: " << outfile_suffix.substr(1)
                  << std::endl;
    }
    std::vector<std::string> ind_set = stratums.index_list[outfile_suffix];
    int ct = 0;
    std::vector<int> ind_set_idx(ind_set.size());
    for (std::string ind : ind_set) {
      if (ind.empty()) {
        Rcpp::Rcout << "Found empty ind value: " << ind << std::endl;
        continue;
      }
      int idx = pheno_df.get_row_idx(ind);
      if (idx == -1) {
        Rcpp::Rcout << "Ind value not found in pheno_df: " << ind << std::endl;
        continue;
      }
      ind_set_idx[ct] = idx;
      ct++;
    }

    // Rcpp::Rcout << "Init matrices" << std::endl;
    // initialize matrices
    FRMatrix pheno_matrix = pheno_df; // check col names
    // n individuals x # covariates
    FRMatrix covar_matrix = covar_df;

    FRMatrix covar_poi_interaction_matrix;
    covar_poi_interaction_matrix.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    covar_poi_interaction_matrix.row_names = covar_matrix.row_names;
    covar_poi_interaction_matrix.col_names = {{"poi", 0}};
    covar_poi_interaction_matrix.col_names_arr.push_back("poi");

    FRMatrix no_interaction_matrix;
    no_interaction_matrix.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    no_interaction_matrix.row_names = covar_matrix.row_names;
    no_interaction_matrix.col_names = {{"poi", 0}};
    no_interaction_matrix.col_names_arr.push_back("poi");

    // n individuals x 1 or 1 + num interacting poi covars
    create_Z_matrix(covar_matrix, config.POI_covar_interactions,
                    covar_poi_interaction_matrix, !config.no_intercept);

    std::vector<int> nan_idx;
    std::vector<std::string> ind_set_filtered;

    for (size_t i = 0; i < covar_matrix.data.n_rows; i++) {
      auto idx = std::find(ind_set.begin(), ind_set.end(),
                           covar_matrix.row_names_arr[i]);
      if (idx == ind_set.end()) {
        nan_idx.push_back(i); // Missing individual
        continue;
      }
      arma::uvec covar_nan_idx = arma::find_nonfinite(covar_matrix.data.row(i));
      arma::uvec pheno_nan_idx = arma::find_nonfinite(pheno_matrix.data.row(i));
      if (covar_nan_idx.size() > 0 || pheno_nan_idx.size() > 0) {
        nan_idx.push_back(i);
      } else {
        if (ind_set[i].empty()) {
          Rcpp::Rcout << "Found empty ind string in missing: " << ind_set[i]
                      << std::endl;
          continue;
        }
        if (ind_set[i].c_str() == nullptr) {
          Rcpp::Rcout << "Null string found at idx: " << i << std::endl;
          continue;
        }
        ind_set_filtered.push_back(ind_set[i]);
      }
    }

    // remove from covar, pheno, and interactions
    covar_matrix.shed_rows(nan_idx, ind_set_filtered);
    pheno_matrix.shed_rows(nan_idx, ind_set_filtered);
    covar_poi_interaction_matrix.shed_rows(nan_idx, ind_set_filtered);
    no_interaction_matrix.shed_rows(nan_idx, ind_set_filtered);

    std::vector<std::string> strat_individuals(ind_set_filtered.size());
    std::transform(
        ind_set_filtered.begin(), ind_set_filtered.end(),
        strat_individuals.begin(), [&intersected_ind](const std::string &elem) {
          return intersected_ind[std::distance(
              intersected_ind.begin(),
              std::find(intersected_ind.begin(), intersected_ind.end(), elem))];
        });
    double nonconvergence_status = 0.0;
    double filtered_pois = 0.0;

    int num_parallel_poi_blocks =
        (int)std::ceil((double)num_poi / (double)chunk_size);

    auto start_time = std::chrono::high_resolution_clock::now();
    // allocate memory space for H5 file to read into
    for (int block = 0; block < num_parallel_poi_blocks; block++) {
      int start_chunk = block * chunk_size;
      int end_chunk = start_chunk + chunk_size;
      if (end_chunk >= num_poi) {
        end_chunk = num_poi;
      }
      std::vector<std::string> poi_names_chunk(poi_names.begin() + start_chunk,
                                               poi_names.begin() + end_chunk);

      POIMatrix poi(&bed_reader, config.maf_threshold, config.hwe_threshold,
                    config.POI_effect_type, config.POI_type);
      FRMatrix poi_matrix = poi.get_chunk(poi_individuals, poi_names_chunk);
      poi.filter_rows(poi_matrix, strat_individuals);
      auto end_time = std::chrono::high_resolution_clock::now();
      poi_reading_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();

      FRMatrix filtered = poi.filter_genotype(poi_matrix);

      if (filtered.data.size() > 0) {
        start_time = std::chrono::high_resolution_clock::now();
        filtered.write_summary(config.output_dir, "POI_Summary", stratum,
                               process_id);
        end_time = std::chrono::high_resolution_clock::now();
        file_writing_time +=
            (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time)
                .count();
      }

      filtered_pois += poi_matrix.data.n_cols;

      start_time = std::chrono::high_resolution_clock::now();

      FRResult result(covar_matrix, poi_matrix, no_interaction_matrix,
                      covar_poi_interaction_matrix);
      end_time = std::chrono::high_resolution_clock::now();
      memory_allocation_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();

#if !defined(__APPLE__) && !defined(__MACH__)
      omp_set_num_threads(num_threads);
#endif
      start_time = std::chrono::high_resolution_clock::now();
      std::unique_ptr<RegressionBase> regression;
      regression.reset(new LogisticRegression());
      regression->run_vla(covar_matrix, pheno_matrix, poi_matrix, result,
                          config.max_iter, config.p_value_type == "t.dist");
      end_time = std::chrono::high_resolution_clock::now();
      regression_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();
      start_time = std::chrono::high_resolution_clock::now();
      FRMatrix::write_vla_results(result, config.output_dir, "Results", stratum,
                                  config.output_exclude_covar, process_id + 1);

      end_time = std::chrono::high_resolution_clock::now();
      file_writing_time +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                start_time)
              .count();
      poi_matrix.col_names.clear();

      arma::fcolvec convergence = arma::conv_to<fcolvec>::from(
          (result.beta_rel_errs > config.rel_conv_tolerance) &&
          (result.beta_abs_errs > config.abs_conv_tolerance));
      nonconvergence_status = arma::sum(convergence);
      // proc_res.print_convergence_percentage(nonconvergence_status,
      // filtered_pois);
      proc_res.process_nonconvergence_status += nonconvergence_status;
      proc_res.process_total_filtered_pois += filtered_pois;
    }
  }

  proc_res.timing_results[0] = poi_reading_time;
  proc_res.timing_results[1] = file_writing_time;
  proc_res.timing_results[2] = memory_allocation_time;
  proc_res.timing_results[3] = regression_time;
  proc_res.print_timing_summary(process_id);
}

// [[Rcpp::export]]
Rcpp::List
FastRegVLA(const std::string phenotype, const std::string regression_type,
           const std::string pvalue_dist, bool output_exclude_covar,
           double maf_threshold, double hwe_threshold, bool no_intercept,
           double colinearity_rsq, int poi_block_size, int max_iter,
           double rel_conv_tolerance, double abs_conv_tolderance,
           int max_openmp_threads, const std::string pheno_file,
           const std::string pheno_rowname_cols,
           const std::string pheno_file_delim, const std::string covar_file,
           const std::string covar_rowname_cols,
           const std::string covar_file_delim, const std::string poi_file_dir,
           const std::string poi_file_delim, const std::string poi_file_format,
           const std::string poi_type, const std::string poi_effect_type,
           const Rcpp::StringVector covariates,
           const Rcpp::StringVector covariate_type,
           const Rcpp::LogicalVector covariate_standardize,
           const Rcpp::StringVector covariate_levels,
           const Rcpp::StringVector covariate_ref_level,
           const Rcpp::StringVector POI_covar_interactions_str,
           const Rcpp::StringVector split_by_str, const std::string output_dir,
           bool compress_results, int max_workers) {
  Config config(
      phenotype, regression_type, pvalue_dist, output_exclude_covar,
      maf_threshold, hwe_threshold, no_intercept, colinearity_rsq,
      poi_block_size, max_iter, rel_conv_tolerance, abs_conv_tolderance,
      max_openmp_threads, pheno_file, pheno_rowname_cols, pheno_file_delim,
      covar_file, covar_rowname_cols, covar_file_delim, poi_file_dir,
      poi_file_delim, poi_file_format, poi_type, poi_effect_type, covariates,
      covariate_type, covariate_standardize, covariate_levels,
      covariate_ref_level, POI_covar_interactions_str, split_by_str, output_dir,
      compress_results, max_workers);

  config.print();
  if (dir_exists(config.output_dir)) {
    delete_dir(config.output_dir);
  }
  CovariateMatrix cov_mat = CovariateMatrix(
      config.covar_file, config.covar_file_delim, config.covar_rowname_cols,
      config.covs, config.colinearity_rsq, config.no_intercept);

  FRMatrix covar_df = cov_mat.create_design_matrix();

  PhenoMatrix pheno_matrix(config.pheno_file, config.pheno_file_delim,
                           config.pheno_rowname_cols, config.phenotype);

  FRMatrix pheno_df = pheno_matrix.create_matrix();
  BEDReader bed_reader(config.poi_files[0]);

  std::vector<std::string> poi_names = bed_reader.get_names();
  std::vector<std::string> poi_individuals = bed_reader.get_individuals();
  std::vector<std::string> common_ind =
      intersect_row_names(pheno_df.sort_map(true), covar_df.sort_map(true));
  std::vector<std::string> intersected_ind =
      intersect_row_names(common_ind, poi_individuals);
  // FRMatrix Z = bed_reader.read_chunk(intersected_ind, poi_names);

  Rcout << intersected_ind.size()
        << " common unique subjects in pheno.file, "
           "covar.file, and POI.file"
        << std::endl;
  if (intersected_ind.empty()) {
    stop("No overlapping individuals found in POI, pheno, and covar files");
  }

  int num_poi = poi_names.size();
  if (num_poi == 0) {
    stop("No overlapping individuals found in POI, pheno, covar files");
  }

  // Stratify data
  int num_poi_files = config.poi_files.size();
  Strata stratums;
  stratums.stratify(config.split_by, covar_df, intersected_ind);
  Chunker chunker =
      Chunker(num_poi, intersected_ind.size(), config.max_openmp_threads,
              config.poi_block_size, num_poi_files, config.max_workers);

  // setup parallel processing
  // total_num_chunks
  int parallel_chunk_size =
      chunker.get_chunk_size() / 2; // twice the memory usage due to 2 fits
  int num_threads = chunker.get_openmp_threads();
  const int timing_results_size = 4;
  double concatenation_time, compression_time = 0.0;
  ProcResult total_proc_res;
  int total_timing_results[timing_results_size] = {0, 0, 0, 0};
#if !defined(__APPLE__) && !defined(__MACH__)
  omp_set_num_threads(num_threads);
#endif
#ifdef _WIN32
  for (int i = 0; i < num_poi_files; i++) {
    ProcResult proc_res;
    int timing_results[] = {0, 0, 0, 0};
    process_chunk_vla(i, config, pheno_df, covar_df, config.poi_files[i],
                      parallel_chunk_size, num_threads, false, proc_res);

    total_proc_res.accumulate(proc_res);
  }
#else

  int num_processes_total = chunker.get_total_workers();
  int max_processes = chunker.get_num_workers();
  std::vector<int> pipe_file_descriptors(num_processes_total * 2);
  std::vector<pid_t> process_ids(num_processes_total);
  std::vector<bool> has_completed(num_processes_total, false);

  int num_processes_started = 0;
  int num_processes_completed = 0;

  while (num_processes_completed < num_processes_total) {
    while ((num_processes_started - num_processes_completed) < max_processes) {
      checkInterrupt();
      if (num_processes_started == num_processes_total) {
        break;
      }
      int i = num_processes_started;
      if (pipe(&pipe_file_descriptors[i * 2]) == -1) {
        perror("pipe");
      }
      process_ids[i] = fork();
      if (process_ids[i] == -1) {
        perror("fork");
      }
      std::string poi_file_path = config.poi_files[i];
      if (process_ids[i] == 0) {             // child process
        close(pipe_file_descriptors[i * 2]); // close read pipe

        ProcResult proc_res;
        process_chunk_vla(i, config, pheno_df, covar_df, config.poi_files[i],
                          parallel_chunk_size, num_threads, false, proc_res);
        ssize_t res = write(pipe_file_descriptors[i * 2 + 1], &proc_res,
                            sizeof(proc_res));
        close(pipe_file_descriptors[i * 2 + 1]);
        _exit(EXIT_SUCCESS);
      } else {
        close(pipe_file_descriptors[i * 2 + 1]);
        num_processes_started++;
      }
    }
    // Check for finished processes
    for (int i = 0; i < num_processes_started; i++) {
      checkInterrupt();
      if (process_ids[i] != 0) { // parent process
        if (!has_completed[i] && waitpid(process_ids[i], NULL, WNOHANG) > 0) {
          has_completed[i] = true;

          ProcResult proc_res;
          ssize_t res =
              read(pipe_file_descriptors[i * 2], &proc_res, sizeof(proc_res));
          close(pipe_file_descriptors[i * 2]);
          total_proc_res.accumulate(proc_res);
          num_processes_completed++;
        }
      }
    }
  }
#endif
  auto start_time = std::chrono::high_resolution_clock::now();
  FRMatrix::concatenate_results(config.output_dir, "Results", "Full");
  FRMatrix::concatenate_results(config.output_dir, "Convergence", "Full");
  if (config.POI_type == "genotype") {
    FRMatrix::concatenate_results(config.output_dir, "POI_Summary", "Full");
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  concatenation_time =
      (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                    start_time)
          .count();

  if (config.compress_results) {
    start_time = std::chrono::high_resolution_clock::now();
    FRMatrix::zip_results(config.output_dir);
    end_time = std::chrono::high_resolution_clock::now();
    compression_time +=
        (double)std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time)
            .count();
    Rcpp::Rcout << "Results compression: " << compression_time / 1000.0 << "s"
                << std::endl;
  }
  total_proc_res.print_nonconvergence_summary();
  total_proc_res.print_totals_summary(
      concatenation_time, compression_time, config.regression_type,
      intersected_ind.size(), num_poi * num_poi_files, num_threads);
  // Rcpp::Rcout << "-----------------------------------------" << std::endl;
  // Rcpp::Rcout << "Convergence Summary: " << std::endl;
  // Rcpp::Rcout << "Reading HDF5: " << total_timing_results[0] / 1000.0 << "s";
  // Rcpp::Rcout << "-----------------------------------------" << std::endl;

  // Rcpp::List covar =
  //     arma_2_df(covar_df.data, covar_df.row_names_arr,
  //     covar_df.col_names_arr);
  // Rcpp::List phen =
  //     arma_2_df(pheno_df.data, pheno_df.row_names_arr,
  //     pheno_df.col_names_arr);
  // Rcpp::List poi = arma_2_df(Z.data, Z.row_names_arr, Z.col_names_arr);
}
