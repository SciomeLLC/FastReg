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
    FRMatrix interactions;
    interactions.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    interactions.row_names = covar_matrix.row_names;
    interactions.col_names = {{"poi", 0}};
    interactions.col_names_arr.push_back("poi");

    FRMatrix interactions_sqrd;
    interactions_sqrd.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    interactions_sqrd.row_names = covar_matrix.row_names;
    interactions_sqrd.col_names = {{"poi^2", 0}};
    interactions_sqrd.col_names_arr.push_back("poi^2");

    FRMatrix no_interaction_matrix;
    no_interaction_matrix.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    no_interaction_matrix.row_names = covar_matrix.row_names;
    no_interaction_matrix.col_names = {{"poi", 0}};
    no_interaction_matrix.col_names_arr.push_back("poi");

    // n individuals x 1 or 1 + num interacting poi covars
    create_interactions(covar_matrix, config.POI_covar_interactions,
                        interactions, !config.no_intercept, "poi");
    create_interactions(covar_matrix, config.POI_covar_interactions,
                        interactions_sqrd, !config.no_intercept, "poi^2");

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
    interactions.shed_rows(nan_idx, ind_set_filtered);
    interactions_sqrd.shed_rows(nan_idx, ind_set_filtered);
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

    // Rcpp::Rcout << "starting block loop" << std::endl;
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
                      interactions, interactions_sqrd);
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
      result.write_to_file(config.output_dir, "Results", stratum,
                           process_id + 1);
      // FRMatrix::write_vla_results(result, config.output_dir, "Results",
      // stratum,
      //                             config.output_exclude_covar, process_id +
      //                             1);

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
  FRResult::concatenate(config.output_dir, "Results", "Full");
  if (config.POI_type == "genotype") {
    FRResult::concatenate(config.output_dir, "POI_Summary", "Full");
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

  Rcpp::DataFrame covar =
      arma_2_df(covar_df.data, covar_df.row_names_arr, covar_df.col_names_arr);
  Rcpp::DataFrame phen =
      arma_2_df(pheno_df.data, pheno_df.row_names_arr, pheno_df.col_names_arr);
  // Rcpp::List poi = arma_2_df(Z.data, Z.row_names_arr, Z.col_names_arr);
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("covar") = covar,
                                      Rcpp::Named("pheno") = phen);
  return res;
}
