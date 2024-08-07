#define ARMA_WARN_LEVEL 0
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <chunker.h>
#include <config.h>
#include <covariate.h>
#include <ctime>
#include <fr_matrix.h>
#include <h5file.h>
#include <iostream>
#include <iterator>
#include <regression.h>
#include <strata.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

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

struct ProcResult
{
  int timing_results[4] = {0, 0, 0, 0};
  double process_nonconvergence_status = 0.0;
  double process_total_filtered_pois = 0.0;
};

void process_chunk(int process_id, Config &config, FRMatrix &pheno_df,
                   FRMatrix &covar_df, std::string poi_file_path,
                   int chunk_size, int num_threads, int timing_results[4])
{
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
  if (intersected_ind.empty())
  {
    stop("No overlapping individuals found in POI, pheno, and covar files");
  }

  // if (!config.POI_subset_file.empty()) {

  //     FRMatrix poi_subset(config.POI_subset_file,
  //     config.POI_subset_file_delim, config.POI_subset_rowname_col);
  //     std::vector<std::string> poi_names =
  //     intersect_row_names(poi_subset.str_data[0], poi.names);
  // }

  int num_poi = poi_names.size();
  if (num_poi == 0)
  {
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
  for (int stratum = 0; stratum < stratums.nstrata; ++stratum)
  {
    std::string outfile_suffix = stratums.ids[stratum];
    if (!config.split_by[0].empty())
    {
      Rcpp::Rcout << "Processing stratum: " << outfile_suffix.substr(1)
                  << std::endl;
    }
    std::vector<std::string> ind_set = stratums.index_list[outfile_suffix];
    int ct = 0;
    std::vector<int> ind_set_idx(ind_set.size());
    for (std::string ind : ind_set)
    {
      ind_set_idx[ct] = pheno_df.get_row_idx(ind);
      ct++;
    }
    
    // Rcpp::Rcout << "Init matrices" << std::endl;
    // initialize matrices
    FRMatrix pheno_matrix = pheno_df; // check col names

    // n individuals x # covariates
    FRMatrix covar_matrix = create_design_matrix(
        covar_df, config.covs, config.no_intercept, config.colinearity_rsq);

    FRMatrix covar_poi_interaction_matrix;
    covar_poi_interaction_matrix.data =
        arma::fmat(covar_matrix.data.n_rows, 1, arma::fill::ones);
    covar_poi_interaction_matrix.row_names = covar_matrix.row_names;
    covar_poi_interaction_matrix.col_names = {{"poi", 0}};
    // n individuals x 1 or 1 + num interacting poi covars
    create_Z_matrix(covar_matrix, config.POI_covar_interactions,
                    covar_poi_interaction_matrix);
    std::vector<int> nan_idx;
    std::vector<std::string> ind_set_filtered;
    
    // Rcpp::Rcout << "Identify missing" << std::endl;
    // identify missing values for covar, pheno matrix
    for (size_t i = 0; i < covar_matrix.data.n_rows; i++)
    {
      arma::uvec covar_nan_idx = arma::find_nonfinite(covar_matrix.data.row(i));
      arma::uvec pheno_nan_idx = arma::find_nonfinite(pheno_matrix.data.row(i));
      if (covar_nan_idx.size() > 0 || pheno_nan_idx.size() > 0)
      {
        nan_idx.push_back(i);
      }
      else
      {
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
    for (size_t j = 0; j < covar_matrix.data.n_rows; j++)
    {
      covar_matrix.row_names[ind_set_filtered[j]] = j;
      pheno_matrix.row_names[ind_set_filtered[j]] = j;
      covar_poi_interaction_matrix.row_names[ind_set_filtered[j]] = j;
    }
    std::vector<std::string> strat_individuals(ind_set_filtered.size());
    std::transform(
        ind_set_filtered.begin(), ind_set_filtered.end(),
        strat_individuals.begin(), [&intersected_ind](const std::string &elem)
        { return intersected_ind[std::distance(
              intersected_ind.begin(),
              std::find(intersected_ind.begin(), intersected_ind.end(), elem))]; });
    double nonconvergence_status = 0.0;
    double total_filtered_pois = 0.0;

    int num_parallel_poi_blocks =
        (int)std::ceil((double)num_poi / (double)chunk_size);
    // int total_nonconvergence_status = 0;
    // double sum_total_filtered_pois = 0.0;
    
    FRMatrix poi_matrix;
    auto start_time = std::chrono::high_resolution_clock::now();
    // allocate memory space for H5 file to read into
    for (int block = 0; block < num_parallel_poi_blocks; block++)
    {
      int start_chunk = block * chunk_size;
      int end_chunk = start_chunk + chunk_size;
      if (end_chunk >= num_poi)
      {
        end_chunk = num_poi;
      }
      
      // Rcpp::Rcout << "Reading poi chunk" << std::endl;
      std::vector<std::string> poi_names_chunk(poi_names.begin() + start_chunk,
                                               poi_names.begin() + end_chunk);
      // Rcpp::Rcout << "start_chunk: " << start_chunk << "\nend_chunk: " << end_chunk << "\nblock: " << block << "/" << num_parallel_poi_blocks << std::endl;
      poi.load_data_chunk(poi_matrix, poi.individuals, poi_names_chunk);
      // Rcpp::Rcout << "loaded chunk" << std::endl;
      std::vector<std::string> srt_cols_2 = poi_matrix.sort_map(false);
      std::vector<std::string> drop_rows =
          set_diff(poi.individuals, strat_individuals);

      int num_dropped = poi.individuals.size() - strat_individuals.size();
      arma::uvec drop_row_idx(drop_rows.size());
      for (size_t i = 0; i < drop_rows.size(); i++)
      {
        drop_row_idx[i] = poi_matrix.row_names[drop_rows[i]];
      }

      std::unordered_map<std::string, int> new_row_names(
          strat_individuals.size());
      for (auto &ind : strat_individuals)
      {
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

      if (config.POI_type == "genotype")
      {
        FRMatrix filtered =
            filter_poi(poi_matrix, config.maf_threshold, config.hwe_threshold);
        arma::uvec filtered_col = arma::find(filtered.data.row(5) == 0);

        if (filtered.data.n_cols == 0 ||
            filtered_col.n_elem == poi_matrix.data.n_cols)
        {
          Rcpp::Rcout << "no POI passed filtering" << std::endl;
          return;
        }

        std::vector<std::string> poi_col_names = filtered.sort_map(false);
        int cols_erased = 0;

        for (unsigned int i = 0; i < poi_col_names.size(); i++)
        {
          if ((unsigned)cols_erased < filtered_col.n_elem &&
              filtered_col[cols_erased] == i)
          {
            poi_matrix.col_names.erase(poi_col_names[i]);
            cols_erased++;
          }
          else
          {
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

      total_filtered_pois += poi_matrix.data.n_cols;

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
      se_beta.data =
          arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);

      FRMatrix neglog10_pvl;
      neglog10_pvl.data =
          arma::fmat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);

      for (auto &col_name : covar_matrix.col_names)
      {
        beta_est.row_names[col_name.first] = col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] =
            beta_est.row_names[col_name.first];
      }
      for (auto &col_name : covar_poi_interaction_matrix.col_names)
      {
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
      for (arma::uword v = 0; v < poi_matrix.data.n_cols; v++)
      {
        arma::uvec G_na = arma::find_nonfinite(poi_matrix.data.col(v));
        for (arma::uword i = 0; i < G_na.n_elem; i++)
        {
          W2(G_na(i), v) = 0;
          poi_matrix.data(G_na(i), v) = 0;
        }
      }

      end_time = std::chrono::high_resolution_clock::now();
      memory_allocation_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();

      // Rcpp::Rcout << "init matrices" << std::endl;
      start_time = std::chrono::high_resolution_clock::now();
      std::unique_ptr<RegressionBase> regression;
      // pheno_matrix.data.print();
      // int covar_nans = arma::sum(arma::find_nonfinite(covar_matrix.data));
      
      // int pheno_nans = arma::sum(arma::find_nonfinite(pheno_matrix.data));
      
      // int poi_nans =arma::sum( arma::find_nonfinite(poi_matrix.data));
      // int w2_nans = arma::sum(arma::find_nonfinite(W2));
      // int inter_nans = arma::sum(arma::find_nonfinite(covar_poi_interaction_matrix.data));

      // Rcpp::Rcout << "covar_nans: " << covar_nans << std::endl;
      // Rcpp::Rcout << "pheno_nans: " << pheno_nans << std::endl;
      // Rcpp::Rcout << "poi_nans: " << poi_nans << std::endl;
      // Rcpp::Rcout << "w2_nans: " << w2_nans << std::endl;
      // Rcpp::Rcout << "inter_nans: " << inter_nans << std::endl;


#if !defined(__APPLE__) && !defined(__MACH__)
      omp_set_num_threads(num_threads);
#endif
      if (config.regression_type == "logistic")
      {
        regression.reset(new LogisticRegression());
      }
      else
      {
        regression.reset(new LinearRegression());
      }

      // Rcpp::Rcout << "start regression" << std::endl;
      regression->run(covar_matrix, pheno_matrix, poi_matrix,
                      covar_poi_interaction_matrix, W2, beta_est, se_beta,
                      neglog10_pvl, beta_rel_errs, beta_abs_errs,
                      config.max_iter, config.p_value_type == "t.dist");
      end_time = std::chrono::high_resolution_clock::now();
      regression_time +=
          (double)std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time - start_time)
              .count();
      start_time = std::chrono::high_resolution_clock::now();
      FRMatrix::write_results(beta_est, se_beta, neglog10_pvl, W2, srt_cols,
                              config.output_dir, "Results", stratum,
                              config.output_exclude_covar, process_id + 1);
      end_time = std::chrono::high_resolution_clock::now();
      file_writing_time +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                start_time)
              .count();
      poi_matrix.col_names.clear();

      if (config.regression_type == "logistic")
      {
        start_time = std::chrono::high_resolution_clock::now();
        std::string convergence_file_prefix = "Convergence";
        FRMatrix::write_convergence_results(
            beta_est, srt_cols, config.output_dir, convergence_file_prefix,
            beta_rel_errs, beta_abs_errs, stratum, process_id + 1);
        end_time = std::chrono::high_resolution_clock::now();
      }

      arma::fcolvec convergence = arma::conv_to<fcolvec>::from(
          (beta_rel_errs > config.rel_conv_tolerance) &&
          (beta_abs_errs > config.abs_conv_tolerance));
      nonconvergence_status = arma::sum(convergence);
      double noncovergence_percent =
          (nonconvergence_status / total_filtered_pois) * 100;
      if (noncovergence_percent > 0.0)
      {
        Rcpp::Rcout << nonconvergence_status << " out of "
                    << total_filtered_pois << " (" << std::setprecision(2)
                    << noncovergence_percent
                    << "%) POIs did not meet relative and absolute convergence "
                       "threshold."
                    << std::endl;
        Rcpp::Rcout << "See convergence_" << stratum
                    << ".tsv for additional details." << std::endl;
      }
    }
  }

  poi.close_all();

  auto end =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  timing_results[0] = poi_reading_time;
  timing_results[1] = file_writing_time;
  timing_results[2] = memory_allocation_time;
  timing_results[3] = regression_time;
  Rcpp::Rcout << "Timing Summary for process: " << process_id + 1 << std::endl;
  Rcpp::Rcout << "Reading HDF5: " << poi_reading_time / 1000.0 << "s"
              << std::endl;
  Rcpp::Rcout << "Writing results: " << file_writing_time / 1000.0 << "s"
              << std::endl;
  Rcpp::Rcout << "Memory allocation: " << memory_allocation_time / 1000.0 << "s"
              << std::endl;
  Rcpp::Rcout << "Regression: " << regression_time / 1000.0 << "s" << std::endl;
  Rcpp::Rcout << "Completed process " << process_id + 1
              << " at: " << std::ctime(&end) << std::endl;
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
    bool compress_results, int max_workers)
{
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
  Rcpp::Rcout << "-----------------------------------------" << std::endl;
  Rcpp::Rcout << "Running FastReg with configuration: " << std::endl;
  Rcpp::Rcout << "phenotype: " << phenotype << std::endl;
  Rcpp::Rcout << "regression_type: " << regression_type << std::endl;
  Rcpp::Rcout << "pvalue_dist: " << pvalue_dist << std::endl;
  Rcpp::Rcout << "output_exclude_covar: " << output_exclude_covar << std::endl;
  Rcpp::Rcout << "maf_threshold: " << maf_threshold << std::endl;
  Rcpp::Rcout << "hwe_threshold: " << hwe_threshold << std::endl;
  Rcpp::Rcout << "no_intercept: " << no_intercept << std::endl;
  Rcpp::Rcout << "colinearity_rsq: " << colinearity_rsq << std::endl;
  Rcpp::Rcout << "poi_block_size: " << poi_block_size << std::endl;
  Rcpp::Rcout << "max_iter: " << max_iter << std::endl;
  Rcpp::Rcout << "rel_conv_tolerance: " << rel_conv_tolerance << std::endl;
  Rcpp::Rcout << "abs_conv_tolderance: " << abs_conv_tolderance << std::endl;
  Rcpp::Rcout << "max_openmp_threads: " << max_openmp_threads << std::endl;
  Rcpp::Rcout << "pheno_file: " << pheno_file << std::endl;
  Rcpp::Rcout << "pheno_rowname_cols: " << pheno_rowname_cols << std::endl;
  Rcpp::Rcout << "pheno_file_delim: " << pheno_file_delim << std::endl;
  Rcpp::Rcout << "covar_file: " << covar_file << std::endl;
  Rcpp::Rcout << "covar_rowname_cols: " << covar_rowname_cols << std::endl;
  Rcpp::Rcout << "covar_file_delim: " << covar_file_delim << std::endl;
  Rcpp::Rcout << "poi_file_dir: " << poi_file_dir << std::endl;
  Rcpp::Rcout << "poi_file_delim: " << poi_file_delim << std::endl;
  Rcpp::Rcout << "poi_file_format: " << poi_file_format << std::endl;
  Rcpp::Rcout << "poi_type: " << poi_type << std::endl;
  Rcpp::Rcout << "poi_effect_type: " << poi_effect_type << std::endl;
  Rcpp::Rcout << "covariates: " << covariates << std::endl;
  Rcpp::Rcout << "covariate_type: " << covariate_type << std::endl;
  Rcpp::Rcout << "covariate_standardize: " << covariate_standardize
              << std::endl;
  Rcpp::Rcout << "covariate_levels: " << covariate_levels << std::endl;
  Rcpp::Rcout << "covariate_ref_level: " << covariate_ref_level << std::endl;
  Rcpp::Rcout << "POI_covar_interactions_str: " << POI_covar_interactions_str
              << std::endl;
  Rcpp::Rcout << "split_by_str: " << split_by_str << std::endl;
  Rcpp::Rcout << "output_dir: " << output_dir << std::endl;
  Rcpp::Rcout << "compress_results: " << compress_results << std::endl;
  Rcpp::Rcout << "max_workers: " << max_workers << std::endl;
  Rcpp::Rcout << "-----------------------------------------" << std::endl;
  // Clean up previous run
  if (dir_exists(config.output_dir))
  {
    delete_dir(config.output_dir);
  }

  FRMatrix pheno_df(config.pheno_file, config.pheno_file_delim,
                    config.pheno_rowname_cols, config.phenotype);
  FRMatrix covar_df(config.covar_file, config.covar_file_delim,
                    config.covar_rowname_cols, config.covariates, config.covariate_type);

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
  if (intersected_ind.empty())
  {
    stop("No overlapping individuals found in POI, pheno, and covar files");
  }

  int num_poi = poi_names.size();
  if (num_poi == 0)
  {
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
  int total_timing_results[timing_results_size] = {0, 0, 0, 0};
#ifdef _WIN32
  for (int i = 0; i < num_poi_files; i++)
  {
    int timing_results[] = {0, 0, 0, 0};
    process_chunk(i, config, pheno_df, covar_df, poi_file_path,
                  parallel_chunk_size, num_threads, timing_results);
    for (int j = 0; j < timing_results_size; j++)
    {
      total_timing_results[j] += timing_results[j];
    }
  }
#else

  int num_processes_total = chunker.get_total_workers();
  int max_processes = chunker.get_num_workers();
  std::vector<int> pipe_file_descriptors(num_processes_total * 2);
  std::vector<pid_t> process_ids(num_processes_total);
  std::vector<bool> has_completed(num_processes_total, false);

  int num_processes_started = 0;
  int num_processes_completed = 0;

  while (num_processes_completed < num_processes_total)
  {
    while ((num_processes_started - num_processes_completed) < max_processes)
    {
      checkInterrupt();
      if (num_processes_started == num_processes_total)
      {
        break;
      }
      int i = num_processes_started;
      if (pipe(&pipe_file_descriptors[i * 2]) == -1)
      {
        perror("pipe");
        return;
      }
      process_ids[i] = fork();
      if (process_ids[i] == -1)
      {
        perror("fork");
        return;
      }
      std::string poi_file_path = config.poi_files[i];
      if (process_ids[i] == 0)
      {                                      // child process
        close(pipe_file_descriptors[i * 2]); // close read pipe
        int timing_results[timing_results_size] = {0, 0, 0, 0};
        // Rcpp::Rcout << "Started processing for " << i + 1 << std::endl;
        process_chunk(i, config, pheno_df, covar_df, poi_file_path,
                      parallel_chunk_size, num_threads, timing_results);
        ssize_t res = write(pipe_file_descriptors[i * 2 + 1], timing_results, sizeof(timing_results));
        close(pipe_file_descriptors[i * 2 + 1]);
        _exit(EXIT_SUCCESS);
        return;
      }
      else
      {
        close(pipe_file_descriptors[i * 2 + 1]);
        num_processes_started++;
      }
    }
    // Check for finished processes
    for (int i = 0; i < num_processes_started; i++)
    {
      checkInterrupt();
      if (process_ids[i] != 0)
      { // parent process
        if (!has_completed[i] && waitpid(process_ids[i], NULL, WNOHANG) > 0)
        {
          has_completed[i] = true;
          int timing_results[timing_results_size] = {0, 0, 0, 0};
          ssize_t res = read(pipe_file_descriptors[i * 2], timing_results, sizeof(timing_results));
          close(pipe_file_descriptors[i * 2]);

          for (int j = 0; j < timing_results_size; j++)
          {
            total_timing_results[j] += timing_results[j];
          }
          num_processes_completed++;
        }
      }
    }
  }
#endif
  auto start_time = std::chrono::high_resolution_clock::now();
  FRMatrix::concatenate_results(config.output_dir, "Results", "Full");
  FRMatrix::concatenate_results(config.output_dir, "Convergence", "Full");
  if (config.POI_type == "genotype")
  {
    FRMatrix::concatenate_results(config.output_dir, "POI_Summary", "Full");
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  concatenation_time = (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_time - start_time)
                           .count();

  // timing_results[0] = poi_reading_time;
  // timing_results[1] = file_writing_time;
  // timing_results[2] = memory_allocation_time;
  // timing_results[3] = regression_time;

  Rcpp::Rcout << "-----------------------------------------" << std::endl;
  Rcpp::Rcout << "Timing Summary: " << std::endl;
  Rcpp::Rcout << "Reading HDF5: " << total_timing_results[0] / 1000.0 << "s"
              << std::endl;
  Rcpp::Rcout << "Writing results: " << total_timing_results[1] / 1000.0 << "s"
              << std::endl;
  Rcpp::Rcout << "Memory allocation: " << total_timing_results[2] / 1000.0 << "s"
              << std::endl;
  Rcpp::Rcout << "Regression: " << total_timing_results[3] / 1000.0 << "s" << std::endl;
  Rcpp::Rcout << "Results concatenation: " << concatenation_time / 1000.0 << "s" << std::endl;

  if (config.compress_results)
  {
    start_time = std::chrono::high_resolution_clock::now();
    FRMatrix::zip_results(config.output_dir);
    end_time = std::chrono::high_resolution_clock::now();
    compression_time += (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                            end_time - start_time)
                            .count();
    Rcpp::Rcout << "Results compression: " << concatenation_time / 1000.0 << "s" << std::endl;
  }
  auto end =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  Rcpp::Rcout << "Completed " << config.regression_type << " regression -" << std::endl;
  Rcpp::Rcout << "\t\tnum individuals: " << intersected_ind.size() << std::endl;
  Rcpp::Rcout << "\t\t~num pois: " << num_poi * num_poi_files << std::endl;
  Rcpp::Rcout << "\t\twith openmp thread(s): " << num_threads << std::endl;
  Rcpp::Rcout << "at: " << std::ctime(&end) << std::endl;
  Rcpp::Rcout << "-----------------------------------------" << std::endl;
}
