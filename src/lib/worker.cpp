#include <worker.h>

using namespace Rcpp;
using namespace arma;

void Worker::operator()() {
    while(true) {
        Job job;
        {
            std::unique_lock<std::mutex> lock(_mtx);
            _cv.wait(lock, [&] { return !_job_queue.empty();});
            if(_job_queue.empty()) {
                return;
            }
            job = _job_queue.front();
            _job_queue.pop();
        }
        process_job(job);
        _file_writing_time += job.file_writing_time;
        _memory_allocation_time += job.memory_allocation_time;
        _poi_reading_time += job.poi_reading_time;
        _regression_time += job.regression_time;
        // _nonconvergence_status_accumulator += job.nonconvergence_status;
        // _total_filtered_pois_accumulator += job.total_filtered_pois;
    }
}

void Worker::process_job(Job& job) {
    FRMatrix poi_matrix;
    H5File poi(job.config.POI_file);
    std::cout << "Processing POIs block: " << job.id + 1 << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    //{
        // reading HDF5 file in critical section
    // std::unique_lock<std::mutex> lock(_mtx);
    poi.open_file(true);
    poi.get_values_dataset_id();
    poi.get_POI_individuals();
    poi.get_POI_names();
    poi.set_memspace(poi.individuals.size(), job.chunk_size);
    int start_chunk = job.id*job.chunk_size;
    int end_chunk = start_chunk + job.chunk_size;
    if (end_chunk > job.num_poi) {
        end_chunk = job.num_poi;
        poi.set_memspace(poi.individuals.size(), end_chunk - start_chunk);
    }
    std::vector<std::string> poi_names_chunk(job.poi_names.begin() + start_chunk, job.poi_names.begin() + end_chunk);
    poi.get_POI_matrix(poi_matrix, poi.individuals, poi_names_chunk, job.chunk_size);
    poi.close_all();
    //}
    
    #if !defined(__APPLE__) && !defined(__MACH__)
    omp_set_num_threads(job.num_threads);
    #endif
    std::vector<std::string> srt_cols_2 = poi_matrix.sort_map(false);
    std::vector<std::string> drop_rows = set_diff(poi.individuals, job.strat_individuals);

    int num_dropped = poi.individuals.size() - job.strat_individuals.size();
    arma::uvec drop_row_idx(drop_rows.size());
    for (size_t i = 0; i < drop_rows.size(); i++) {
        drop_row_idx[i] = poi_matrix.row_names[drop_rows[i]];
    }
    
    std::unordered_map<std::string, int> new_row_names(job.strat_individuals.size());
    for(auto& ind : job.strat_individuals) {
        new_row_names[ind] = poi_matrix.row_names[ind] - num_dropped;
    }
    poi_matrix.data.shed_rows(drop_row_idx);
    poi_matrix.row_names = new_row_names;
    
    srt_cols_2 = poi_matrix.sort_map(false);
    auto end_time = std::chrono::high_resolution_clock::now();
    job.poi_reading_time = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    std::cout << "Reading POI timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";

    if (job.config.POI_type == "genotype") {
        std::cout << "Filtering MAF and HWE" << std::endl;

        FRMatrix filtered = filter_poi(poi_matrix, job.config.maf_threshold, job.config.hwe_threshold);
        arma::uvec filtered_col = arma::find(filtered.data.row(5) == 0);

        if (filtered.data.n_cols == 0 || filtered_col.n_elem == poi_matrix.data.n_cols) {
            std::cout << "no POI passed filtering" << std::endl;
            return;
        }
        
        std::vector<std::string> poi_col_names = filtered.sort_map(false);
        int cols_erased = 0;

        for(unsigned int i = 0; i < poi_col_names.size(); i++) {
            if ((unsigned) cols_erased < filtered_col.n_elem && filtered_col[cols_erased] == i) {
                poi_matrix.col_names.erase(poi_col_names[i]);
                cols_erased++;
            }
            else {
                poi_matrix.col_names[poi_col_names[i]] = poi_matrix.col_names[poi_col_names[i]] - cols_erased;
            }
        }

        poi_matrix.data.shed_cols(filtered_col);

        std::cout << "transforming POI and writing statistics" << std::endl;
        srt_cols_2 = poi_matrix.sort_map(false);
        transform_poi(poi_matrix, job.config.POI_effect_type);
        start_time = std::chrono::high_resolution_clock::now();
        std::string summary_name = "POI_Summary_" + std::to_string(job.id);
        filtered.write_summary(job.config.output_dir, summary_name, job.stratum);
        end_time = std::chrono::high_resolution_clock::now();
        job.file_writing_time = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
        std::cout << "Wrting POI Summary timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
        std::cout << "Effective block size after filtering: " << poi_matrix.data.n_cols << std::endl;
    }

    job.total_filtered_pois = poi_matrix.data.n_cols;
    
    start_time = std::chrono::high_resolution_clock::now();
    // 1 phenotype 
    FRMatrix beta_est;
    FRMatrix se_beta;
    int num_parms = job.interactions_matrix.data.n_cols + job.covar_matrix.data.n_cols;
    beta_est.data = arma::mat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
    arma::colvec beta_rel_errs = arma::colvec(poi_matrix.data.n_cols, arma::fill::zeros);
    arma::colvec beta_abs_errs = arma::colvec(poi_matrix.data.n_cols, arma::fill::zeros);
    se_beta.data = arma::mat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);
    
    FRMatrix neglog10_pvl;
    neglog10_pvl.data = arma::mat(num_parms, poi_matrix.data.n_cols, arma::fill::zeros);


    for (auto& col_name : job.covar_matrix.col_names) {
        beta_est.row_names[col_name.first] = col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] = beta_est.row_names[col_name.first];
    }
    for (auto& col_name : job.interactions_matrix.col_names) {
        beta_est.row_names[col_name.first] = job.covar_matrix.col_names.size() + col_name.second;
        se_beta.row_names[col_name.first] = beta_est.row_names[col_name.first];
        neglog10_pvl.row_names[col_name.first] = beta_est.row_names[col_name.first];
    }

    beta_est.col_names = job.covar_matrix.row_names;
    se_beta.col_names = beta_est.col_names;
    neglog10_pvl.col_names = beta_est.col_names;
    std::vector<std::string> srt_cols = poi_matrix.sort_map(false);


    // create weight mask matrix    
    arma::umat W2 = arma::umat(poi_matrix.data.n_rows, poi_matrix.data.n_cols, arma::fill::ones);
    for (arma::uword v = 0; v < poi_matrix.data.n_cols; v++) {
        arma::uvec G_na = arma::find_nonfinite(poi_matrix.data.col(v));
        for (arma::uword i = 0; i < G_na.n_elem; i++) {
            W2(G_na(i), v) = 0;
            poi_matrix.data(G_na(i), v) = 0;
        }
    }
    
    end_time = std::chrono::high_resolution_clock::now();
    job.memory_allocation_time = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    std::cout << "Memory allocation timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";

    start_time = std::chrono::high_resolution_clock::now();
    std::unique_ptr<RegressionBase> regression;


    if (job.config.regression_type == "logistic") {
        std::cout << "Started logistic regression" << std::endl;
        regression.reset(new LogisticRegression());
    }
    else {
        std::cout << "Started linear regression" << std::endl;
        regression.reset(new LinearRegression ());
    }
    regression->run(
        job.covar_matrix, 
        job.pheno_matrix, 
        poi_matrix, 
        job.interactions_matrix, 
        W2, 
        beta_est, 
        se_beta, 
        neglog10_pvl, 
        beta_rel_errs,
        beta_abs_errs,
        job.config.max_iter, 
        job.config.p_value_type == "t.dist"
    );
    end_time = std::chrono::high_resolution_clock::now();
    Rcpp::Rcout << "regression timing for block " << job.id + 1 <<  "/" << job.num_poi_blocks << ": " <<std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
    job.regression_time = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    start_time = std::chrono::high_resolution_clock::now();
    std::string results_file_prefix = "Results_" + std::to_string(job.id);
    FRMatrix::write_results(beta_est, se_beta, neglog10_pvl, W2, srt_cols, job.config.output_dir, results_file_prefix, job.stratum, job.config.output_exclude_covar);
    end_time = std::chrono::high_resolution_clock::now();
    job.file_writing_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    std::cout << "Writing results for block " << job.id + 1 <<  "/" << job.num_poi_blocks << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " milliseconds\n";
    poi_matrix.col_names.clear();

    if (job.config.regression_type == "logistic") {
        start_time = std::chrono::high_resolution_clock::now();
        std::string convergence_file_prefix = "convergence_" + std::to_string(job.id);
        FRMatrix::write_convergence_results(beta_est, srt_cols, job.config.output_dir, convergence_file_prefix, beta_rel_errs, beta_abs_errs, job.stratum);
        end_time = std::chrono::high_resolution_clock::now();
        job.file_writing_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count();
    }

    arma::colvec convergence = conv_to<colvec>::from((beta_rel_errs > job.config.rel_conv_tolerance) && (beta_abs_errs > job.config.abs_conv_tolerance));
    job.nonconvergence_status = arma::sum(convergence);
    auto end = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Completed process " << job.id << " at: " << std::ctime(&end) << std::endl;
}