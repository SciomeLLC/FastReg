#include <fr_matrix.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <list>
#include <algorithm>
#include <regex>
#include <sys/stat.h>
#include <chrono>

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#  if __cplusplus >= 201703L && __has_include(<filesystem>)
#    include <filesystem>
namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#  endif
#endif

int FRMatrix::get_row_idx(const std::string& row_name) {
    bool exists = row_names.find(row_name) != row_names.end();
    if(!exists) {
      return -1;
    }
    return row_names[row_name];
}

int FRMatrix::get_col_idx(const std::string& col_name) {
    bool exists = col_names.find(col_name) != col_names.end();
    if(!exists) {
      return -1;
    }
    return col_names[col_name];
}

bool FRMatrix::validate_cols(std::string &names) {
    char delim = ';';

    std::vector<std::string> id_cols = split(names, delim);
    for(std::string &col : id_cols) {
        if(col_names.find(col) == col_names.end()) {
            return false;
        }
    }
    return true;
}

void FRMatrix::load_from_csv(std::string& filename, std::string& delim, std::string& id) {
    fs::path df_file_path(filename);
    if (!fs::exists(df_file_path)) {
        Rcpp::stop("%s does not exist.", filename);
    }
    std::ifstream file_stream(filename);
    std::string line;

    std::getline(file_stream, line);
    // remove carriage returns for windows
    if (line[line.size() - 1] == '\r') {
        line.erase(line.size() - 1);
    }
    std::stringstream ss(line);

    char delim_char;
    if (delim == "tab") {
        delim_char = '\t';
    } else if (delim == "comma") {
        delim_char = ',';
    } else if (delim == "semicolon") {
        delim_char = ';';
    } else {
        Rcpp::stop("Invalid delim! delim for %s must be 'tab', 'comma', or 'semicolon'", filename);
    }

    std::vector<std::string> col_headers = split(line, delim_char);

    int col_idx = 0;
    for (const std::string &col_name : col_headers) {
        if(col_name != id) {
            col_names[col_name] = col_idx;
            col_idx++;
        }
    }

    // Count lines
    size_t num_lines = 0;
    while (std::getline(file_stream, line)) {
        num_lines++;
    }

    // Reset file pointer to the beginning
    file_stream.clear();
    file_stream.seekg(0, std::ios::beg);

    // Skip header line
    std::getline(file_stream, line);

    // Initialize Armadillo matrix
    data.set_size(num_lines, col_idx + 1);
    str_data.resize(num_lines);
    std::vector<std::string> row_items;
    // Read file content
    size_t row_count = 0;
    while (std::getline(file_stream, line)) {
        row_items = split(line, delim_char);

        std::string row_name = row_items.front();
        if (row_names.find(row_name) != row_names.end()) {
            Rcpp::Rcout << "Found duplicate row name: " << row_name << std::endl;
            continue; // skip rows with duplicate names
        }
        row_names[row_name] = row_count;
        
        int str_col_count = 0;
        for (size_t i = 1; i < row_items.size(); ++i) {
            try {
                data(row_count, i - 1) = std::stod(row_items[i]);
            } catch (const std::invalid_argument &e) {
                str_data[row_count].push_back(row_items[i]);
                col_names_str[col_headers.at(i)] = str_col_count;
                str_col_count++; 
            }
        }

        row_count++;
    }
    file_stream.close();
}

FRMatrix FRMatrix::get_submat_by_cols(const std::vector<int>& row_idx, const std::vector<std::string>& names) {
    std::vector<int> col_indices;
    FRMatrix subm;

    // Iterate over the names vector and retrieve the corresponding column indices
    for (const std::string& item : names) {
        auto itr = col_names.find(item);
        if (itr != col_names.end()) {
            col_indices.push_back(itr->second);
            subm.col_names[item] = itr->second;
        }
    }

    // Extract the submatrix using the row and column indices
    subm.data = data.submat(arma::conv_to<arma::uvec>::from(row_idx), arma::conv_to<arma::uvec>::from(col_indices));

    // Assign the row names of the submatrix
    for (const int idx : row_idx) {
        auto itr = std::find_if(row_names.begin(), row_names.end(), [&](const std::pair<std::string, int>& p) {
            return p.second == idx;
        });

        if (itr != row_names.end()) {
            subm.row_names[itr->first] = idx;
        }
    }

    return subm;
}

FRMatrix FRMatrix::get_submat_by_cols(const std::vector<int>& row_idx, const std::unordered_map<std::string, int>& names) {
    std::vector<int> col_indices;
    FRMatrix subm;

    // Iterate over the names vector and retrieve the corresponding column indices
    for (auto& item : names) {
        auto itr = col_names.find(item.first);
        if (itr != col_names.end()) {
            col_indices.push_back(itr->second);
            subm.col_names[item.first] = itr->second;
        }
    }

    // Extract the submatrix using the row and column indices
    subm.data = data.submat(arma::conv_to<arma::uvec>::from(row_idx), arma::conv_to<arma::uvec>::from(col_indices));

    // Assign the row names of the submatrix
    for (const int idx : row_idx) {
        auto itr = std::find_if(row_names.begin(), row_names.end(), [&](const std::pair<std::string, int>& p) {
            return p.second == idx;
        });

        if (itr != row_names.end()) {
            subm.row_names[itr->first] = idx;
        }
    }

    return subm;
}

std::vector<std::string> FRMatrix::split(const std::string& str_tokens, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str_tokens);
    std::string item;
    while (std::getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

std::vector<std::string> FRMatrix::get_col_str(const std::string& col_name) {
    bool exists = col_names_str.find(col_name) != col_names_str.end();
    if(!exists) {
      Rcpp::stop("Column name " + col_name + " doesn't match any non-numeric columns.");
    }
    int col_idx = col_names_str[col_name];
    std::vector<std::string> col_vals;
    int count = 0;
    for (std::vector<std::string> row : str_data) {
        col_vals.push_back(row[col_idx]);
        count++;
    }

    return col_vals;
}


std::vector<std::string> FRMatrix::sort_map(bool rows) {
    std::unordered_map<std::string, int> temp;
    if (rows) {
        temp = row_names;
    }
    else {
        temp = col_names;
    }
    
    std::vector<std::string> sorted_arr(temp.size());
    int count = 0;
    for (auto item : temp) {
        sorted_arr[item.second] = item.first;
        count++;
    }
    return sorted_arr;
}

void FRMatrix::write_summary(std::string dir, std::string name, int stratum) {
    fs::create_directory(dir);
    std::stringstream ss;
    ss << dir << "/" << name << "_stratum_" << stratum + 1 << ".tsv";
    std::string file_name = ss.str();

    arma::mat temp = data.t();
    std::ofstream outfile;

    if(fs::exists(file_name)) {
        outfile = std::ofstream(file_name, std::ios::app);
    } else {
        outfile = std::ofstream(file_name);
        std::vector<std::string> ordered_row_names = sort_map(true);
        for (const auto& row_name : ordered_row_names) {
            outfile << '\t' << row_name;
        }
        outfile << std::endl;
    }
    
    outfile << std::fixed << std::setprecision(17);
    std::stringstream buffer;
    std::vector<std::string> ordered_col_names = sort_map(false);
    int row_idx = 0;
    for (const auto& col_name : ordered_col_names) {
        buffer << col_name;
        for (size_t col_idx = 0; col_idx < temp.n_cols; ++col_idx) {
            buffer << '\t' << temp(row_idx, col_idx);
        }
        buffer << std::endl;
        row_idx++;
    }
    outfile << buffer.str();
    outfile.close();
}

void FRMatrix::print() {
    std::vector<std::string> sort_cols = sort_map(false);
    std::vector<std::string> sort_rows = sort_map(true);
    for (auto& col : sort_cols) {
        Rcpp::Rcout << col << "\t";
    }
    Rcpp::Rcout << std::endl;

    for (size_t i = 0; i < data.n_rows; i++) {
        Rcpp::Rcout << sort_rows[i] << "\t";
        data.row(i).print();
    }
}

void FRMatrix::write_results(
    FRMatrix& beta, 
    FRMatrix& se_beta, 
    FRMatrix& neglog10, 
    arma::umat& W2, 
    std::vector<std::string> poi_names, 
    std::string dir, 
    std::string file_name, 
    int stratum, 
    bool exclude_covars,
    int process_id) 
{
    int n_parms = beta.data.n_rows;

    std::vector<std::string> sorted_row_names = beta.sort_map(true);
    // create dir if it doesn't exist
    fs::create_directory(dir);
    // Rcpp::Rcout << "Dir created or already exists" << std::endl;
    if (poi_names.size() != beta.data.n_cols) {
        Rcpp::Rcout << "Error: The size of poi_names does not match the number of columns in the beta matrix." << std::endl;
        // return;
    }

    std::stringstream ss;
    ss << dir << "/" << file_name << "_stratum_" << stratum + 1 << "_" << process_id << ".tsv";
    std::string result_file = ss.str();
    std::ofstream outfile;

    if(fs::exists(result_file)) {
        outfile.open(result_file, std::ios::app);
        if (!outfile.is_open()) {
            Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file << std::endl;
            return;
        }
        // Rcpp::Rcout << "File already exists and opened for writing/appending." << std::endl;
    } else {
        outfile.open(result_file);
        if (!outfile.is_open()) {
            Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file << std::endl;
            return;
        }
        outfile << "POI\tN\tDF\tEffect\tEstimate\tStd Error\tNegLog10 P-val" << std::endl;
        // Rcpp::Rcout << "File created for writing." << std::endl;
    }

    outfile << std::fixed << std::setprecision(17);
    
    std::stringstream buffer;
    for (int col = 0; col < (int)beta.data.n_cols; col++) {
        std::string poi_name = poi_names[col];
        arma::uvec w2_col = W2.col(col);
        int N = arma::as_scalar(arma::sum(w2_col, 0));
        int df = N - n_parms;
        int adder = 1;
        int row = 0;
        if (exclude_covars) {
            adder = n_parms;
            row = n_parms - 1;
        }

        for (; row < (int)beta.data.n_rows; row += adder) {
            std::string effect_name = sorted_row_names[row];
            double estimate = beta.data.at(row, col);
            double std_error = se_beta.data.at(row, col);
            double neglog10_pval = neglog10.data.at(row, col);
            
            buffer << poi_name << "\t" << N << "\t" << df << "\t" << effect_name << "\t" << estimate << "\t" << std_error << "\t" << neglog10_pval << "\t" << std::endl;
        }
    }
    outfile << buffer.str();
    outfile.close();
}

void FRMatrix::write_convergence_results(
    FRMatrix& beta, 
    std::vector<std::string> poi_names, 
    std::string dir, 
    std::string file_name, 
    arma::colvec& rel_err,
    arma::colvec& abs_err,
    int stratum,
    int process_id
    ) 
{
    std::vector<std::string> sorted_row_names = beta.sort_map(true);
    // create dir if it doesn't exist
    fs::create_directory(dir);
    // Rcpp::Rcout << "Dir created or already exists" << std::endl;
    if (poi_names.size() != beta.data.n_cols) {
        Rcpp::Rcout << "Error: The size of poi_names does not match the number of columns in the beta matrix." << std::endl;
        // return;
    }

    std::stringstream ss;
    ss << dir << "/" << file_name << "_stratum_" << stratum + 1 << "_" << process_id << ".tsv";
    std::string result_file = ss.str();
    std::ofstream outfile;

    if(fs::exists(result_file)) {
        outfile.open(result_file, std::ios::app);
        if (!outfile.is_open()) {
            Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file << std::endl;
            return;
        }
        // Rcpp::Rcout << "File already exists and opened for writing/appending." << std::endl;
    } else {
        outfile.open(result_file);
        if (!outfile.is_open()) {
            Rcpp::Rcout << "Error: Unable to open file for writing: " << result_file << std::endl;
            return;
        }
        outfile << "POI\tAbs Err\tRel Err" << std::endl;
        // Rcpp::Rcout << "File created for writing." << std::endl;
    }

    outfile << std::fixed << std::setprecision(17);
    
    std::stringstream buffer;
    for (int col = 0; col < (int)beta.data.n_cols; col++) {
        std::string poi_name = poi_names[col];
        
        double abs_err_val = abs_err.at(col);
        double rel_err_val = rel_err.at(col);
        
        buffer << poi_name << "\t" << abs_err_val << "\t" << rel_err_val << std::endl;
    }
    outfile << buffer.str();
    outfile.close();
}

void FRMatrix::zip_results(std::string output_dir) {
    Rcpp::Environment utils_env("package:utils");
    Rcpp::Function zip = utils_env["zip"];
    if(fs::exists(output_dir)) {
        std::string parent_path = fs::path(output_dir).parent_path().string();
        const auto time_now = std::chrono::system_clock::now();
        const auto time_secs = std::chrono::duration_cast<std::chrono::seconds>(time_now.time_since_epoch()).count();
        Rcpp::Rcout << parent_path << std::endl;
        std::string archive_name = "results_" + std::to_string(time_secs) + ".zip";
        zip(archive_name, output_dir);
    }
}

void FRMatrix::concatenate_results(std::string output_dir, std::string file_name_prefix, int stratum) {
    std::string outputFile = output_dir + "/" + file_name_prefix + "_stratum_" + std::to_string(stratum + 1) + ".tsv";

    std::ofstream out(outputFile, std::ios::binary);

    if (!out.is_open()) {
        std::cerr << "Failed to open the output file: " << outputFile << std::endl;
        return;
    }

    for (const auto& entry : fs::directory_iterator(output_dir)) {
        if (entry.path().extension() == ".tsv" && entry.path().filename().string().rfind(file_name_prefix, 0) == 0) {
            std::ifstream in(entry.path(), std::ios::binary);

            if (!in.is_open()) {
                std::cerr << "Failed to open " << entry.path() << std::endl;
                continue; 
            }

            out << in.rdbuf();  // Stream the file content directly
        }
    }

    out.close();
}