#ifndef CSVREADER_H
#define CSVREADER_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

using namespace Rcpp;
/**
 * @brief The CSVReader class reads and processes CSV files.
 *
 * This class handles reading, tokenizing, and extracting columns and rows
 * from a CSV file. It supports delimiters such as commas, tabs, and semicolons.
 */
class CSVReader
{
public:
  std::string file_name;
  char delim_char = ',';
  std::vector<std::vector<std::string>> values;
  std::vector<std::string> row_names;
  std::vector<std::string> headers;
  CSVReader() {};
  /**
   * @brief Constructs a CSVReader with the specified filename and delimiter.
   *
   * @param filename The name of the CSV file to read.
   * @param delim The delimiter used in the file (e.g., "comma", "tab", "semicolon").
   * @param rowname_id The name of the column used as row names (default is "id").
   */
  CSVReader(std::string filename, std::string delim,
            std::string rowname_id = "id")
  {
    fs::path file_path = filename;
    if (!fs::exists(file_path))
    {
      Rcpp::stop("Covariate file: %s does not exist.", filename);
    }
    if (delim == "tab")
    {
      delim_char = '\t';
    }
    else if (delim == "comma")
    {
      delim_char = ',';
    }
    else if (delim == "semicolon")
    {
      delim_char = ';';
    }
    else
    {
      Rcpp::stop(
          "Invalid delim! delim for %s must be 'tab', 'comma', or 'semicolon'",
          filename);
    }

    file_name = std::move(filename);
    read_file();
    headers = get_headers();
    set_row_names(rowname_id);
  };
  /**
   * @brief Tokenizes a string into values based on the delimiter character.
   *
   * @param line The string to tokenize.
   * @return A vector of strings representing the tokenized values.
   */
  std::vector<std::string> tokenize(std::string &line);
  /**
   * @brief Retrieves the headers from the CSV file.
   *
   * This assumes that the first row of the file contains the headers.
   *
   * @return A vector of strings representing the headers of the CSV file.
   */
  std::vector<std::string> get_headers() const;
  /**
   * @brief Retrieves a column by its index.
   *
   * @param col_idx The index of the column to retrieve.
   * @return A vector of strings containing the values from the specified column.
   */
  std::vector<std::string> get_col(const int col_idx) const;
  /**
   * @brief Retrieves a column by its name.
   *
   * @param col_name The name of the column to retrieve.
   * @return A vector of strings containing the values from the specified column.
   */
  std::vector<std::string> get_col(const std::string &col_name) const;
  /**
   * @brief Sets the row names using the values from the specified column.
   *
   * @param rowname_id The name of the column to use as row names.
   */
  void set_row_names(const std::string &rowname_id);
  /**
   * @brief Reads the entire CSV file and stores the values in memory.
   */
  void read_file();
};
#endif