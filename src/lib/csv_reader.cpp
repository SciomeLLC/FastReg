#include "csv_reader.h"

////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief tokenize - split a string into tokens assuming based upon the
/// 'delim_char' character
/// @param line string to split
/// @return returns a standard vector of strings where each one is represents
/// a value between the tokens
std::vector<std::string> CSVReader::tokenize(std::string &line) {
  std::vector<std::string> tokens(0);
  line.erase(std::remove_if(line.begin(), line.end(),
                            [](char i) { return (i == '\r'); }),
             line.end()); // remove the carriage return if it is there
  std::stringstream stream(line);
  std::string temp;
  // Loop over the stringstream until newline '\n' is hit
  while (!stream.eof()) {
    std::getline(stream, temp, delim_char);
    temp.erase(std::remove_if(temp.begin(), temp.end(), ::isspace),
               temp.end()); // remove all whitespace from the value
    tokens.push_back(temp);
  }

  return tokens;
}

/// @brief read_file() - reads the file line by line and tokenizes each line.
///                      The result is saved into 'values'.
/// @return - none
void CSVReader::read_file() {
  std::ifstream file;
  file.open(file_name);
  std::string in_line;
  std::vector<std::string> line;

  // go through the file line-by-line and tokenize it
  while (std::getline(file, in_line)) {
    line = tokenize(in_line);
    values.push_back(line);
  }
  file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief get_headers - returns the first line of the file. This assumes that
/// the file has headers
/// @return returns a standard vector of strings where each element is the
/// name of the column at that index
std::vector<std::string> CSVReader::get_headers() const {
  if (values.size() == 0) {
    Rcpp::stop("Unable to read file: %s", file_name);
  }
  // return the values first row of the file
  return values.at(0);
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set_row_names - loads the row_names field with values from the column
/// specified by 'col_idx'. This skips the first row as it assumes that it
/// contains the header names.
/// @return none
void CSVReader::set_row_names(const std::string& rowname_id) {
  auto idx = std::find(headers.begin(), headers.end(), rowname_id);
  if (idx == headers.end()) {
    Rcpp::stop("Column name %s not found in %s.", rowname_id, file_name);
  }
  int col_idx = idx - headers.begin();
  for (size_t i = 1; i < values.size(); i++) {
    row_names.push_back(values[i][col_idx]);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief get_col - gets the column values specified by the column index -
/// 'col_idx'. This skips the first row as it assumes that it contains the
/// header names.
/// @param col_idx - index of the column to retrieve
/// @return returns a standard vector of strings containing the values from the
/// column specified by 'col_idx'
std::vector<std::string> CSVReader::get_col(const int col_idx) const {
  std::vector<std::string> column;
  if (col_idx > (int)headers.size() || col_idx < 0) {
    Rcpp::stop("Column index %d is greater than the number of columns in %s.",
               col_idx, file_name);
  }
  for (size_t i = 1; i < values.size(); i++) {
    column.push_back(values[i][col_idx]);
  }

  return column;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief get_col - gets the column values specified by the column index -
/// 'col_idx'. This skips the first row as it assumes that it contains the
/// headers.
/// @param col_idx - index of the column to retrieve
/// @return returns a standard vector of strings containing the values from the
/// column specified by 'col_idx'
std::vector<std::string> CSVReader::get_col(const std::string& col_name) const {
  std::vector<std::string> column;
  auto idx = std::find(headers.begin(), headers.end(), col_name);
  if (idx == headers.end()) {
    Rcpp::stop("Column name %s not found in %s.", col_name, file_name);
  }

  int col_idx = idx - headers.begin();
  for (size_t i = 1; i < values.size(); i++) {
    column.push_back(values[i][col_idx]);
  }

  return column;
}