#include "csv_reader.h"

/**
 * @brief Tokenizes a string into values based on the delimiter character.
 *
 * This function splits a line into tokens based on the delimiter (`delim_char`)
 * and removes any extraneous whitespace or carriage returns.
 *
 * @param line The string to tokenize.
 * @return A vector of strings representing the tokenized values.
 */
std::vector<std::string> CSVReader::tokenize(std::string &line)
{
    std::vector<std::string> tokens;
    line.erase(std::remove_if(line.begin(), line.end(),
                              [](char i)
                              { return (i == '\r'); }),
               line.end()); // remove the carriage return if it is there

    std::stringstream stream(line);
    std::string temp;
    // Loop over the stringstream until newline '\n' is hit
    while (std::getline(stream, temp, delim_char))
    {
        // Remove leading and trailing whitespace
        temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        temp.erase(std::find_if(temp.rbegin(), temp.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), temp.end());

        tokens.push_back(temp);
    }

    return tokens;
}
/**
 * @brief Reads the entire CSV file and stores the values in memory.
 *
 * This function reads the CSV file line by line, tokenizes each line,
 * and saves the result into the `values` member.
 */
void CSVReader::read_file()
{
  std::ifstream file;
  file.open(file_name);
  std::string in_line;
  std::vector<std::string> line;

  // go through the file line-by-line and tokenize it
  while (std::getline(file, in_line))
  {
    line = tokenize(in_line);
    values.push_back(line);
  }
  file.close();
}

/**
 * @brief Retrieves the headers from the CSV file.
 *
 * This function assumes the first row of the CSV file contains the headers.
 *
 * @return A vector of strings representing the headers.
 */
std::vector<std::string> CSVReader::get_headers() const
{
  if (values.size() == 0)
  {
    Rcpp::stop("Unable to read file: %s", file_name);
  }
  // return the values first row of the file
  return values.at(0);
}

/**
 * @brief Sets the row names using the values from the specified column.
 *
 * This function sets the row names by using the values from the column specified by `rowname_id`.
 * It assumes that the first row contains the headers.
 *
 * @param rowname_id The name of the column to use as row names.
 */
void CSVReader::set_row_names(const std::string &rowname_id)
{
  auto idx = std::find(headers.begin(), headers.end(), rowname_id);
  if (idx == headers.end())
  {
    Rcpp::stop("Column name %s not found in %s.", rowname_id, file_name);
  }
  int col_idx = idx - headers.begin();
  for (size_t i = 1; i < values.size(); i++)
  {
    row_names.push_back(values[i][col_idx]);
  }
}

/**
 * @brief Retrieves a column by its index.
 *
 * This function returns the values of the column specified by its index.
 * It skips the first row as it assumes it contains the headers.
 *
 * @param col_idx The index of the column to retrieve.
 * @return A vector of strings containing the values from the specified column.
 */
std::vector<std::string> CSVReader::get_col(const int col_idx) const
{
  std::vector<std::string> column;
  if (col_idx > (int)headers.size() || col_idx < 0)
  {
    Rcpp::stop("Column index %d is greater than the number of columns in %s.",
               col_idx, file_name);
  }
  for (size_t i = 1; i < values.size(); i++)
  {
    column.push_back(values[i][col_idx]);
  }

  return column;
}

/**
 * @brief Retrieves a column by its name.
 *
 * This function returns the values of the column specified by its name.
 * It skips the first row as it assumes it contains the headers.
 *
 * @param col_name The name of the column to retrieve.
 * @return A vector of strings containing the values from the specified column.
 */
std::vector<std::string> CSVReader::get_col(const std::string &col_name) const
{
  std::vector<std::string> column;
  auto idx = std::find(headers.begin(), headers.end(), col_name);
  if (idx == headers.end())
  {
    Rcpp::stop("Column name %s not found in %s.", col_name, file_name);
  }

  int col_idx = idx - headers.begin();
  for (size_t i = 1; i < values.size(); i++)
  {
    column.push_back(values[i][col_idx]);
  }

  return column;
}