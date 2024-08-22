#include "csv_reader.h"

////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief tokenize - split a string into tokens assuming based upon the 'delim_char' character
/// @param line   -- string to split
/// @return -- returns a standard vector of strings where each one is represents a value between the tokens
std::vector<std::string> CSVReader::tokenize(std::string line)
{
  std::vector<std::string> tokens(0);
  line.erase(std::remove_if(line.begin(), line.end(),
                            [](char i)
                            { return (i == '\r'); }),
             line.end()); // remove the carriage return if it is there
  std::stringstream stream(line);
  std::string temp;
  // Loop over the stringstream until newline '\n' is hit
  while (!stream.eof())
  {
    std::getline(stream, temp, delim_char);
    temp.erase(std::remove_if(temp.begin(), temp.end(), ::isspace), temp.end()); // remove all whitespace from the value
    tokens.push_back(temp);
  }

  return tokens;
}

/// @brief read_file() - reads the file line by line and tokenizes each line.
///                      The result is saved into 'tokenized'.
/// @return - none
void CSVReader::read_file()
{
  std::ifstream file;
  file.open(file_name);
  std::string in_line;
  std::vector<std::string> line;                        // each line of the file
  std::vector<std::vector<std::string>> tokenized_file; // a vector of vectors representing the tokenized lines of the file

  // go through the file line-by-line and tokenize it
  while (std::getline(file, in_line))
  {
    line = tokenize(in_line);
    tokenized.push_back(line);
  }
  file.close();
}

std::vector<std::string> CSVReader::get_headers()
{
  if (tokenized.size() == 0)
  {
    Rcpp::stop("Unable to read file: %s", file_name);
  }
  // return the tokenized first row of the file
  return tokenized.at(0);
}

void CSVReader::set_row_names()
{
  for(size_t i = 1; i < tokenized.size(); i++) {
    row_names.push_back(tokenized[i][0]);
  }
}