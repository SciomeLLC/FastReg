#ifndef CSVREADER_H
#define CSVREADER_H
#pragma once
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <list>
#include <regex>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
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

#include <utils.h>
using namespace arma;

class CSVReader
{
public:
  std::string file_name;
  char delim_char;
  std::vector<std::vector<std::string>> tokenized;
  std::vector<std::string> row_names;
  CSVReader() {};
  CSVReader(std::string filename, std::string delim)
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

    file_name = filename;
    read_file();
    set_row_names();
  };
  std::vector<std::string> tokenize(std::string line);
  std::vector<std::string> get_headers();
  void set_row_names();
  void read_file();
};
#endif