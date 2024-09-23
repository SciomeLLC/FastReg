#ifndef BEDREADER_H
#define BEDREADER_H
#pragma once
#include "BEDMatrix.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fr_matrix.h>
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
class BEDReader {
public:
  std::string file_name, base_filename;
  arma::fmat mat;
  std::vector<std::string> row_names;
  std::vector<std::string> headers;
  std::vector<int> row_idx, col_idx;
  BEDReader(){};
  BEDReader(std::string &filename) {
    fs::path file_path = filename;
    if (!fs::exists(file_path)) {
      Rcpp::stop("BED file: %s does not exist.", filename);
    }

    file_name = filename;
    // Assume that the BED file is "filename.bed", and the BIM and FAM files are "filename.bim" and "filename.fam"
    base_filename = filename.substr(0, filename.find_last_of("."));
    get_individuals();
    get_names();
  }
  std::vector<std::string> get_individuals();
  std::vector<std::string> get_names();
  FRMatrix read_chunk(std::vector<std::string> &rows,
                           std::vector<std::string> &cols);

  private:
  float recode_genotype2(int genotype);
};
#endif