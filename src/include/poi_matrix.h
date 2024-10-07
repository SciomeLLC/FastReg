#ifndef POIMATRIX_H
#define POIMATRIX_H
#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <reader.h>
#include <fr_matrix.h>
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

class POIMatrix {
public:
  std::vector<std::string> names;
  std::vector<std::string> individuals;
  bool squared = false;
  arma::fmat data;
  Reader *reader;
  POIMatrix(){};
  POIMatrix(Reader *rdr, bool sqrd = false) {
    reader = rdr;
    names = reader->get_names();
    individuals = reader->get_individuals();
    squared = sqrd;
  }
  FRMatrix get_chunk(const std::vector<std::string>& rows, const std::vector<std::string>& cols);
};
#endif