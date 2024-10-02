#ifndef POIMATRIX_H
#define POIMATRIX_H
#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <fr_matrix.h>
#include <fstream>
#include <iostream>
#include <reader.h>
#include <sstream>
#include <string>
#include <utils.h>
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
  arma::fmat data;
  double maf_threshold, hwe_threshold;
  std::string effect_type;
  std::string poi_type;
  Reader *reader;
  POIMatrix(){};
  POIMatrix(Reader *rdr, double maf, double hwe, std::string poi_effect_type,
            std::string poi_type) {
    reader = rdr;
    names = reader->get_names();
    individuals = reader->get_individuals();
    maf_threshold = maf;
    hwe_threshold = hwe;
    effect_type = poi_effect_type;
    this->poi_type = poi_type;
  }
  FRMatrix get_chunk(const std::vector<std::string> &rows,
                     const std::vector<std::string> &cols);
  void filter_rows(FRMatrix &chunk, const std::vector<std::string> r_names);
  FRMatrix filter_genotype(FRMatrix &poi_matrix);
  void set_mask();

private:
  void check_threshold(FRMatrix &res, FRMatrix &poi_matrix);
  void check_dominant_recessive(FRMatrix &poi_matrix);
};
#endif