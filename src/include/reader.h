#ifndef READER_H
#define READER_H
#pragma once
#include "BEDMatrix.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
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

class Reader {
public:
  std::string file_name;
  virtual std::vector<std::string> get_individuals() = 0;
  virtual std::vector<std::string> get_names() = 0;
  virtual FRMatrix read_chunk(const std::vector<std::string> &rows,
                              const std::vector<std::string> &cols) = 0;
  virtual ~Reader(){};
};

class BEDReader : public Reader {
public:
  std::string base_filename;
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
    // Assume that the BED file is "filename.bed", and the BIM and FAM files are
    // "filename.bim" and "filename.fam"
    base_filename = filename.substr(0, filename.find_last_of("."));
    get_individuals();
    get_names();
  }
  ~BEDReader(){};
  std::vector<std::string> get_individuals() override;
  std::vector<std::string> get_names() override;
  FRMatrix read_chunk(const std::vector<std::string> &rows,
                      const std::vector<std::string> &cols) override;

private:
  float recode_genotype2(int genotype);
};

class POI : public Reader {
public:
  std::vector<std::string> names;
  std::vector<std::string> individuals;
  std::unordered_map<std::string, int> individuals_map;
  std::unordered_map<std::string, int> names_map;
  hid_t file_id = -1;
  hid_t values_dataset_id = -1;
  hid_t values_dataspace_id = -1;
  hid_t values_datatype = -1;
  H5T_class_t values_type_class;
  int rank = 2;
  bool transpose = false;
  hid_t memspace_id = -1;
  hsize_t hyperslab_dims[2];

  POI(const std::string &filename) { file_name = filename; }
  POI() {}
  ~POI() {
    H5garbage_collect();
    H5close();
  }
  void open(bool read_only = false);
  void get_values_dataset_id();
  void close_all();
  std::vector<std::string> get_individuals() override;
  std::vector<std::string> get_names() override;
  void get_data_type();
  void set_memspace(size_t rows, size_t cols);
  FRMatrix read_chunk(const std::vector<std::string> &rows,
                      const std::vector<std::string> &cols) override;
  void load_int_data_chunk(FRMatrix &G, hsize_t *hyperslab_dims,
                           hsize_t *src_offset);
  void load_float_data_chunk(FRMatrix &G, hsize_t *hyperslab_dims,
                             hsize_t *src_offset);
};
#endif