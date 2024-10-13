#ifndef READER_H
#define READER_H
#pragma once
// #include "BEDMatrix.h"
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

/**
 * @brief The Reader class is an abstract base class for reading data from various file types.
 *
 * This class defines a common interface for accessing datasets such as individuals and names
 * (predictors of interest), and reading chunks of data. It is designed to be extended by specific
 * implementations, like HDF5-based readers.
 */
class Reader
{
public:
  std::string file_name;
  /**
   * @brief Retrieves a list of individuals from the dataset.
   *
   * This function is pure virtual and must be implemented by any class inheriting from Reader.
   *
   * @return A vector of strings representing individual identifiers.
   */
  virtual std::vector<std::string> get_individuals() = 0;
  /**
   * @brief Retrieves the names (predictors of interest) from the dataset.
   *
   * This function is pure virtual and must be implemented by any class inheriting from Reader.
   *
   * @return A vector of strings representing the names of the predictors.
   */
  virtual std::vector<std::string> get_names() = 0;
  /**
   * @brief Reads a chunk of data from the dataset for specified rows and columns.
   *
   * This function is pure virtual and must be implemented by any class inheriting from Reader.
   *
   * @param rows A vector of row identifiers to read.
   * @param cols A vector of column identifiers to read.
   * @return An FRMatrix containing the requested chunk of data.
   */
  virtual FRMatrix read_chunk(const std::vector<std::string> &rows,
                              const std::vector<std::string> &cols) = 0;
  /**
   * @brief Virtual destructor to ensure proper cleanup.
   */
  virtual ~Reader() {};
};

/**
 * @brief The POI class provides an interface to read predictors of interest (POI) from HDF5 files.
 *
 * This class extends the Reader class and provides the functionality to read individuals, names,
 * and chunks of data from HDF5 files. It handles various HDF5 types and ensures memory is managed correctly.
 */
class POI : public Reader
{
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
  /**
   * @brief Constructs a POI object with a given filename.
   *
   * @param filename Name of the HDF5 file to open.
   */
  POI(const std::string &filename) { file_name = filename; }
  POI() {}
  ~POI()
  {
    H5garbage_collect();
    H5close();
  }
  /**
   * @brief Opens the HDF5 file.
   *
   * Opens the HDF5 file with either read-only or write permissions.
   *
   * @param read_only Flag to indicate if the file should be opened in read-only mode.
   */
  void open(bool read_only = false);
  /**
   * @brief Retrieves the dataset identifier for the values dataset.
   */
  void get_values_dataset_id();
  /**
   * @brief Closes all open HDF5 objects and file handles.
   */
  void close_all();
  /**
   * @brief Retrieves the list of individuals from the HDF5 file.
   *
   * @return A vector of strings representing individual identifiers.
   */
  std::vector<std::string> get_individuals() override;
  /**
   * @brief Retrieves the list of names (predictors) from the HDF5 file.
   *
   * @return A vector of strings representing predictor names.
   */
  std::vector<std::string> get_names() override;
  /**
   * @brief Retrieves the data type of the values dataset in the HDF5 file.
   */
  void get_data_type();
  /**
   * @brief Sets the memory space for reading a chunk of data.
   *
   * Defines the size of the memory space to be used when reading chunks of data from the file.
   *
   * @param rows Number of rows to read.
   * @param cols Number of columns to read.
   */
  void set_memspace(size_t rows, size_t cols);
  /**
   * @brief Reads a chunk of data from the HDF5 file for specified rows and columns.
   *
   * @param rows A vector of row identifiers to read.
   * @param cols A vector of column identifiers to read.
   * @return An FRMatrix containing the requested chunk of data.
   */
  FRMatrix read_chunk(const std::vector<std::string> &rows,
                      const std::vector<std::string> &cols) override;
  /**
   * @brief Loads an integer chunk of data from the HDF5 dataset.
   *
   * @param G The matrix to store the data.
   * @param hyperslab_dims The dimensions of the hyperslab.
   * @param src_offset The offset for the data in the HDF5 file.
   */
  void load_int_data_chunk(FRMatrix &G, hsize_t *hyperslab_dims,
                           hsize_t *src_offset);
  /**
   * @brief Loads a floating-point chunk of data from the HDF5 dataset.
   *
   * @param G The matrix to store the data.
   * @param hyperslab_dims The dimensions of the hyperslab.
   * @param src_offset The offset for the data in the HDF5 file.
   */
  void load_float_data_chunk(FRMatrix &G, hsize_t *hyperslab_dims,
                             hsize_t *src_offset);
};
#endif