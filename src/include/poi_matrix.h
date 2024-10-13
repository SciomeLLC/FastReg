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
/**
 * @brief The POIMatrix class handles the matrix representation of points of interest (POI).
 *
 * This class uses a Reader object to read and manipulate chunks of POI data, with the option
 * to return the squared values of the matrix if specified.
 */
class POIMatrix
{
public:
  std::vector<std::string> names;
  std::vector<std::string> individuals;
  bool squared = false;
  arma::fmat data;
  Reader *reader;
  POIMatrix() {};
  /**
   * @brief Constructs a POIMatrix object using a Reader.
   *
   * Initializes the POIMatrix by reading the names and individuals from the provided Reader object.
   * Optionally, the matrix data can be squared.
   *
   * @param rdr Pointer to the Reader object used to read the data.
   * @param sqrd Boolean flag to indicate if the data should be squared (default: false).
   */
  POIMatrix(Reader *rdr, bool sqrd = false)
  {
    reader = rdr;
    names = reader->get_names();
    individuals = reader->get_individuals();
    squared = sqrd;
  }
  /**
   * @brief Retrieves a chunk of the POI matrix for specific rows and columns.
   *
   * This method reads a subset of the POI matrix corresponding to the given rows and columns.
   * If the `squared` flag is set, the values in the matrix will be squared.
   *
   * @param rows Vector of strings representing the row names to retrieve.
   * @param cols Vector of strings representing the column names to retrieve.
   * @return An FRMatrix object containing the requested chunk of data.
   */
  FRMatrix get_chunk(const std::vector<std::string> &rows, const std::vector<std::string> &cols);
};
#endif