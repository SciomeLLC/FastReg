#ifndef COV_MATRIX_H
#define COV_MATRIX_H
#pragma once
#include <fr_matrix.h>
#include <csv_reader.h>
#include <covariate.h>

/**
 * @brief Manages covariates and creates a design matrix for statistical analysis.
 *
 * The CovariateMatrix class reads a CSV file, manages covariates, and generates a design matrix
 * based on the provided covariates. It handles cases with or without intercepts and performs
 * colinearity filtering on the covariates.
 */
class CovariateMatrix
{
public:
  std::string file_name;
  std::string id_col;
  int id_col_idx;
  double colinearity_rsq;
  bool has_intercept;
  std::vector<Covariate> covariates;
  std::vector<std::string> headers;
  std::vector<int> covariate_cols_idx;
  std::vector<std::string> covariate_col_names;
  CSVReader reader;
  FRMatrix design_matrix;

  /**
   * @brief Constructs a CovariateMatrix object.
   *
   * This constructor initializes the CovariateMatrix by reading the CSV file and setting
   * the covariates, headers, and intercept information.
   *
   * @param filename Name of the CSV file containing covariate data.
   * @param delim The delimiter used in the CSV file (e.g., comma, tab).
   * @param id The name of the ID column in the CSV file.
   * @param covariates A vector of Covariate objects to be included in the matrix.
   * @param colinearity_rsq R-squared threshold for colinearity filtering.
   * @param no_intercept If true, the design matrix will not include an intercept column.
   */
  CovariateMatrix(std::string &filename, std::string &delim, std::string &id, std::vector<Covariate> covariates, double colinearity_rsq, bool no_intercept)
  {
    file_name = filename;
    reader = CSVReader(filename, delim, id);
    this->covariates = covariates;
    headers = reader.get_headers();
    id_col = id;
    this->colinearity_rsq = colinearity_rsq;
    has_intercept = !no_intercept;
    set_cov_col_idx();
  };

  FRMatrix create_design_matrix();

private:
  /**
   * @brief Retrieves column names based on a vector of indices.
   *
   * Given a list of column indices, this function returns the corresponding column names
   * from the CSV file headers.
   *
   * @param idx A vector of column indices.
   * @return A vector of strings representing the column names.
   */
  std::vector<std::string> get_col_names(std::vector<int> idx);
    /**
   * @brief Finds the index of a given column name.
   *
   * Searches for the specified column name in the headers and returns its index.
   * If the column name is not found, a warning is displayed.
   *
   * @param col_name The name of the column to find.
   * @return The index of the column, or -1 if not found.
   */
  int find_col_idx(std::string col_name);
    /**
   * @brief Sets the column indices for the covariates.
   *
   * This function searches the CSV file headers for the covariate names and
   * sets their corresponding column indices. If a covariate name is not found,
   * the function throws an error.
   */
  void set_cov_col_idx();
};
#endif