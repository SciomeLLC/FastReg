#ifndef COV_MATRIX_H
#define COV_MATRIX_H
#pragma once
#include <fr_matrix.h>
#include <csv_reader.h>
#include <covariate.h>

class CovariateMatrix
{
public:
  std::string file_name;
  std::string id_col;
  int id_cold_idx;
  double colinearity_rsq;
  bool has_intercept;
  std::vector<Covariate> covariates;
  std::vector<std::string> headers;
  std::vector<int> covariate_cols_idx;
  std::vector<std::string> covariate_col_names;
  CSVReader reader;
  FRMatrix design_matrix;

  CovariateMatrix(std::string &filename, std::string &delim, std::string &id, std::vector<Covariate> covariates, double colinearity_rsq, bool no_intercept)
  {
    file_name = filename;
    reader = CSVReader(filename, delim);
    this->covariates = covariates;
    headers = reader.get_headers();
    id_col = id;
    this->colinearity_rsq = colinearity_rsq;
    has_intercept = !no_intercept;
    set_cov_col_idx();
  };

  FRMatrix create_design_matrix();

private:
  std::vector<std::string> get_col_names(std::vector<int> idx);
  int find_col_idx(std::string col_name);
  void set_cov_col_idx();
};
#endif