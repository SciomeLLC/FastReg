#ifndef PHENO_MATRIX_H
#define PHENO_MATRIX_H
#pragma once
#include <csv_reader.h>
#include <fr_matrix.h>
#include <utils.h>
class PhenoMatrix {
public:
  std::string file_name;
  std::string id_col;
  std::string phenotype;
  int id_col_idx;
  CSVReader reader;
  FRMatrix mat;
  std::vector<std::string> headers;

  PhenoMatrix(std::string &filename, std::string &delim, std::string &id,
              std::string &phenotype) {
    file_name = filename;
    reader = CSVReader(filename, delim, id);
    headers = reader.get_headers();
    id_col = id;
    this->phenotype = phenotype;
  }
  FRMatrix create_matrix();
  FRMatrix get_squared();
};
#endif