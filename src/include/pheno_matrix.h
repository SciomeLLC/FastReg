#ifndef PHENO_MATRIX_H
#define PHENO_MATRIX_H
#pragma once
#include <csv_reader.h>
#include <fr_matrix.h>
/**
 * @brief The PhenoMatrix class handles creating a matrix from phenotype data.
 *
 * This class reads a CSV file containing phenotype data, processes it into an
 * FRMatrix object, and provides methods to manipulate the data, such as squaring the values.
 */
class PhenoMatrix
{
public:
  std::string file_name;
  std::string id_col;
  std::string phenotype;
  int id_col_idx;
  CSVReader reader;
  FRMatrix mat;
  std::vector<std::string> headers;
  /**
   * @brief Constructs a PhenoMatrix object.
   *
   * This constructor initializes the PhenoMatrix by reading the CSV file and setting
   * the ID and phenotype columns.
   *
   * @param filename The name of the phenotype file.
   * @param delim The delimiter used in the file (e.g., comma, tab).
   * @param id The name of the ID column in the file.
   * @param phenotype The name of the phenotype column in the file.
   */
  PhenoMatrix(std::string &filename, std::string &delim, std::string &id,
              std::string &phenotype)
  {
    file_name = filename;
    reader = CSVReader(filename, delim, id);
    headers = reader.get_headers();
    id_col = id;
    this->phenotype = phenotype;
  }
  /**
   * @brief Creates the phenotype matrix.
   *
   * This method generates a matrix from the phenotype data, populating the FRMatrix object.
   *
   * @return The FRMatrix containing the phenotype data.
   */
  FRMatrix create_matrix();
};
#endif