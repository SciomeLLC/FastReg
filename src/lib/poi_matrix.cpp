#include <poi_matrix.h>

FRMatrix POIMatrix::get_chunk(const std::vector<std::string> &rows,
                              const std::vector<std::string> &cols) {
  FRMatrix G = reader->read_chunk(rows, cols);
  if (squared) {
    G.data = arma::square(G.data);
  }
  return G;
}