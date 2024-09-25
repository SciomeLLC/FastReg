#include "reader.h"

float BEDReader::recode_genotype2(int genotype) {
  float coding = arma::datum::nan; // missing
  if (genotype == 0) {
    coding = 2.0; // two copies of A1
  } else if (genotype == 3) {
    coding = 0.0; // zero copies of A1
  } else if (genotype == 2) {
    coding = 1.0; // one copy of A1
  }
  return coding;
}

std::vector<std::string> BEDReader::get_names() {
  if (headers.size() > 0) {
    return headers;
  }
  // Read BIM file to get num_snps and SNP IDs
  std::ifstream bim_file((base_filename + ".bim").c_str());
  if (!bim_file) {
    Rcpp::stop("Cannot open BIM file: %s", (base_filename + ".bim"));
  }
  std::string line;
  while (std::getline(bim_file, line)) {
    if (!line.empty()) {
      std::istringstream iss(line);
      std::string chrom, snp_id;
      iss >> chrom >> snp_id;
      headers.push_back(snp_id);
    }
  }
  bim_file.close();
  return headers;
}

std::vector<std::string> BEDReader::get_individuals() {
  if (row_names.size() > 0) {
    return row_names;
  }

  // Read FAM file to get num_samples and sample IDs
  std::ifstream fam_file((base_filename + ".fam").c_str());
  if (!fam_file) {
    Rcpp::stop("Cannot open FAM file: %s", (base_filename + ".fam"));
  }

  std::string line;
  while (std::getline(fam_file, line)) {
    if (!line.empty()) {
      std::istringstream iss(line);
      std::string fid, iid;
      iss >> fid >> iid;
      row_names.push_back(iid);
    }
  }
  fam_file.close();
  return row_names;
}

FRMatrix BEDReader::read_chunk(const std::vector<std::string> &rows,
                               const std::vector<std::string> &cols) {
  // Open the BED file in binary mode
  std::ifstream bed_file(file_name.c_str(), std::ios::in | std::ios::binary);
  if (!bed_file) {
    Rcpp::stop("Cannot open BED file: %s", file_name);
  }

  // Read the first three bytes (magic numbers and mode)
  unsigned char header[3];
  bed_file.read(reinterpret_cast<char *>(header), 3);
  if (bed_file.gcount() != 3) {
    Rcpp::stop("Failed to read the BED file header.");
  }

  // Check the magic numbers
  if (header[0] != 0x6C || header[1] != 0x1B) {
    Rcpp::stop("Invalid BED file magic numbers.");
  }

  // Check the mode
  unsigned char mode = header[2];
  if (mode != 0x01) {
    Rcpp::stop("Unsupported BED file mode. Only SNP-major mode (mode=1) is "
               "supported.");
  }

  // Determine the number of samples (individuals) and SNPs (variants)
  int num_samples = row_names.size();

  // Compute number of bytes per SNP
  int num_bytes_per_snp = (num_samples + 3) / 4;

  std::vector<int> row_ids;
  for (size_t i = 0; i < rows.size(); i++) {
    auto idx = std::find(row_names.begin(), row_names.end(), rows[i]);
    if (idx == row_names.end()) {
      Rcpp::Rcout << "Cannot find: " << rows[i] << " in BED file" << std::endl;
      continue;
    }
    row_ids.push_back(idx - row_names.begin());
  }

  std::vector<int> col_ids;
  for (size_t i = 0; i < cols.size(); i++) {
    auto idx = std::find(headers.begin(), headers.end(), cols[i]);
    if (idx == headers.end()) {
      Rcpp::Rcout << "Cannot find: " << cols[i] << " in BED file" << std::endl;
      continue;
    }
    col_ids.push_back(idx - headers.begin());
  }

  int ni = row_ids.size();
  int nj = col_ids.size();

  FRMatrix snps_mat;
  snps_mat.data = arma::fmat(ni, nj, arma::fill::zeros);
  for (int cj = 0; cj < nj; cj++) {
    int snp_idx = col_ids[cj];

    // Compute the file offset for this SNP
    std::streampos offset =
        3 + static_cast<std::streampos>(snp_idx) * num_bytes_per_snp;

    // Seek to the offset
    bed_file.seekg(offset);
    if (!bed_file) {
      Rcpp::stop("Failed to seek to SNP position in BED file.");
    }

    // Read the genotype data for this SNP
    std::vector<unsigned char> genotype_bytes(num_bytes_per_snp);
    bed_file.read(reinterpret_cast<char *>(&genotype_bytes[0]),
                  num_bytes_per_snp);
    if (bed_file.gcount() != num_bytes_per_snp) {
      Rcpp::stop("Failed to read genotype data for SNP %d.", snp_idx + 1);
    }

    arma::fcolvec geno_vec(ni);
    // For each sample in i
    for (int ci = 0; ci < ni; ci++) {
      int sample_idx = row_ids[ci];
      // Compute the genotype for this sample
      int byte_idx = sample_idx / 4;
      int within_byte_idx = sample_idx % 4;

      unsigned char byte = genotype_bytes[byte_idx];

      // Extract the two bits corresponding to the genotype
      int genotype_code = (byte >> (2 * within_byte_idx)) & 0x03;

      // Recode the genotype
      float genotype = recode_genotype2(genotype_code);
      geno_vec[ci] = genotype;
    }

    std::string snp_id = headers[snp_idx];
    snps_mat.data.col(cj) = geno_vec;
    snps_mat.col_names[snp_id] = cj;
    snps_mat.col_names_arr.push_back(snp_id);
  }

  snps_mat.row_names_arr = row_names;
  for (size_t i = 0; i < row_names.size(); i++) {
    snps_mat.row_names[row_names[i]] = i;
  }
  return snps_mat;
}