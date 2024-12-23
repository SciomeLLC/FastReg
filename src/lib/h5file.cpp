#include <reader.h>

std::vector<std::string> POI::get_individuals() {
  const char dname[] = "individuals";
  hid_t ind_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

  if (ind_dataset < 0) {
    Rcpp::stop("Failed to open individuals dataset.");
  }

  hid_t datatype = H5Dget_type(ind_dataset);
  hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

  if (H5Tget_class(native_type) != H5T_STRING) {
    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(ind_dataset);
    Rcpp::stop("Dataset does not have the expected string data type.");
  }
  if (H5Tget_class(native_type) != H5T_STRING) {
    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(ind_dataset);
    Rcpp::stop(
        "Dataset does not have the expected variable-length string data type.");
  }

  if (H5Tis_variable_str(native_type)) {
    hid_t space = H5Dget_space(ind_dataset);
    hsize_t num_ind;
    H5Sget_simple_extent_dims(space, &num_ind, NULL);

    char **rdata = new char *[num_ind];

    if (H5Dread(ind_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                rdata) < 0) {
      H5Tclose(native_type);
      H5Tclose(datatype);
      H5Dclose(ind_dataset);
      delete[] rdata;
      Rcpp::stop("Failed to read individuals dataset");
    }

    for (hsize_t i = 0; i < num_ind; i++) {
      std::string name = rdata[i];
      individuals.push_back(name);
      individuals_map[name] = i;
      // Freeing the variable-length data
      free(rdata[i]);
    }

    delete[] rdata;
    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(ind_dataset);
    H5Sclose(space);
  } else {
    size_t datatype_size = H5Tget_size(native_type);
    hid_t space = H5Dget_space(ind_dataset);
    hsize_t num_ind;
    H5Sget_simple_extent_dims(space, &num_ind, NULL);
    individuals.reserve(num_ind);

    char *buffer = new char[num_ind * datatype_size];

    if (H5Dread(ind_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                buffer) < 0) {
      delete[] buffer;
      H5Sclose(space);
      H5Tclose(native_type);
      H5Tclose(datatype);
      H5Dclose(ind_dataset);
      Rcpp::stop("Failed to read individuals dataset");
    }
    for (hsize_t i = 0; i < num_ind; i++) {
      std::string name = std::string(&buffer[i * datatype_size], datatype_size);
      // Remove potential trailing white spaces due to fixed-length string
      // format
      name.erase(name.find_last_not_of(' ') + 1);
      individuals.push_back(name);
      individuals_map[name] = i;
    }

    delete[] buffer;
    H5Sclose(space);
    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(ind_dataset);
  }
  return individuals;
}

std::vector<std::string> POI::get_names() {
  const char dname[] = "predictors_of_interest";
  hid_t poi_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

  if (poi_dataset < 0) {
    Rcpp::stop("Failed to open predictors_of_interest dataset.");
  }

  hid_t datatype = H5Dget_type(poi_dataset);
  hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

  if (H5Tget_class(native_type) != H5T_STRING) {
    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(poi_dataset);
    Rcpp::stop("Dataset does not have the expected string data type.");
  }

  if (H5Tis_variable_str(native_type)) {
    hid_t space = H5Dget_space(poi_dataset);
    hsize_t num_poi;
    H5Sget_simple_extent_dims(space, &num_poi, NULL);
    char **rdata = new char *[num_poi];

    if (H5Dread(poi_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                rdata) < 0) {
      delete[] rdata;
      H5Tclose(native_type);
      H5Tclose(datatype);
      H5Dclose(poi_dataset);
      H5Sclose(space);
      Rcpp::stop("Failed to read individuals dataset");
    }

    for (hsize_t i = 0; i < num_poi; i++) {
      std::string name = rdata[i];
      names.push_back(name);
      names_map[name] = i;
      // Freeing the variable-length data
      free(rdata[i]);
    }

    delete[] rdata;
    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(poi_dataset);
    H5Sclose(space);
  } else {
    size_t datatype_size = H5Tget_size(native_type);
    hid_t space = H5Dget_space(poi_dataset);
    hsize_t num_poi;
    H5Sget_simple_extent_dims(space, &num_poi, NULL);
    names.reserve(num_poi);

    char *buffer = new char[num_poi * datatype_size];
    if (H5Dread(poi_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                buffer) < 0) {
      delete[] buffer;
      H5Sclose(space);
      H5Tclose(native_type);
      H5Tclose(datatype);
      H5Dclose(poi_dataset);
      Rcpp::stop("Failed to read POI dataset");
    }

    for (hsize_t i = 0; i < num_poi; i++) {
      std::string name = std::string(&buffer[i * datatype_size], datatype_size);
      // Remove potential trailing white spaces due to fixed-length string
      // format
      name.erase(name.find_last_not_of(' ') + 1);
      names.push_back(name);
      names_map[name] = i;
    }

    H5Tclose(native_type);
    H5Tclose(datatype);
    H5Dclose(poi_dataset);
    delete[] buffer;
    H5Sclose(space);
  }

  return names;
}

void POI::set_memspace(size_t rows, size_t cols) {
  hyperslab_dims[0] = rows;
  hyperslab_dims[1] = cols;
  if (memspace_id > -1) {
    H5Sclose(memspace_id);
    memspace_id = -1;
  }
  memspace_id = H5Screate_simple(rank, hyperslab_dims, NULL);
}

void POI::get_data_type() {

  // Get the dataspace
  values_dataspace_id = H5Dget_space(values_dataset_id);
  values_datatype = H5Dget_type(values_dataset_id);
  values_type_class = H5Tget_class(values_datatype);
}

FRMatrix POI::read_chunk(const std::vector<std::string> &rows,
                         const std::vector<std::string> &cols) {
  FRMatrix G;
  G.row_names.reserve(rows.size());
  G.col_names.reserve(cols.size());
  G.row_names = std::unordered_map<std::string, int>(rows.size());
  G.col_names = std::unordered_map<std::string, int>(cols.size());
  std::vector<hsize_t> row_indices(rows.size());
  std::vector<hsize_t> col_indices(cols.size());

  for (size_t i = 0; i < rows.size(); i++) {
    G.row_names[rows[i]] = i;
    row_indices[i] = individuals_map[rows[i]];
  }
  for (size_t i = 0; i < cols.size(); i++) {
    G.col_names[cols[i]] = i;
    col_indices[i] = names_map[cols[i]];
  }

  // Get dimensions of the dataspace
  values_dataspace_id = H5Dget_space(values_dataset_id);
  values_datatype = H5Dget_type(values_dataset_id);
  values_type_class = H5Tget_class(values_datatype);
  int ndims = H5Sget_simple_extent_ndims(values_dataspace_id);
  std::vector<hsize_t> dims(ndims);
  H5Sget_simple_extent_dims(values_dataspace_id, dims.data(), NULL);
  if (individuals.size() != dims[0] && individuals.size() == dims[1]) {
    transpose = true;
  }

  if (individuals.size() != dims[0] && individuals.size() != dims[1]) {
    close_all();
    Rcpp::stop(
        "Dimensions of the dataset do not match the sizes of "
        "poi_individuals. Please check the hdf5 file dataset dimensions.");
  }

  hsize_t hyperslab_dims[2] = {rows.size(), cols.size()};
  // Define the hyperslab for the entire range of columns needed
  hsize_t src_offset[2] = {0, col_indices[0]};
  if (transpose) {
    // Swap hyperslab dimensions and src offsets
    std::swap(hyperslab_dims[0], hyperslab_dims[1]);
    std::swap(src_offset[0], src_offset[1]);
  }

  if (values_type_class == H5T_FLOAT) {
    load_float_data_chunk(G, hyperslab_dims, src_offset);
  } else {
    close_all();
    Rcpp::stop("HDF5 dataset type class is not float or int");
  }
  return G;
}

void POI::load_float_data_chunk(FRMatrix &G, hsize_t *hyperslab_dims,
                                hsize_t *src_offset) {
  hsize_t memspace_dims[2] = {hyperslab_dims[0], hyperslab_dims[1]};
  G.data.resize(memspace_dims[1], memspace_dims[0]).fill(0.0);
  hsize_t dst_offset[2] = {0, 0};
  // n x m
  herr_t status = H5Sselect_hyperslab(values_dataspace_id, H5S_SELECT_SET,
                                      src_offset, NULL, hyperslab_dims, NULL);
  if (status < 0) {
    Rcpp::stop("Error selecting hyperslab in file space.");
  }

  // Create a memory dataspace
  memspace_id = H5Screate_simple(2, memspace_dims, NULL);
  status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, dst_offset, NULL,
                               memspace_dims, NULL);
  if (status < 0) {
    Rcpp::stop("Error selecting hyperslab in memory space.");
  }
  // Reading 64-bit double
  status = H5Dread(values_dataset_id, H5T_NATIVE_FLOAT, memspace_id,
                   values_dataspace_id, H5P_DEFAULT, G.data.memptr());
  if (status < 0) {
    Rcpp::stop("Error reading float data from values dataset.");
  }

  if (!transpose) {
    arma::inplace_trans(G.data);
  }
#if defined(_DEBUG)
  Rcpp::Rcout << "Post read G data size: " << G.data.n_cols << "x"
              << G.data.n_rows << std::endl;
#endif
  if (memspace_id > 0) {
    H5Sclose(memspace_id);
    memspace_id = -1;
  }
}
void POI::close_all() {
  if (memspace_id > 0) {
    H5Sclose(memspace_id);
    memspace_id = -1;
  }

  if (values_dataspace_id > 0) {
    H5Sclose(values_dataspace_id);
    values_dataspace_id = -1;
  }

  if (values_dataset_id > 0) {
    H5Dclose(values_dataset_id);
    values_dataset_id = -1;
  }

  if (values_datatype > 0) {
    H5Tclose(values_datatype);
    values_datatype = -1;
  }

  if (file_id > 0) {
    H5Fclose(file_id);
    file_id = -1;
  }
}

void POI::open(bool read_only) {
  file_id = H5Fopen(file_name.c_str(), H5F_ACC_SWMR_READ | H5F_ACC_RDONLY,
                    H5P_DEFAULT);
  if (file_id < 0) {
    Rcpp::stop("Failed to open HDF5 file.");
  }
}

void POI::get_values_dataset_id() {
  values_dataset_id = H5Dopen(file_id, "values", H5P_DEFAULT);
}
