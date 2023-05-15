// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "hdf5.h"
#include <fr_matrix.h>
#include <omp.h>
#include <h5file.h>

void H5File::get_POI_individuals() {
    const char dname[] = "individuals";
    hid_t ind_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

    if (ind_dataset < 0) {
        Rcpp::Rcerr << "Failed to open individuals dataset." << std::endl;
    }

    hid_t datatype = H5Dget_type(ind_dataset);
    hid_t space = H5Dget_space(ind_dataset);
    hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

    hsize_t num_ind;
    H5Sget_simple_extent_dims(space, &num_ind, NULL);

    individuals.reserve(num_ind);

    std::vector<char*> ind_names(num_ind);
    H5Dread(ind_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ind_names.data());

    for (unsigned int i = 0; i < num_ind; i++) {
        std::string name = std::string(ind_names[i], ind_names[i] + strlen(ind_names[i]));
        individuals.push_back(name);
        individuals_map[name] = i;
        free(ind_names[i]);
    }
    H5Tclose(datatype);
    H5Sclose(space);
    H5Dclose(ind_dataset);
}

void H5File::get_POI_names() {
    const char dname[] = "predictors_of_interest";
    hid_t poi_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

    if (poi_dataset < 0) {
        Rcpp::Rcerr << "Failed to open predictors_of_interest dataset." << std::endl;
        return;
    }

    hid_t datatype = H5Dget_type(poi_dataset);
    hid_t space = H5Dget_space(poi_dataset);
    hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

    hsize_t num_ind;
    H5Sget_simple_extent_dims(space, &num_ind, NULL);
    names.reserve(num_ind);

    std::vector<char*> poi_names(num_ind);
    H5Dread(poi_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, poi_names.data());

    for (unsigned int i = 0; i < num_ind; i++) {
        std::string poi_name = std::string(poi_names[i], poi_names[i] + strlen(poi_names[i]));
        names.push_back(poi_name);
        names_map[poi_name] = i;
        free(poi_names[i]);
    }

    H5Tclose(datatype);
    H5Sclose(space);
    H5Dclose(poi_dataset);
    
}

void H5File::set_memspace(size_t rows, size_t cols) {
    hyperslab_dims[0] = rows;
    hyperslab_dims[1] = cols;
    if (memspace_id > -1) {
        H5Sclose(memspace_id);
    }
    memspace_id = H5Screate_simple(rank, hyperslab_dims, NULL);
}
void H5File::get_POI_matrix(
    FRMatrix& G,
    const std::vector<std::string>& poi_individuals,
    const std::vector<std::string>& poi_names,
    const int chunk_size
) {

    G.row_names = std::unordered_map<std::string, int>(poi_individuals.size());
    G.col_names = std::unordered_map<std::string, int>(poi_names.size());
    std::vector<hsize_t> row_indices(poi_individuals.size());
    std::vector<hsize_t> col_indices(poi_names.size());

    for (size_t i = 0; i < poi_individuals.size(); i++) {
        G.row_names[poi_individuals[i]] = i;
        row_indices[i] = individuals_map[poi_individuals[i]];
    }
    for (size_t i = 0; i < poi_names.size(); i++) {
        G.col_names[poi_names[i]] = i;
        col_indices[i] = names_map[poi_names[i]];
    }

    // Open the HDF5 file
    hid_t file_id = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Open the dataset
    hid_t dataset_id = H5Dopen(file_id, "values", H5P_DEFAULT);

    // Get the dataspace
    hid_t dataspace_id = H5Dget_space(dataset_id);
    // check for hdf5 memory buffer size - can this be optimized?
    G.data = arma::mat(poi_individuals.size(), chunk_size, arma::fill::zeros);

    // Define the hyperslab for the entire range of columns needed
    hsize_t src_offset[2] = {0, col_indices[0]};
    hsize_t dst_offset[2] = {0, 0};
    // n x m
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, src_offset, NULL, hyperslab_dims, NULL);

    // Create a memory dataspace
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, dst_offset, NULL, hyperslab_dims, NULL);

    // Read the data
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, G.data.memptr());

    // Close resources
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
}