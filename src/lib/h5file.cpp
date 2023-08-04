// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "hdf5.h"
#include <fr_matrix.h>
#if !defined(__APPLE__) && !defined(__MACH__)
  #include <omp.h>
#endif
#include <h5file.h>

void H5File::get_POI_individuals() {
    const char dname[] = "individuals";
    hid_t ind_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

    if (ind_dataset < 0) {
        Rcpp::Rcerr << "Failed to open individuals dataset." << std::endl;
        return;
    }

    hid_t datatype = H5Dget_type(ind_dataset);
    hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

    if (H5Tget_class(native_type) != H5T_STRING) {
        Rcpp::Rcerr << "Dataset does not have the expected string data type." << std::endl;
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(ind_dataset);
        return;
    }
    
    size_t datatype_size = H5Tget_size(native_type);
    hid_t space = H5Dget_space(ind_dataset);
    hsize_t num_ind;
    H5Sget_simple_extent_dims(space, &num_ind, NULL);
    individuals.reserve(num_ind);

    char* buffer = new char[num_ind * datatype_size];

    if (H5Dread(ind_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
        Rcpp::Rcerr << "Failed to read individuals dataset" << std::endl;
        delete[] buffer;
        H5Sclose(space);
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(ind_dataset);
        return;
    }

    for (hsize_t i = 0; i < num_ind; i++) {
        std::string name = std::string(&buffer[i * datatype_size], datatype_size);
        // Remove potential trailing white spaces due to fixed-length string format
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

void H5File::get_POI_names() {
    const char dname[] = "predictors_of_interest";
    hid_t poi_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

    if (poi_dataset < 0) {
        Rcpp::Rcerr << "Failed to open predictors_of_interest dataset." << std::endl;
        return;
    }

    hid_t datatype = H5Dget_type(poi_dataset);
    hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

    if (H5Tget_class(native_type) != H5T_STRING) {
        Rcpp::Rcerr << "Dataset does not have the expected string data type." << std::endl;
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(poi_dataset);
        return;
    }
    size_t datatype_size = H5Tget_size(native_type);
    hid_t space = H5Dget_space(poi_dataset);
    hsize_t num_poi;
    H5Sget_simple_extent_dims(space, &num_poi, NULL);
    names.reserve(num_poi);

    char* buffer = new char[num_poi * datatype_size];
    if (H5Dread(poi_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
        Rcpp::Rcerr << "Failed to read POI dataset" << std::endl;
        delete[] buffer;
        H5Sclose(space);
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(poi_dataset);
        return;
    }

    for (hsize_t i = 0; i < num_poi; i++) {
        std::string name = std::string(&buffer[i * datatype_size], datatype_size);
        // Remove potential trailing white spaces due to fixed-length string format
        name.erase(name.find_last_not_of(' ') + 1);
        names.push_back(name);
        names_map[name] = i;
    }

    delete[] buffer;
    H5Sclose(space);
    H5Tclose(native_type);
    H5Tclose(datatype);
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
// void H5File::get_POI_matrix(
//     FRMatrix& G,
//     const std::vector<std::string>& poi_individuals,
//     const std::vector<std::string>& poi_names,
//     const int chunk_size
// ) {

//     G.row_names = std::unordered_map<std::string, int>(poi_individuals.size());
//     G.col_names = std::unordered_map<std::string, int>(poi_names.size());
//     std::vector<hsize_t> row_indices(poi_individuals.size());
//     std::vector<hsize_t> col_indices(poi_names.size());

//     for (size_t i = 0; i < poi_individuals.size(); i++) {
//         G.row_names[poi_individuals[i]] = i;
//         row_indices[i] = individuals_map[poi_individuals[i]];
//     }
//     for (size_t i = 0; i < poi_names.size(); i++) {
//         G.col_names[poi_names[i]] = i;
//         col_indices[i] = names_map[poi_names[i]];
//     }

//     // Open the HDF5 file
//     hid_t file_id = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

//     // Open the dataset
//     hid_t dataset_id = H5Dopen(file_id, "values", H5P_DEFAULT);

//     // Get the dataspace
//     hid_t dataspace_id = H5Dget_space(dataset_id);
//     // check for hdf5 memory buffer size - can this be optimized?
//     G.data = arma::mat(poi_individuals.size(), chunk_size, arma::fill::zeros);

//     // Define the hyperslab for the entire range of columns needed
//     hsize_t src_offset[2] = {0, col_indices[0]};
//     hsize_t dst_offset[2] = {0, 0};
//     // n x m
//     H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, src_offset, NULL, hyperslab_dims, NULL);

//     // Create a memory dataspace
//     H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, dst_offset, NULL, hyperslab_dims, NULL);

//     // Read the data
//     H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, G.data.memptr());

//     // Close resources
//     H5Sclose(dataspace_id);
//     H5Dclose(dataset_id);
//     H5Fclose(file_id);
// }

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

    // Get the dataspace
    hid_t dataspace_id = H5Dget_space(values_dataset_id);
    hid_t datatype = H5Dget_type(values_dataset_id);
    H5T_class_t type_class = H5Tget_class(datatype);
    switch (type_class)
    {
    case H5T_INTEGER:
        Rcpp::Rcout << "type is Integer" << std::endl;
        break;
    case H5T_FLOAT: 
        Rcpp::Rcout << "type is float" << std::endl;
        break;
    case H5T_COMPOUND: 
        Rcpp::Rcout << "type is compound" << std::endl;
        break;
    case H5T_ARRAY: 
        Rcpp::Rcout << "type is array" << std::endl;
        break;
    case H5T_BITFIELD: 
        Rcpp::Rcout << "type is bitfield" << std::endl;
        break;
    case H5T_ENUM: 
        Rcpp::Rcout << "type is enum" << std::endl;
        break;
    case H5T_STRING: 
        Rcpp::Rcout << "type is string" << std::endl;
        break;
    default:
        Rcpp::Rcout << "type is unknown" << std::endl;
        break;
    }
    

    // Get dimensions of the dataspace
    int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(dataspace_id, dims.data(), NULL);

    bool transpose = false;
    if(poi_individuals.size() != dims[0] && poi_individuals.size() == dims[1]) {
        transpose = true;
    }

    if(poi_individuals.size() != dims[0] && poi_individuals.size() != dims[1]) {
        Rcpp::Rcerr << "Dimensions of the dataset do not match the sizes of poi_individuals. Please check the hdf5 file dataset dimensions." << std::endl;
        return;
    }

    // check for hdf5 memory buffer size - can this be optimized?
    G.data = arma::mat(poi_individuals.size(), chunk_size, arma::fill::zeros);

    hsize_t hyperslab_dims[2] = {poi_individuals.size(), poi_names.size()};
    
    // Define the hyperslab for the entire range of columns needed
    hsize_t src_offset[2] = {0, col_indices[0]};
    if (transpose) {
        std::swap(hyperslab_dims[0], hyperslab_dims[1]);
        std::swap(src_offset[0], src_offset[1]);
    }
    hsize_t dst_offset[2] = {0, 0};
    // n x m
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, src_offset, NULL, hyperslab_dims, NULL);

    // Create a memory dataspace
    hid_t memspace_id = H5Screate_simple(2, hyperslab_dims, NULL);
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, dst_offset, NULL, hyperslab_dims, NULL);

    // Read the data
    if (type_class == H5T_INTEGER) {
        size_t type_size = H5Tget_size(datatype);

        if (type_size == 4) {
            // Reading 32-bit integers
            arma::Mat<int32_t> tmp(poi_individuals.size(), poi_names.size(), arma::fill::zeros);
            // Reading the data directly into the matrix
            H5Dread(values_dataset_id, H5T_NATIVE_INT32, memspace_id, dataspace_id, H5P_DEFAULT, tmp.memptr());
            // Convert to arma::mat
            G.data = arma::conv_to<arma::mat>::from(tmp);
            G.data.replace(-2147483648, arma::datum::nan);
        }
    }
    else if (type_class == H5T_FLOAT) {
        size_t type_size = H5Tget_size(datatype);
        Rcpp::Rcout << "type size: " << type_size << std::endl;
        // Reading 64-bit double
        arma::Mat<double> tmp(poi_individuals.size(), poi_names.size(), arma::fill::zeros);
        // Reading the data directly into the matrix
        H5Dread(values_dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, tmp.memptr());
        G.data = arma::conv_to<arma::mat>::from(tmp);
        G.data(arma::span(0, 3), arma::span(0, 3)).print();
    }
    G.data.replace(-2.1475e+09, arma::datum::nan);
    // Close resources
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
}