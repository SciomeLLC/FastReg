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

void POI::get_individuals() {
    const char dname[] = "individuals";
    hid_t ind_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

    if (ind_dataset < 0) {
        Rcpp::stop("Failed to open individuals dataset.");
        return;
    }

    hid_t datatype = H5Dget_type(ind_dataset);
    hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

    if (H5Tget_class(native_type) != H5T_STRING) {
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(ind_dataset);
        Rcpp::stop("Dataset does not have the expected string data type.");
        return;
    }
    if (H5Tget_class(native_type) != H5T_STRING) {
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(ind_dataset);
        Rcpp::stop("Dataset does not have the expected variable-length string data type.");
        return;
    }

    if (H5Tis_variable_str(native_type)) {
        hid_t space = H5Dget_space(ind_dataset);
        hsize_t num_ind;
        H5Sget_simple_extent_dims(space, &num_ind, NULL);

        char** rdata = new char*[num_ind];

        if (H5Dread(ind_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata) < 0) {
            H5Tclose(native_type);
            H5Tclose(datatype);
            H5Dclose(ind_dataset);
            delete[] rdata;
            Rcpp::stop("Failed to read individuals dataset");
            return;
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

        char* buffer = new char[num_ind * datatype_size];

        if (H5Dread(ind_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
            delete[] buffer;
            H5Sclose(space);
            H5Tclose(native_type);
            H5Tclose(datatype);
            H5Dclose(ind_dataset);
            Rcpp::stop("Failed to read individuals dataset");
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
}

void POI::get_names() {
    const char dname[] = "predictors_of_interest";
    hid_t poi_dataset = H5Dopen(file_id, dname, H5P_DEFAULT);

    if (poi_dataset < 0) {
        Rcpp::stop("Failed to open predictors_of_interest dataset.");
        return;
    }

    hid_t datatype = H5Dget_type(poi_dataset);
    hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

    if (H5Tget_class(native_type) != H5T_STRING) {
        H5Tclose(native_type);
        H5Tclose(datatype);
        H5Dclose(poi_dataset);
        Rcpp::stop("Dataset does not have the expected string data type.");
        return;
    }

    if (H5Tis_variable_str(native_type)) {
        hid_t space = H5Dget_space(poi_dataset);
        hsize_t num_poi;
        H5Sget_simple_extent_dims(space, &num_poi, NULL);
        char** rdata = new char*[num_poi];

        if (H5Dread(poi_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata) < 0) {
            delete[] rdata;
            H5Tclose(native_type);
            H5Tclose(datatype);
            H5Dclose(poi_dataset);
            H5Sclose(space);
            Rcpp::stop("Failed to read individuals dataset");
            return;
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

        char* buffer = new char[num_poi * datatype_size];
        if (H5Dread(poi_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0) {
            delete[] buffer;
            H5Sclose(space);
            H5Tclose(native_type);
            H5Tclose(datatype);
            H5Dclose(poi_dataset);
            Rcpp::stop("Failed to read POI dataset");
            return;
        }

        for (hsize_t i = 0; i < num_poi; i++) {
            std::string name = std::string(&buffer[i * datatype_size], datatype_size);
            // Remove potential trailing white spaces due to fixed-length string format
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
    
    return;
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

void POI::load_data_chunk(
    FRMatrix& G,
    const std::vector<std::string>& poi_individuals,
    const std::vector<std::string>& poi_names,
    const int chunk_size
) {
    G.row_names.reserve(poi_individuals.size());
    G.col_names.reserve(poi_names.size());
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
    
    // Get dimensions of the dataspace
    values_dataspace_id = H5Dget_space(values_dataset_id);
    values_datatype = H5Dget_type(values_dataset_id);
    values_type_class = H5Tget_class(values_datatype);
    int ndims = H5Sget_simple_extent_ndims(values_dataspace_id);
    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(values_dataspace_id, dims.data(), NULL);
    bool transpose = false;
    if(poi_individuals.size() != dims[0] && poi_individuals.size() == dims[1]) {
        transpose = true;
    }

    if(poi_individuals.size() != dims[0] && poi_individuals.size() != dims[1]) {
        close_all();
        Rcpp::stop("Dimensions of the dataset do not match the sizes of poi_individuals. Please check the hdf5 file dataset dimensions.");
        return;
    }   

    hsize_t hyperslab_dims[2] = {poi_individuals.size(), poi_names.size()};
    
    // Define the hyperslab for the entire range of columns needed
    hsize_t src_offset[2] = {0, col_indices[0]};
    G.data.set_size(poi_individuals.size(), poi_names.size());
    if (transpose) {
        std::swap(hyperslab_dims[0], hyperslab_dims[1]);
        std::swap(src_offset[0], src_offset[1]);
        G.data.set_size(poi_names.size(), poi_individuals.size());
    }
    hsize_t dst_offset[2] = {0, 0};
    // n x m
    H5Sselect_hyperslab(values_dataspace_id, H5S_SELECT_SET, src_offset, NULL, hyperslab_dims, NULL);

    // Create a memory dataspace
    hid_t memspace_id = H5Screate_simple(2, hyperslab_dims, NULL);
    H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, dst_offset, NULL, hyperslab_dims, NULL);
    // Read the data
    Rcpp::Rcout << "Reading dataset" << std::endl;
    if (values_type_class == H5T_INTEGER) {
        arma::Mat<int32_t> tmp(hyperslab_dims[0], hyperslab_dims[1], arma::fill::zeros);
        // Reading the data directly into the matrix
        H5Dread(values_dataset_id, H5T_NATIVE_INT32, memspace_id, values_dataspace_id, H5P_DEFAULT, tmp.memptr());
        
        // Convert to arma::fmat
        G.data = arma::conv_to<arma::fmat>::from(tmp);
        G.data.replace(-2147483648, arma::datum::nan);
        if (transpose) {
            arma::inplace_trans(G.data); 
        }
    }
    else if (values_type_class == H5T_FLOAT) {
        // Reading 64-bit double
        // Reading the data directly into the matrix
        H5Dread(values_dataset_id, H5T_NATIVE_DOUBLE, memspace_id, values_dataspace_id, H5P_DEFAULT, G.data.memptr());
    }
    else {
        Rcpp::stop("HDF5 dataset type class is not float or int");
    }
    
    Rcpp::Rcout << "Read dataset" << std::endl;
    H5Sclose(memspace_id);
    memspace_id = -1;
}
void POI::close_all() {
    if(memspace_id > 0) { 
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

    if(file_id > 0){
        // H5garbage_collect();
        H5Fclose(file_id);
        file_id = -1;
    }
    // H5garbage_collect();
    // H5close();
}

void POI::open(bool read_only) {
    if  (read_only == false) {
        file_id = H5Fopen(file_path.c_str(), H5F_ACC_SWMR_WRITE, H5P_DEFAULT);
    }
    else {
        file_id = H5Fopen(file_path.c_str(), H5F_ACC_SWMR_READ | H5F_ACC_RDONLY, H5P_DEFAULT);       
    } 
    if (file_id < 0) {
        Rcpp::stop("Failed to open HDF5 file.");
    }
    // Rcpp::Rcout << "poi file opened" << std::endl;
}

void POI::get_values_dataset_id() {
    values_dataset_id = H5Dopen(file_id, "values", H5P_DEFAULT);
    // Rcpp::Rcout << "poi file values dataset read" << std::endl;

}
