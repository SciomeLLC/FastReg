

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

#pragma once
class H5File {
public:
    std::vector<std::string> names;
    std::vector<std::string> individuals;
    std::unordered_map<std::string, int> individuals_map;
    std::unordered_map<std::string, int> names_map;
    hid_t file_id;
    hid_t values_dataset_id;
    hid_t values_dataspace_id;
    hid_t values_datatype;
    H5T_class_t values_type_class;
    int rank = 2;
    hid_t memspace_id = -1;
    hsize_t hyperslab_dims[2];

    std::string file_path;
    
    H5File(const std::string& file_name) {
        file_path = file_name;
        file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        values_dataset_id = H5Dopen(file_id, "values", H5P_DEFAULT);
        if (file_id < 0) {
            Rcpp::stop("Failed to open HDF5 file.");
        }
        get_POI_individuals();
        get_POI_names();
        values_type_class = get_POI_data_type();
    }

    ~H5File() {
        if(memspace_id > -1) { 
            H5Sclose(memspace_id);  
        }
        if(file_id >= 0){
            H5Fclose(file_id);
        }

        if (values_dataset_id >= 0) {
            H5Dclose(values_dataset_id);
        }

        if (values_dataspace_id >= 0) {
            H5Sclose(values_dataspace_id);
        }
    }
    void get_POI_individuals();
    void get_POI_names();
    H5T_class_t get_POI_data_type();
    void set_memspace(size_t rows, size_t cols);
    void get_POI_matrix(
        FRMatrix& G,
        const std::vector<std::string>& poi_individuals,
        const std::vector<std::string>& poi_names,
        const int chunk_size
    );
};
