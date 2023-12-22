

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
    hid_t file_id = -1;
    hid_t values_dataset_id = -1;
    hid_t values_dataspace_id = -1;
    hid_t values_datatype = -1;
    H5T_class_t values_type_class;
    int rank = 2;
    hid_t memspace_id = -1;
    hsize_t hyperslab_dims[2];

    std::string file_path;
    
    H5File(const std::string& file_name) {
        file_path = file_name;
        // H5Pset_fclose_degree(H5F_CLOSE_STRONG);
        // H5set_free_list_limits(0, 0, 0, 0, 0, 0);
        // file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        // values_dataset_id = H5Dopen(file_id, "values", H5P_DEFAULT);
        // if (file_id < 0) {
        //     Rcpp::stop("Failed to open HDF5 file.");
        // }
        // get_POI_individuals();
        // get_POI_names();
        // values_type_class = get_POI_data_type();
    }
    H5File(){}
    ~H5File() {
        H5garbage_collect();
        H5close();
    }
    void open_file(bool read_only = false);
    void get_values_dataset_id();
    void close_all();
    void get_POI_individuals();
    void get_POI_names();
    void get_POI_data_type();
    void set_memspace(size_t rows, size_t cols);
    void get_POI_matrix(
        FRMatrix& G,
        const std::vector<std::string>& poi_individuals,
        const std::vector<std::string>& poi_names,
        const int chunk_size
    );
};


// cmake -G "Unix Makefiles" -DHDF5_ENABLE_THREADSAFE=ON -DCMAKE_INSTALL_PREFIX="mnt/e/Development/hdf5-parallel/hdf5-threadsafe" -DHDF5_BUILD_CPP_LIB=ON -DALLOW_UNSUPPORTED=ON ..

// cmake -G "Visual Studio 17 2022" -DHDF5_ENABLE_THREADSAFE=ON -DCMAKE_INSTALL_PREFIX="mnt/e/Development/hdf5-parallel/hdf5-threadsafe" -DHDF5_BUILD_CPP_LIB=ON -DALLOW_UNSUPPORTED=ON ..
