

// [[Rcpp::depends(RcppArmadillo)]]
#include "hdf5.h"
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
class POI {
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
    bool transpose = false;
    hid_t memspace_id = -1;
    hsize_t hyperslab_dims[2];

    std::string file_path;
    
    POI(const std::string& file_name) {
        file_path = file_name;
    }
    POI(){}
    ~POI() {
        H5garbage_collect();
        H5close();
    }
    void open(bool read_only = false);
    void get_values_dataset_id();
    void close_all();
    void get_individuals();
    void get_names();
    void get_data_type();
    void set_memspace(size_t rows, size_t cols);
    void load_data_chunk(
        FRMatrix& G,
        const std::vector<std::string>& poi_individuals,
        const std::vector<std::string>& poi_names
    );
    void load_int_data_chunk(FRMatrix& G,  hsize_t* hyperslab_dims,  hsize_t* src_offset);
    void load_float_data_chunk(FRMatrix& G,  hsize_t* hyperslab_dims,  hsize_t* src_offset);
};
