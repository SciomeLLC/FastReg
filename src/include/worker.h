#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <ctime>
#include <chrono>
#include <condition_variable>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <h5file.h>
#include <fr_matrix.h>
#include <regression.h>
#include <config.h>
#ifndef UTILS_H
#include <utils.h>
#include <job.h>
#endif
#if !defined(__APPLE__) && !defined(__MACH__)
  #include <omp.h>
#endif


class Worker {
    private:
        int _id;
        std::queue<Job>& _job_queue;
        std::mutex& _mtx;
        std::condition_variable& _cv;
        // std::atomic<double> _nonconvergence_status_accumulator;
        // std::atomic<double> _total_filtered_pois_accumulator;
        std::atomic<int>& _poi_reading_time;
        std::atomic<int>& _file_writing_time;
        std::atomic<int>& _memory_allocation_time;
        std::atomic<int>& _regression_time;
        void process_job(Job& job);
    public:
        Worker(
          int id, 
          std::queue<Job>& job_queue, 
          std::mutex& mtx, 
          std::condition_variable& cv, 
          std::atomic<int>& memory_allocation_time,
          std::atomic<int>& file_writing_time,
          std::atomic<int>& poi_reading_time, 
          std::atomic<int>& regression_time
        ) : _id(id),
          _job_queue(job_queue),
          _mtx(mtx),
          _cv(cv),
          _memory_allocation_time(memory_allocation_time),
          _file_writing_time(file_writing_time),
          _poi_reading_time(poi_reading_time),
          _regression_time(regression_time) {}

        void operator()();
};