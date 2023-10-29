#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>

struct Job {
    int id;
};

class Worker {
    private:
        int _id;
        std::queue<Job>& _job_queue;
        std::mutex& _mtx;
        std::condition_variable& _cv;
        std::ofstream& _output_stream;
        void process_job(const Job& job);
    public:
        Worker(int id, std::queue<Job>& job_queue, std::mutex& mtx, std::ofstream& output_stream, std::condition_variable& cv)
        : _id(id), _job_queue(job_queue), _mtx(mtx), _output_stream(output_stream), _cv(cv) {};

        void operator()(){};
};