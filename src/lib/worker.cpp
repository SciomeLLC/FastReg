#include <worker.h>

void Worker::operator()() {
    while(true) {
        Job job;
        {
            std::unique_lock<std::mutex> lock(_mtx);
            if(_job_queue.empty()) {
                _cv.notify_all();
                break;
            }
            job = _job_queue.front();
            _job_queue.pop();
        }
        process_job(job);
    }
}

void Worker::process_job(const Job& job) {
    // TODO: regression logic for job
    std::unique_lock<std::mutex> lock(_mtx);

}