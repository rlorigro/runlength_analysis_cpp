#include <iterator>
#include <numeric>
#include <future>
#include <atomic>
#include <vector>
#include <tuple>
#include <iostream>
#include <cmath>
#include <map>
#include <chrono>
#include <utility>
#include <thread>

using std::vector;
using std::accumulate;
using std::future;
using std::async;
using std::atomic;
using std::tuple;
using std::make_tuple;
using std::forward_as_tuple;
using std::cout;
using std::pow;
using std::ref;
using std::move;
using std::thread;
using std::get;


void generate_vector(vector <uint64_t>& v, uint64_t length, uint64_t start_value){
    for (uint64_t i=0; i<length; i++){
        v.push_back(start_value+1);
    }
}


void generate_vector_threaded(
        vector <vector <uint64_t> >& v,
        uint16_t thread_index,
        uint64_t length,
        uint64_t start_value,
        atomic<uint16_t>& n_jobs_finished){

    for (uint64_t i=0; i<length; i++){
        v.push_back(start_value+1);
    }

    n_jobs_finished++;
}


uint64_t sum_vector(vector<uint64_t>& data){
    uint64_t sum = 0;

    for (auto& i: data){
        sum += i;
    }

    return sum;
}


uint64_t sum_vector_threaded(vector<uint64_t>& data, atomic<uint16_t>& n_jobs_finished){
    uint64_t sum = 0;

    for (auto& i: data){
        sum += i;
    }

    n_jobs_finished++;
    return sum;
}


void wait_for_jobs_to_finish(atomic<uint16_t>& n_jobs_finished, uint16_t n_jobs){
    uint16_t i = 0;

    // Update stdout job counter every 1000 cycles
    while (n_jobs_finished < n_jobs){
        if (i % 1000){
            cout << "\r" << n_jobs_finished;
        }
        i++;
    }

    // Terminate with a newline
    cout << "\n";
}


void threaded_sum_test(uint64_t job_size, uint16_t n_threads){
    vector <tuple <uint64_t, uint64_t> > producer_job_queue;        // List of bounds that specify chunks of input
    tuple <uint64_t, uint64_t> chunk;                               // The necessary criteria to start processing a chunk
    vector <vector <uint64_t> > consumer_job_queue(n_threads);      // Loaded input data
    vector <vector <uint64_t> > consumer_result_queue(n_threads);   // List of "future"s waiting for .get()
    atomic <uint16_t> n_jobs_finished;                              // Counter to tell how many jobs have finished
    vector <thread> threads;

    uint64_t job_index = 0;

    // Chunk input data
    for (uint64_t i=0; i<n_threads; i++){
        // Prepare chunk boundaries to be produced
        chunk = make_tuple(job_size, job_index);
        producer_job_queue.push_back(chunk);

        job_index++;
    }

    // Produce
    n_jobs_finished = 0;
    for (uint16_t i=0; i<n_jobs; i++){
        threads.push_back(thread(generate_vector_threaded,
                                 ref(consumer_job_queue),

                                 get<0>(producer_job_queue[i]),
                                 get<1>(producer_job_queue[i]),
                                 ref(n_jobs_finished)));
    }

//    // Wait for jobs to finish
//    wait_for_jobs_to_finish(n_jobs_finished, n_jobs);
//
//    // Consume
//    n_jobs_finished = 0;
//    for (auto& data: consumer_job_queue){
//        auto result = async(sum_vector_threaded, ref(data), ref(n_jobs_finished));
//        consumer_result_queue.push_back(move(result));
//    }
//
//    // Wait for jobs to finish
//    wait_for_jobs_to_finish(n_jobs_finished, n_jobs);
//
//    // Fetch and write results
//    for (auto& result: consumer_result_queue){
//        cout << result.get() << "\n";
//    }
}


void unthreaded_sum_test(uint64_t job_size, uint16_t n_threads){
    vector <tuple <vector <uint64_t>, uint64_t, uint64_t> > producer_job_queue;
    tuple <vector <uint64_t>, uint64_t, uint64_t> chunk;
    vector <vector <uint64_t> > consumer_job_queue;
    vector <uint64_t> result_queue;
    uint64_t result;

    uint64_t job_index = 0;

    // Chunk input data
    for (uint64_t i=0; i<n_threads; i++){
        vector <uint64_t> data;
        chunk = make_tuple(data, job_size, job_index);
        producer_job_queue.push_back(move(chunk));
        job_index++;
    }

    // Produce
    for (auto& [data, job_size, job_index]: producer_job_queue){
        generate_vector(data, job_size, job_index);
        consumer_job_queue.push_back(move(data));
    }

    // Consume
    for (auto& item: consumer_job_queue){
        result = sum_vector(ref(item));
        result_queue.push_back(move(result));
    }

    // Write results
    for (auto& item: result_queue){
        cout << item << "\n";
    }
}


int main(){
    uint64_t job_size = 10000000;
    uint16_t n_jobs = 4;

    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    std::chrono::duration<double> elapsed_seconds;
    std::time_t end_time;

    start = std::chrono::system_clock::now();

    unthreaded_sum_test(job_size, n_threads);

    end = std::chrono::system_clock::now();

    elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "\n";

    start = std::chrono::system_clock::now();

    threaded_sum_test(job_size, n_jobs, n_threads);

    end = std::chrono::system_clock::now();

    elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "\n";
}
