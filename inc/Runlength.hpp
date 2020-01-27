#ifndef RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
#define RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP

#include "MarginPolishReader.hpp"
#include "ShastaReader.hpp"
#include "AlignedSegment.hpp"
#include "RunnieReader.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "BedReader.hpp"
#include "BamReader.hpp"
#include "Runlength.hpp"
#include "Matrix.hpp"
#include "Align.hpp"
#include "RunlengthSequenceElement.hpp"
#include "FastaReader.hpp"
#include "Matrix.hpp"
#include <vector>
#include <vector>
#include <thread>
#include <string>
#include <iostream>
#include <mutex>
#include <exception>
#include <atomic>
#include <experimental/filesystem>

using std::vector;
using std::string;
using std::thread;
using std::cout;
using std::cerr;
using std::flush;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::ref;
using std::move;
using std::exception;
using std::atomic;
using std::atomic_fetch_add;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;
using std::vector;


template<class T> void runlength_encode(RunlengthSequenceElement& runlength_sequence, T& sequence);


path runlength_encode_fasta_file(path input_file_path,
                                 unordered_map <string,RunlengthSequenceElement>& runlength_sequences,
                                 path output_dir,
                                 bool store_in_memory,
                                 uint16_t max_threads);

void measure_runlength_distribution_from_marginpolish(path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads,
        path bed_path=path());


void measure_runlength_distribution_from_shasta(path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads,
        path bed_path=path());


void measure_runlength_distribution_from_fasta(path reads_fasta_path,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads);


void measure_runlength_distribution_from_runnie(path runnie_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads);


void get_vector_from_index_map(vector< pair <string,FastaIndex> >& items, unordered_map<string,FastaIndex>& map_object);



/// TEMPLATE BASED METHODS ///

template<class T> void runlength_encode(RunlengthSequenceElement& runlength_sequence, T& sequence){
    runlength_sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    char current_character = 0;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (auto& character: sequence.sequence){
        if (tolower(character) != tolower(current_character)){
            runlength_sequence.sequence += character;
            runlength_sequence.lengths.push_back(1);
        }
        else{
            runlength_sequence.lengths.back()++;
        }

        current_character = character;
    }
}


template<class T> void runlength_decode(RunlengthSequenceElement& runlength_sequence, T& sequence){
    sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (size_t i=0; i<runlength_sequence.sequence.size(); i++){
        sequence.sequence += string(runlength_sequence.lengths[i], runlength_sequence.sequence[i]);
    }
}

template<class T> void runlength_decode(RunlengthSequenceElement& runlength_sequence, T& sequence, size_t start, size_t stop){
    sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (size_t i=start; i<stop; i++){
        sequence.sequence += string(runlength_sequence.lengths[i], runlength_sequence.sequence[i]);
    }
}


template <typename T> void write_segment_consensus_sequence_to_fasta(path& parent_directory,
                                               unordered_map<string,path>& read_paths,
                                               vector<string> read_names,
                                               mutex& file_write_mutex,
                                               FastaWriter& fasta_writer,
                                               atomic<uint64_t>& job_index){

    while (job_index < read_names.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);

        // Initialize containers
        CoverageSegment segment;

        // Fetch Fasta sequence
        T reader = T(parent_directory);
        reader.set_index(read_paths);
        reader.fetch_consensus_sequence(segment, read_names[thread_job_index]);

        // Write RLE sequence to file (no lengths written)
        file_write_mutex.lock();
        fasta_writer.write(segment);
        file_write_mutex.unlock();

        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << segment.name << flush;
    }
}


template <typename T> path write_all_consensus_sequences_to_fasta(T& coverage_reader,
        vector<string>& read_names,
        unordered_map<string,path>& read_paths,
        path input_directory,
        path output_directory,
        uint16_t max_threads){

    // Generate output file path
    path output_fasta_filename = absolute(input_directory);

    if (output_fasta_filename.stem() == ".") {
        output_fasta_filename = output_fasta_filename.parent_path().stem().string() + ".fasta";
    }
    else{
        output_fasta_filename = output_fasta_filename.stem().string() + ".fasta";
    }

    path output_fasta_path = output_directory / output_fasta_filename;

    cerr << "WRITING TO: " << output_fasta_filename <<'\n';

    // Initialize Fasta Writer
    FastaWriter read_fasta_writer = FastaWriter(output_fasta_path);

    // Thread-related variables
    vector<thread> threads;
    mutex file_write_mutex;

    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(write_segment_consensus_sequence_to_fasta<T>,
                                        ref(input_directory),
                                        ref(read_paths),
                                        ref(read_names),
                                        ref(file_write_mutex),
                                        ref(read_fasta_writer),
                                        ref(job_index)));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }
    cerr << "\n" << flush;

    return output_fasta_path;
}


void label_coverage_data_from_shasta(
        path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_threads,
        path bed_path,
        uint16_t insert_cutoff=32768);


#endif //RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
