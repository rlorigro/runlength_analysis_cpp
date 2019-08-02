#include <thread>
#include <vector>
#include <string>
#include <iostream>
#include "Runlength.hpp"
#include <experimental/filesystem>
#include <mutex>
#include <thread>
#include <atomic>
#include <exception>
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Runlength.hpp"

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


void runlength_encode(runlength_sequence_element& runlength_sequence, SequenceElement& sequence){
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


void runlength_encode_sequence_to_file(path& fasta_path, uint64_t read_index, string read_name, mutex& file_write_mutex, FastaWriter& fasta_writer, atomic<uint64_t>& n_jobs_running){
    n_jobs_running++;

    // Initialize containers
    SequenceElement sequence;
    runlength_sequence_element runlength_sequence;

    // Fetch Fasta sequence
    FastaReader fasta_reader = FastaReader(fasta_path);
    fasta_reader.fetch_sequence(sequence, read_name, read_index);

    // Convert to Run-length Encoded sequence element
    runlength_encode(runlength_sequence, sequence);

    // Write RLE sequence to file (no lengths written)
    file_write_mutex.lock();
    fasta_writer.write(runlength_sequence);
    file_write_mutex.unlock();

    // Print status update to stdout
    cerr << "\33[2K\rParsed: " << sequence.name << flush;

    n_jobs_running--;
}


void get_key_vector_from_map(vector<string>& keys, map <string, uint64_t>& map_object){
    for (auto& element: map_object){
        keys.push_back(element.first);
    }
}


path fasta_to_RLE_fasta(path input_file_path, path output_dir, uint64_t max_threads) {
    // Generate parent directories if necessary
    create_directories(output_dir);

    string output_filename;
    output_filename = string(input_file_path.filename());
    output_filename = output_filename.substr(0, output_filename.find_last_of(".")) + "_RLE.fasta";

    path output_file_path = output_dir / output_filename;

    cout << "READING FILE: " << input_file_path.string() << "\n";
    cout << "WRITING FILE: " << output_file_path.string() << "\n";

    // This reader is used to fetch an index
    FastaReader fasta_reader(input_file_path);

    // This writer is mutexed across threads
    FastaWriter fasta_writer(output_file_path);

    // Get index
    auto read_indexes = fasta_reader.get_index();

    vector<string> read_names;
    get_key_vector_from_map(read_names, read_indexes);

    atomic<uint64_t> n_jobs_running = 0;
    uint64_t job_index = 0;

    vector<thread> threads;
    mutex file_write_mutex;

    string read_name;
    uint64_t read_index;

    // Launch threads
    while (job_index < read_names.size()){
        if (n_jobs_running < max_threads) {
            // Get data to send to threads (must not be sent by reference, or will lead to concurrency issues)
            read_name = read_names[job_index];
            read_index = read_indexes[read_name];
            job_index++;

            try {
                // Call thread safe function to RL encode and write to file
                threads.emplace_back(thread(runlength_encode_sequence_to_file,
                                            ref(input_file_path),
                                            read_index,
                                            read_name,
                                            ref(file_write_mutex),
                                            ref(fasta_writer),
                                            ref(n_jobs_running)));
            } catch (const exception &e) {
                cerr << e.what() << "\n";
                exit(1);
            }
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }

    cerr << "\n" << flush;

    return output_file_path;
}
