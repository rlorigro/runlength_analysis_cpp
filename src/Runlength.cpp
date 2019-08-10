#include "MarginPolishReader.hpp"
#include "AlignedSegment.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "BamReader.hpp"
#include "Runlength.hpp"
#include "Align.hpp"
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


void runlength_encode_fasta_sequence_to_file_and_memory(path& fasta_path,
                                                        vector <pair <string, FastaIndex> >& read_index_vector,
                                                        unordered_map<string, RunlengthSequenceElement>& runlength_sequences,
                                                        mutex& map_mutex,
                                                        mutex& file_write_mutex,
                                                        FastaWriter& fasta_writer,
                                                        atomic<uint64_t>& job_index){

    while (job_index < read_index_vector.size()) {
        // Fetch add
        uint64_t thread_job_index = job_index.fetch_add(1);

        // Initialize containers
        SequenceElement sequence;
        RunlengthSequenceElement runlength_sequence;

        // Fetch Fasta sequence
        FastaReader fasta_reader = FastaReader(fasta_path);

        // Get read name and index
        string read_name = read_index_vector[thread_job_index].first;
        uint64_t read_index = read_index_vector[thread_job_index].second.byte_index;

        // First of each element is read_name, second is its index
        fasta_reader.fetch_sequence(sequence, read_name, read_index);

        // Convert to Run-length Encoded sequence element
        runlength_encode(runlength_sequence, sequence);

        // Write RLE sequence to file (no lengths written)
        file_write_mutex.lock();
        fasta_writer.write(runlength_sequence);
        file_write_mutex.unlock();

        // Append the sequence to a map of names:sequence
        map_mutex.lock();
        runlength_sequences[sequence.name] = move(runlength_sequence);
        map_mutex.unlock();

        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << sequence.name << flush;
    }
}


void write_marginpolish_consensus_sequence_to_fasta(path& marginpolish_directory,
                                               vector<string> read_names,
                                               mutex& file_write_mutex,
                                               FastaWriter& fasta_writer,
                                               atomic<uint64_t>& job_index){

    while (job_index < read_names.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);

        // Initialize containers
        MarginPolishSegment marginpolish_segment;

        // Fetch Fasta sequence
        MarginPolishReader marginpolish_reader = MarginPolishReader(marginpolish_directory);
        marginpolish_reader.fetch_consensus_sequence(marginpolish_segment, read_names[thread_job_index]);

        // Write RLE sequence to file (no lengths written)
        file_write_mutex.lock();
        fasta_writer.write(marginpolish_segment);
        file_write_mutex.unlock();

        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << marginpolish_segment.name << flush;
    }
}


//void get_key_vector_from_map(vector<string>& keys, unordered_map<string, uint64_t>& map_object){
//    for (auto& element: map_object){
//        keys.push_back(element.first);
//    }
//}


void get_vector_from_index_map(vector< pair <string,FastaIndex> >& items, unordered_map<string,FastaIndex>& map_object){
    for (auto& item: map_object){
        items.emplace_back(item.first, item.second);
    }
}


path runlength_encode_fasta_file(path input_file_path,
                                 unordered_map <string, RunlengthSequenceElement>& runlength_sequences,
                                 path output_dir,
                                 uint16_t max_threads) {

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

    // Convert the map object into an indexable object
    vector <pair <string, FastaIndex> > read_index_vector;
    get_vector_from_index_map(read_index_vector, read_indexes);

    mutex map_mutex;
    mutex file_write_mutex;

    atomic<uint64_t> job_index = 0;
    vector<thread> threads;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        // Get data to send to threads (must not be sent by reference, or will lead to concurrency issues)
        try {
            // Call thread safe function to RL encode and write to file
            threads.emplace_back(thread(runlength_encode_fasta_sequence_to_file_and_memory,
                                        ref(input_file_path),
                                        ref(read_index_vector),
                                        ref(runlength_sequences),
                                        ref(map_mutex),
                                        ref(file_write_mutex),
                                        ref(fasta_writer),
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

    return output_file_path;
}


path write_all_marginpolish_consensus_sequences_to_fasta(MarginPolishReader& mp_reader,
        vector<string>& read_names,
        path marginpolish_directory,
        path output_directory,
        uint16_t max_threads){

    // Generate output file path
    path output_fasta_filename = marginpolish_directory;
    output_fasta_filename = output_fasta_filename.stem().string() + ".fasta";
    path output_fasta_path = output_directory / output_fasta_filename;

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
            threads.emplace_back(thread(write_marginpolish_consensus_sequence_to_fasta,
                                        ref(marginpolish_directory),
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


void update_runlength_matrix(AlignedSegment& aligned_segment,
        RunlengthSequenceElement& ref_runlength_sequence,
        RunlengthSequenceElement& read_runlength_sequence,
        vector <vector <uint16_t> >& runlength_matrix) {
    ///
    /// Iterate:
    ///   - AlignedSegment (cigars)
    ///   - A vector of lengths associated with the read
    ///   - A vector of lengths associated with the reference
    ///   - A string associated with the read
    ///
    /// to update (not construct) a matrix of observed vs true frequencies
    ///

}


void parse_region(path bam_path,
    unordered_map<string,RunlengthSequenceElement>& ref_runlength_sequences,
    vector <Region>& regions,
    vector <vector <uint16_t> >& runlength_matrix,
    atomic <uint64_t>& job_index){
    ///
    ///
    ///

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

    uint64_t thread_job_index;
    string ref_name;
    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;

    // Only allow matches
    unordered_set<uint8_t> valid_cigar_codes = {Cigar::cigar_code_key.at("=")};

    while (job_index < regions.size()) {
        thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        bam_reader.initialize_region(region.name, region.start, region.stop);

        int i = 0;

        while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
            cigars += cigar.get_cigar_code_as_string();

            // Subset alignment to portions of the read that are within the window/region
            if ((region.start <= uint64_t(coordinate.ref_index - 1)) and (uint64_t(coordinate.ref_index - 1) < region.stop)) {
                if (cigar.is_ref_move()) {
                    ref_alignment += ref_runlength_sequences.at(region.name).sequence[coordinate.ref_index];
                } else {
                    ref_alignment += "-";
                }
                if (cigar.is_read_move()) {
                    read_alignment += aligned_segment.get_read_base(coordinate.read_index);
                    read_alignment_inferred += "_";
                } else {
                    read_alignment += "*";
                    read_alignment_inferred += "*";
                }
            }
            i++;
        }

        if (i > 0) {
            cout << region.to_string() << "\n" << flush;
            cout << aligned_segment.to_string();
            cout << cigars << "\n";
            cout << ref_alignment << "\n";
            cout << read_alignment << "\n";
            cout << read_alignment_inferred << "\n";
            cout << "\n\n";
        }
    }
}


void chunk_sequences_into_regions(vector<Region>& regions, unordered_map<string,RunlengthSequenceElement>& sequences, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& [name, item]: sequences){
        chunk_sequence(regions, name, chunk_size, item.sequence.size());
    }
}


void measure_runlength_distribution(path bam_path,
        unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
        vector <Region>& regions,
        uint16_t max_runlength,
        uint16_t max_threads){
    ///
    ///
    ///

    // Initialize template 2d matrix with 0
    vector <vector <uint16_t> > matrix(max_runlength, vector<uint16_t>(max_runlength, 0));

    // Make as many matrices as there are threads
    vector <vector <vector <uint16_t> > > runlength_matrices_per_thread(max_threads, matrix);

    vector<thread> threads;
    atomic<uint64_t> job_index;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(parse_region,
                    ref(bam_path),
                    ref(ref_runlength_sequences),
                    ref(regions),
                    ref(matrix),
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
}


void measure_runlength_distribution_from_marginpolish(path marginpolish_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads){

    cout << "USING " + to_string(max_threads) + " THREADS\n";

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;

    MarginPolishReader marginpolish_reader = MarginPolishReader(marginpolish_directory);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Runlength encode the reference (store in memory)
    unordered_map<string,RunlengthSequenceElement> ref_runlength_sequences;
    path reference_fasta_path_rle;
    reference_fasta_path_rle = runlength_encode_fasta_file(reference_fasta_path, ref_runlength_sequences, output_directory, max_threads);

    marginpolish_reader.index();
    unordered_map<string,path> read_paths = marginpolish_reader.get_index();

    // Extract read names from index
    vector<string> read_names;
    for (auto& element: read_paths){
        read_names.push_back(element.first);
    }

    // Extract the sequences from MarginPolish TSVs (ignore coverage data for now)
    path reads_fasta_path_rle;
    reads_fasta_path_rle = write_all_marginpolish_consensus_sequences_to_fasta(marginpolish_reader, read_names, marginpolish_directory, output_directory, max_threads);

    // Index all the files in the MP directory
    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    uint16_t k = 19;
    string minimap_preset = "asm20";    //TODO: make command line argument
    bool explicit_mismatch = true;

    // Align Marginpolish reads to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle, reads_fasta_path_rle, output_directory, sort, index, delete_intermediates, k, minimap_preset, explicit_mismatch, max_threads);

    //TODO: cut off 50bp overlaps on MarginPolish TSVs!!

    // Chunk alignment regions
    vector<Region> regions;
    chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);

    // Launch threads for parsing alignments and generating matrices
    measure_runlength_distribution(bam_path, ref_runlength_sequences, regions, max_runlength, 1);

    // Save matrices
}
