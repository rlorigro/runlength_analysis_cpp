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
using std::max;
using std::min;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;


void chunk_sequences_into_regions(vector<Region>& regions, unordered_map<string,RunlengthSequenceElement>& sequences, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& [name, item]: sequences){
        chunk_sequence(regions, name, chunk_size, item.sequence.size());
    }
}


void write_matrix_to_file(path output_directory, runlength_matrix& matrix){
    path directional_matrix_path = absolute(output_directory) / "frequency_matrix_directional.csv";
    ofstream directional_matrix_file = ofstream(directional_matrix_path);

    path nondirectional_matrix_path = absolute(output_directory) / "frequency_matrix_nondirectional.csv";
    ofstream nondirectional_matrix_file = ofstream(nondirectional_matrix_path);

    if (not directional_matrix_file.is_open()){
        throw runtime_error("ERROR: file could not be written: " + directional_matrix_path.string());
    }

    if (not nondirectional_matrix_file.is_open()){
        throw runtime_error("ERROR: file could not be written: " + nondirectional_matrix_path.string());
    }

    cerr << "WRITING: matrix file " + directional_matrix_path.string() << '\n';
    cerr << "WRITING: matrix file " + nondirectional_matrix_path.string() << '\n';

    runlength_matrix nondirectional_matrix = sum_reverse_complements(matrix);

    directional_matrix_file << (matrix_to_string(matrix));
    nondirectional_matrix_file << (matrix_to_string(nondirectional_matrix));
}


void runlength_encode_fasta_sequence_to_file(path& fasta_path,
                                             vector <pair <string, FastaIndex> >& read_index_vector,
                                             unordered_map<string, RunlengthSequenceElement>& runlength_sequences,
                                             mutex& map_mutex,
                                             mutex& file_write_mutex,
                                             FastaWriter& fasta_writer,
                                             bool store_in_memory,
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
        fasta_reader.get_sequence(sequence, read_name, read_index);

        // Convert to Run-length Encoded sequence element
        runlength_encode(runlength_sequence, sequence);

        // Write RLE sequence to file (no lengths written)
        file_write_mutex.lock();
        fasta_writer.write(runlength_sequence);
        file_write_mutex.unlock();

        if (store_in_memory) {
            // Append the sequence to a map of names:sequence
            map_mutex.lock();
            runlength_sequences[sequence.name] = move(runlength_sequence);
            map_mutex.unlock();
        }
        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << sequence.name << flush;
    }
}


void write_runnie_sequence_to_fasta(path& runnie_directory,
        vector<string> read_names,
        unordered_map<string,RunnieIndex> read_indexes,
        mutex& file_write_mutex,
        FastaWriter& fasta_writer,
        atomic<uint64_t>& job_index){

    while (job_index < read_names.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);

        // Initialize containers
        RunnieSequenceElement runnie_sequence;

        // Fetch Fasta sequence
        RunnieReader runnie_reader = RunnieReader(runnie_directory);
        runnie_reader.set_index(read_indexes);
        runnie_reader.fetch_sequence_bases(runnie_sequence, read_names[thread_job_index]);

        // Write RLE sequence to file (no lengths written)
        file_write_mutex.lock();
        fasta_writer.write(runnie_sequence);
        file_write_mutex.unlock();

        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << runnie_sequence.name << flush;
    }
}


path runlength_encode_fasta_file(path input_file_path,
                                 unordered_map <string, RunlengthSequenceElement>& runlength_sequences,
                                 path output_dir,
                                 bool store_in_memory,
                                 uint16_t max_threads) {

    // Generate parent directories if necessary
    create_directories(output_dir);

    string output_filename;
    output_filename = string(input_file_path.filename());
    output_filename = output_filename.substr(0, output_filename.find_last_of(".")) + "_RLE.fasta";

    path output_file_path = output_dir / output_filename;

    cerr << "READING FILE: " << input_file_path.string() << "\n";
    cerr << "WRITING FILE: " << output_file_path.string() << "\n";

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
            threads.emplace_back(thread(runlength_encode_fasta_sequence_to_file,
                                        ref(input_file_path),
                                        ref(read_index_vector),
                                        ref(runlength_sequences),
                                        ref(map_mutex),
                                        ref(file_write_mutex),
                                        ref(fasta_writer),
                                        store_in_memory,
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


path write_all_runnie_sequences_to_fasta(RunnieReader& runnie_reader,
        vector<string>& read_names,
        unordered_map<string,RunnieIndex>& read_indexes,
        path runnie_directory,
        path output_directory,
        uint16_t max_threads){

    // Generate output file path
    path output_fasta_filename = runnie_directory;
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
            threads.emplace_back(thread(write_runnie_sequence_to_fasta,
                                        ref(runnie_directory),
                                        ref(read_names),
                                        ref(read_indexes),
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


ofstream& operator<<(ofstream& output_file, vector<CoverageElement>& coverage_data){
    for (auto &coverage_element: coverage_data) {
        output_file << coverage_element.base +
        to_string(coverage_element.length) +
        CoverageElement::reversal_string_plus_minus[coverage_element.reversal] + ' ' +
        to_string(uint32_t(coverage_element.weight)) + ',';
    }

    return output_file;
}


char bool_to_vertex_label(bool is_vertex){
    if (is_vertex){
        return 'v';
    }
    else{
        return 'e';
    }
}


char bool_to_reversal_char(bool reversal){
    if (reversal){
        return '-';
    }
    else{
        return '+';
    }
}


void write_labeled_coverage_data(
        ofstream& output_file,
        string true_base,
        uint16_t true_length,
        char consensus_base,
        uint16_t consensus_length,
        vector<CoverageElement>& coverage_data,
        bool is_vertex,
        bool reversal){

    if (reversal and consensus_base != '_'){
        consensus_base = complement_base(consensus_base);
    }

    output_file << bool_to_reversal_char(reversal) << ',' << true_base << ',' << true_length << ',' << consensus_base << ',' << consensus_length << ',';
    output_file << coverage_data << " ";
    output_file << bool_to_vertex_label(is_vertex);
    output_file << '\n';
}


template<typename T> void label_aligned_coverage(path bam_path,
                                                 path parent_directory,
                                                 unordered_map <string,path>& read_paths,
                                                 unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
                                                 vector <Region>& regions,
                                                 path output_directory,
                                                 atomic <uint64_t>& job_index){
    ///
    /// Create a duplicate series of CSVs which contain the true base and length as aligned to reference
    ///

    // Initialize Reader and relevant containers
    T reader = T(parent_directory);
    reader.store_length_consensus = true;
    reader.store_vertex_labels = true;
    CoverageSegment segment;
    reader.set_index(read_paths);

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

//    uint64_t thread_job_index;
    string ref_name;

    bool filter_secondary = true;
    uint16_t map_quality_cutoff = 5;

    // Volatiles
    bool in_left_bound;
    bool in_right_bound;
    string true_base;
    char consensus_base = '_';
    uint16_t consensus_length = -1;
    uint16_t true_length = -1;
    bool is_vertex = false;
    vector<CoverageElement> coverage_data;

    uint8_t match_code = Cigar::cigar_code_key.at("=");
    uint8_t mismatch_code = Cigar::cigar_code_key.at("X");
    uint8_t insert_code = Cigar::cigar_code_key.at("I");
    uint8_t delete_code = Cigar::cigar_code_key.at("D");

    // Only allow matches
    unordered_set<uint8_t> valid_cigar_codes = {match_code,
                                                mismatch_code,
                                                insert_code,
                                                delete_code};

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        // BAM coords are 1 based
        bam_reader.initialize_region(region.name, region.start+1, region.stop+1);

        int i = 0;

        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary)) {
            reader.fetch_read(segment, aligned_segment.read_name);

            int64_t read_start = max(int64_t(aligned_segment.ref_start_index), int64_t(region.start));
            int64_t read_stop = min(int64_t(aligned_segment.infer_reference_stop_position_from_alignment()), int64_t(region.stop));
            string read_region = region.name + "_" + to_string(read_start) + "-" + to_string(read_stop);

            path output_path = output_directory / (aligned_segment.read_name + "_" + read_region + ".csv");
            output_path = absolute(output_path);
            ofstream output_file(output_path);

            // Iterate cigars that match the criteria
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
                in_left_bound = (int64_t(region.start) <= coordinate.ref_index - 1);
                in_right_bound = (coordinate.ref_index - 1 < int64_t(region.stop));

                // Subset alignment to portions of the read that are within the window/region
                if (in_left_bound and in_right_bound) {
                    coverage_data = segment.coverage_data[coordinate.read_true_index];

                    /// MATCH OR MISMATCH
                    if (cigar.code == match_code or cigar.code == mismatch_code) {
                        true_base = ref_runlength_sequences.at(aligned_segment.ref_name).sequence[coordinate.ref_index];
                        true_length = ref_runlength_sequences.at(aligned_segment.ref_name).lengths[coordinate.ref_index];
                        consensus_base = segment.sequence[coordinate.read_true_index];
                        consensus_length = segment.lengths[coordinate.read_true_index];
                        is_vertex = segment.is_vertex[coordinate.read_true_index];
                    }
                    /// INSERT
                    else if (cigar.code == insert_code) {
                        true_base = '_';
                        true_length = 0;
                        consensus_base = segment.sequence[coordinate.read_true_index];
                        consensus_length = segment.lengths[coordinate.read_true_index];
                        is_vertex = segment.is_vertex[coordinate.read_true_index];
                    }
                    /// DELETE
                    else if (cigar.code == delete_code) {
                        true_base = ref_runlength_sequences.at(aligned_segment.ref_name).sequence[coordinate.ref_index];
                        true_length = ref_runlength_sequences.at(aligned_segment.ref_name).lengths[coordinate.ref_index];
                        consensus_base = '_';
                        consensus_length = 0;
                        coverage_data = {};
                        is_vertex = false;
                    }

                    write_labeled_coverage_data(
                            output_file,
                            true_base,
                            true_length,
                            consensus_base,
                            consensus_length,
                            coverage_data,
                            is_vertex,
                            aligned_segment.reversal);
                }
            }

            cerr << "\33[2K\rParsed: " << output_path << flush;

            i++;
        }
    }
}



void parse_aligned_runnie(path bam_path,
        path runnie_directory,
        unordered_map <string,RunnieIndex>& read_indexes,
        unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
        vector <Region>& regions,
        runlength_matrix& runlength_matrix,
        atomic <uint64_t>& job_index){
    ///
    ///
    ///

    // Initialize MarginPolishReader and relevant containers
    RunnieReader runnie_reader = RunnieReader(runnie_directory);
    RunnieSequenceElement runnie_sequence;
    runnie_reader.set_index(read_indexes);

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

    auto matrix_shape = runlength_matrix.shape();
    size_t max_true_length = matrix_shape[2];

    string ref_name;

    bool filter_secondary = true;
    uint16_t map_quality_cutoff = 5;

    // Volatiles
    bool in_left_bound;
    bool in_right_bound;
    string true_base;
    string observed_base;
    string consensus_base;
    float scale;
    float shape;
    uint16_t true_length = -1;
    uint8_t observed_base_index;
    vector<CoverageElement> coverage_data;

    // Only allow matches
    unordered_set<uint8_t> valid_cigar_codes = {Cigar::cigar_code_key.at("=")};

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        // BAM coords are 1 based
        bam_reader.initialize_region(region.name, region.start+1, region.stop+1);

        int i = 0;

        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary)) {
            runnie_reader.fetch_sequence(runnie_sequence, aligned_segment.read_name);

            // Iterate cigars that match the criteria (must be '=')
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
                in_left_bound = (int64_t(region.start) <= coordinate.ref_index - 1);
                in_right_bound = (coordinate.ref_index - 1 < int64_t(region.stop));

                // Subset alignment to portions of the read that are within the window/region
                if (in_left_bound and in_right_bound) {
                    true_base = ref_runlength_sequences.at(aligned_segment.ref_name).sequence[coordinate.ref_index];

                    // Skip anything other than ACTG
                    if (not is_valid_base(true_base)){
                        continue;
                    }

                    true_length = ref_runlength_sequences.at(aligned_segment.ref_name).lengths[coordinate.ref_index];
                    observed_base = runnie_sequence.sequence[coordinate.read_true_index];
                    observed_base_index = base_to_index(true_base);
                    scale = runnie_sequence.scales[coordinate.read_true_index];
                    shape = runnie_sequence.shapes[coordinate.read_true_index];

                    if (true_length >= max_true_length){
                        continue;
                    }

                    update_runlength_matrix_with_weibull_probabilities(runlength_matrix, aligned_segment.reversal, observed_base_index, true_length, scale, shape);
                }
            }

            i++;
        }

        cerr << "\33[2K\rParsed: " << region.to_string() << flush;
    }

}


template<typename T> void parse_aligned_coverage(path bam_path,
                                                 path parent_directory,
                                                 unordered_map <string,path>& read_paths,
                                                 unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
                                                 vector <Region>& regions,
                                                 runlength_matrix& runlength_matrix,
                                                 atomic <uint64_t>& job_index){
    ///
    ///
    ///

    // Initialize Reader and relevant containers
    T reader = T(parent_directory);
    CoverageSegment segment;
    reader.set_index(read_paths);

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

    auto matrix_shape = runlength_matrix.shape();
    size_t max_true_length = matrix_shape[2];
    size_t max_observed_length = matrix_shape[3];

//    uint64_t thread_job_index;
    string ref_name;

    bool filter_secondary = true;
    uint16_t map_quality_cutoff = 5;

    // Volatiles
    bool in_left_bound;
    bool in_right_bound;
    string true_base;
    char consensus_base;
    uint16_t true_length = -1;
    uint16_t observed_length = -1;
    uint8_t observed_base_index;
    vector<CoverageElement> coverage_data;

    // Only allow matches
    unordered_set<uint8_t> valid_cigar_codes = {Cigar::cigar_code_key.at("=")};

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        // BAM coords are 1 based
        bam_reader.initialize_region(region.name, region.start+1, region.stop+1);

        int i = 0;

        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary)) {
            reader.fetch_read(segment, aligned_segment.read_name);

            // Iterate cigars that match the criteria (must be '=')
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
                in_left_bound = (int64_t(region.start) <= coordinate.ref_index - 1);
                in_right_bound = (coordinate.ref_index - 1 < int64_t(region.stop));

                true_base = ref_runlength_sequences.at(aligned_segment.ref_name).sequence[coordinate.ref_index];
                consensus_base = segment.sequence[coordinate.read_true_index];

                // Subset alignment to portions of the read that are within the window/region
                if (in_left_bound and in_right_bound) {
                    true_length = ref_runlength_sequences.at(aligned_segment.ref_name).lengths[coordinate.ref_index];
                    coverage_data = segment.coverage_data[coordinate.read_true_index];

                    // Walk through all the coverage data for this position and update the matrix for each observation
                    // Only matches are allowed so checking the read base effectively checks the true base
                    for (auto& coverage_element: coverage_data){

                        // Skip anything other than ACTG
                        if (not coverage_element.is_conventional_base()){
                            continue;
                        }

                        // Do not record bases that were mismatches within the coverage data
                        if (coverage_element.base != consensus_base){
                            continue;
                        }

                        observed_length = coverage_element.length;
                        observed_base_index = coverage_element.get_base_index();

                        if (observed_length >= max_observed_length or true_length >= max_true_length){
                            continue;
                        }

                        runlength_matrix[coverage_element.reversal][observed_base_index][true_length][observed_length] += coverage_element.weight;
                    }
                }
            }

            i++;
        }

        cerr << "\33[2K\rParsed: " << region.to_string() << flush;
    }
}


void parse_aligned_fasta(path bam_path,
        path reads_fasta_path,
        unordered_map <string,FastaIndex>& read_indexes,
        unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
        vector <Region>& regions,
        runlength_matrix& runlength_matrix,
        atomic <uint64_t>& job_index){
    ///
    ///
    ///

    // Initialize FastaReader and relevant containers
    FastaReader fasta_reader = FastaReader(reads_fasta_path);
    SequenceElement sequence;
    RunlengthSequenceElement runlength_sequence;
    fasta_reader.set_index(read_indexes);

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

    auto matrix_shape = runlength_matrix.shape();
    size_t max_true_length = matrix_shape[2];
    size_t max_observed_length = matrix_shape[3];

    string ref_name;

    bool filter_secondary = true;
    uint16_t map_quality_cutoff = 5;

    // Volatiles
    bool in_left_bound;
    bool in_right_bound;
    string true_base;
    string observed_base;
    uint16_t true_length = -1;
    uint16_t observed_length = -1;
    uint8_t observed_base_index;
    vector<CoverageElement> coverage_data;

    // Only allow matches
    unordered_set<uint8_t> valid_cigar_codes = {Cigar::cigar_code_key.at("=")};

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        // BAM coords are 1 based
        bam_reader.initialize_region(region.name, region.start+1, region.stop+1);

        int i = 0;

        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary)) {
            fasta_reader.get_sequence(sequence, aligned_segment.read_name);
            runlength_encode(runlength_sequence, sequence);

            // Iterate cigars that match the criteria (must be '=')
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
                in_left_bound = (int64_t(region.start) <= coordinate.ref_index - 1);
                in_right_bound = (coordinate.ref_index - 1 < int64_t(region.stop));

                // Subset alignment to portions of the read that are within the window/region
                if (in_left_bound and in_right_bound) {
                    true_base = ref_runlength_sequences.at(aligned_segment.ref_name).sequence[coordinate.ref_index];
                    observed_base = runlength_sequence.sequence[coordinate.read_true_index];
                    true_length = ref_runlength_sequences.at(aligned_segment.ref_name).lengths[coordinate.ref_index];

                    // Skip anything other than ACTG
                    if (not is_valid_base(true_base)){
                        continue;
                    }

                    observed_length = runlength_sequence.lengths[coordinate.read_true_index];
                    observed_base_index = base_to_index(true_base);

                    if (observed_length >= max_observed_length or true_length >= max_true_length){
                        continue;
                    }

                    runlength_matrix[aligned_segment.reversal][observed_base_index][true_length][observed_length] += 1;
                }
            }

            i++;
        }

        cerr << "\33[2K\rParsed: " << region.to_string() << flush;
    }
}


template <typename T> void get_coverage_labels(path bam_path,
                                       path input_directory,
                                       path output_directory,
                                       unordered_map <string,path>& read_paths,
                                       unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
                                       vector <Region>& regions,
                                       uint16_t max_threads){
    ///
    ///
    ///

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(label_aligned_coverage<T>,
                                        ref(bam_path),
                                        ref(input_directory),
                                        ref(read_paths),
                                        ref(ref_runlength_sequences),
                                        ref(regions),
                                        ref(output_directory),
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


runlength_matrix get_runnie_runlength_matrix(path bam_path,
        path runnie_directory,
        unordered_map <string,RunnieIndex>& read_indexes,
        unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
        vector <Region>& regions,
        uint16_t max_runlength,
        uint16_t max_threads){
    ///
    ///
    ///

    runlength_matrix template_matrix(boost::extents[2][4][max_runlength+1][max_runlength+1]);   // 0 length included
    vector<runlength_matrix> matrices_per_thread(max_threads, template_matrix);

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(parse_aligned_runnie,
                                        ref(bam_path),
                                        ref(runnie_directory),
                                        ref(read_indexes),
                                        ref(ref_runlength_sequences),
                                        ref(regions),
                                        ref(matrices_per_thread[i]),
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

    cerr << "Summing matrices from " << max_threads << " threads...\n";

//    int i = 0;
//    for (auto& m: matrices_per_thread){
//        i++;
//        cout << "MATRIX " << i << ":\n" << matrix_to_string(m, 4) << "\n";
//    }

    runlength_matrix matrix_sum = sum_matrices(matrices_per_thread);

    return matrix_sum;
}


template <typename T> runlength_matrix get_runlength_matrix(path bam_path,
                                       path input_directory,
                                       unordered_map <string,path>& read_paths,
                                       unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
                                       vector <Region>& regions,
                                       uint16_t max_runlength,
                                       uint16_t max_threads){
    ///
    ///
    ///

    runlength_matrix template_matrix(boost::extents[2][4][max_runlength+1][max_runlength+1]);   // 0 length included
    vector<runlength_matrix> matrices_per_thread(max_threads, template_matrix);

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(parse_aligned_coverage<T>,
                                        ref(bam_path),
                                        ref(input_directory),
                                        ref(read_paths),
                                        ref(ref_runlength_sequences),
                                        ref(regions),
                                        ref(matrices_per_thread[i]),
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

    cerr << "Summing matrices from " << max_threads << " threads...\n";

//    int i = 0;
//    for (auto& matrix: matrices_per_thread){
//        i++;
//        cout << i << '\n' << matrix_to_string(matrix, 10) << '\n';
//    }

    runlength_matrix matrix_sum = sum_matrices(matrices_per_thread);

    return matrix_sum;
}


runlength_matrix get_fasta_runlength_matrix(path bam_path,
                                       path reads_fasta_path,
                                       unordered_map <string,FastaIndex>& read_indexes,
                                       unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
                                       vector <Region>& regions,
                                       uint16_t max_runlength,
                                       uint16_t max_threads){
    ///
    ///
    ///

    runlength_matrix template_matrix(boost::extents[2][4][max_runlength+1][max_runlength+1]);   // 0 length included
    vector<runlength_matrix> matrices_per_thread(max_threads, template_matrix);

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(parse_aligned_fasta,
                                        ref(bam_path),
                                        ref(reads_fasta_path),
                                        ref(read_indexes),
                                        ref(ref_runlength_sequences),
                                        ref(regions),
                                        ref(matrices_per_thread[i]),
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

    cerr << "Summing matrices from " << max_threads << " threads...\n";

    runlength_matrix matrix_sum = sum_matrices(matrices_per_thread);

    return matrix_sum;
}


template <typename T> void measure_runlength_distribution_from_coverage_data(path input_directory,
                                                       path reference_fasta_path,
                                                       path output_directory,
                                                       uint16_t max_runlength,
                                                       uint16_t max_threads,
                                                       path bed_path){

    cerr << "Using " + to_string(max_threads) + " threads\n";

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;

    T reader = T(input_directory);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Runlength encode the reference (store in memory)
    unordered_map<string,RunlengthSequenceElement> ref_runlength_sequences;
    path reference_fasta_path_rle;
    bool store_in_memory = true;
    reference_fasta_path_rle = runlength_encode_fasta_file(reference_fasta_path,
            ref_runlength_sequences,
            output_directory,
            store_in_memory,
            max_threads);

    reader.index();
    unordered_map<string,path> read_paths = reader.get_index();

    // Extract read names from index
    vector<string> read_names;
    for (auto& element: read_paths){
        read_names.push_back(element.first);
    }

    // Extract the sequences from MarginPolish TSVs (ignore coverage data for now)
    path reads_fasta_path_rle;
    reads_fasta_path_rle = write_all_consensus_sequences_to_fasta(reader,
            read_names,
            read_paths,
            input_directory,
            output_directory,
            max_threads);

    // Index all the files in the MP directory
    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    uint16_t k = 19;
    string minimap_preset = "asm20";    //TODO: make command line argument?
    bool explicit_mismatch = true;

    // Align coverage segments to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle,
            reads_fasta_path_rle,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    //TODO: cut off 50bp overlaps on MarginPolish TSVs?


    // If a BED file was provided, only iterate the regions of the BAM that may be found in the reference provided.
    // Otherwise, iterate the entire BAM.
    vector<Region> regions;
    if (bed_path.empty()){
        // Chunk alignment regions
        chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);
    }
    else{
        // Load the BED regions, and then find the union of the regions in BED and reference FASTA
        set<string> names;
        for (auto& item: ref_runlength_sequences){
            names.insert(item.first);
        }

        BedReader bed_reader(bed_path);
        bed_reader.read_regions(regions);
        bed_reader.subset_by_regions_name(regions, names);
    }

    cerr << "Iterating alignments...\n" << std::flush;

    // Launch threads for parsing alignments and generating matrices
    runlength_matrix matrix = get_runlength_matrix<T>(bam_path,
            input_directory,
            read_paths,
            ref_runlength_sequences,
            regions,
            max_runlength,
            max_threads);

    cerr << '\n';

    // Write output
    write_matrix_to_file(output_directory, matrix);
}


void measure_runlength_distribution_from_runnie(path runnie_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads){

    cerr << "Using " + to_string(max_threads) + " threads\n";

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;

    RunnieReader runnie_reader = RunnieReader(runnie_directory);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Runlength encode the reference (store in memory)
    unordered_map<string,RunlengthSequenceElement> ref_runlength_sequences;
    path reference_fasta_path_rle;
    bool store_in_memory = true;
    reference_fasta_path_rle = runlength_encode_fasta_file(reference_fasta_path,
            ref_runlength_sequences,
            output_directory,
            store_in_memory,
            max_threads);

    runnie_reader.index();
    unordered_map<string,RunnieIndex> read_indexes = runnie_reader.get_index();

    // Extract read names from index
    vector<string> read_names;
    for (auto& element: read_indexes){
        read_names.push_back(element.first);
    }

    // Extract the sequences from MarginPolish TSVs (ignore coverage data for now)
    path reads_fasta_path_rle;
    reads_fasta_path_rle = write_all_runnie_sequences_to_fasta(runnie_reader,
            read_names,
            read_indexes,
            runnie_directory,
            output_directory,
            max_threads);

    // Index all the files in the MP directory
    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    uint16_t k = 19;
    string minimap_preset = "asm20";    //TODO: make command line argument?
    bool explicit_mismatch = true;

    // Align Marginpolish reads to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle,
            reads_fasta_path_rle,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    //TODO: cut off 50bp overlaps on MarginPolish TSVs?

    // Chunk alignment regions
    vector<Region> regions;
    chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);

    cerr << "Iterating alignments...\n" << std::flush;

    // Launch threads for parsing alignments and generating matrices
    runlength_matrix matrix = get_runnie_runlength_matrix(bam_path,
            runnie_directory,
            read_indexes,
            ref_runlength_sequences,
            regions,
            max_runlength,
            max_threads);

    cerr << '\n';

    // Write output
    write_matrix_to_file(output_directory, matrix);
}


void measure_runlength_distribution_from_fasta(path reads_fasta_path,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads){

    cerr << "Using " + to_string(max_threads) + " threads\n";

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;

    // Initialize readers
    FastaReader reads_fasta_reader = FastaReader(reads_fasta_path);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Flag that decides whether RLE sequences should be added to a hash map in memory
    bool store_in_memory;

    // Runlength encode the REFERENCE, rewrite to another FASTA, and store in memory
    unordered_map<string,RunlengthSequenceElement> ref_runlength_sequences;
    path reference_fasta_path_rle;
    store_in_memory = true;
    reference_fasta_path_rle = runlength_encode_fasta_file(reference_fasta_path,
            ref_runlength_sequences,
            output_directory,
            store_in_memory,
            max_threads);

    // Runlength encode the READ SEQUENCES, rewrite to another FASTA, and DON'T store in memory
    unordered_map<string, RunlengthSequenceElement> _;
    path reads_fasta_path_rle;
    store_in_memory = false;
    reads_fasta_path_rle = runlength_encode_fasta_file(reads_fasta_path,
            _,
            output_directory,
            store_in_memory,
            max_threads);

    // Index all the files in the MP directory
    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    uint16_t k = 19;
    string minimap_preset = "asm20";    //TODO: make command line argument?
    bool explicit_mismatch = true;

    // Align reads to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle,
            reads_fasta_path_rle,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    // Chunk alignment regions
    vector<Region> regions;
    chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);

    cerr << "Iterating alignments...\n" << std::flush;

    reads_fasta_reader.index();
    unordered_map<string,FastaIndex> read_indexes = reads_fasta_reader.get_index();

    // Launch threads for parsing alignments and generating matrices
    runlength_matrix matrix = get_fasta_runlength_matrix(bam_path,
            reads_fasta_path,
            read_indexes,
            ref_runlength_sequences,
            regions,
            max_runlength,
            max_threads);

    cerr << '\n';

    // Write output
    write_matrix_to_file(output_directory, matrix);
}


template <typename T> void label_coverage_data(
        path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_threads,
        path bed_path){

    cerr << "Using " + to_string(max_threads) + " threads\n";

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;

    T reader = T(input_directory);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Runlength encode the reference (store in memory)
    unordered_map<string,RunlengthSequenceElement> ref_runlength_sequences;
    path reference_fasta_path_rle;
    bool store_in_memory = true;
    reference_fasta_path_rle = runlength_encode_fasta_file(reference_fasta_path,
            ref_runlength_sequences,
            output_directory,
            store_in_memory,
            max_threads);

    reader.index();
    unordered_map<string,path> read_paths = reader.get_index();

    // Extract read names from index
    vector<string> read_names;
    for (auto& element: read_paths){
        read_names.push_back(element.first);
    }

    // Extract the sequences from MarginPolish TSVs (ignore coverage data for now)
    path reads_fasta_path_rle;
    reads_fasta_path_rle = write_all_consensus_sequences_to_fasta(reader,
            read_names,
            read_paths,
            input_directory,
            output_directory,
            max_threads);

    // Index all the files in the MP directory
    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    uint16_t k = 19;
    string minimap_preset = "asm20";    //TODO: make command line argument?
    bool explicit_mismatch = true;

    // Align coverage segments to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle,
            reads_fasta_path_rle,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    //TODO: cut off 50bp overlaps on MarginPolish TSVs?


    // If a BED file was provided, only iterate the regions of the BAM that may be found in the reference provided.
    // Otherwise, iterate the entire BAM.
    vector<Region> regions;
    if (bed_path.empty()){
        // Chunk alignment regions
        chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);
    }
    else{
        // Load the BED regions, and then find the union of the regions in BED and reference FASTA
        set<string> names;
        for (auto& item: ref_runlength_sequences){
            names.insert(item.first);
        }

        BedReader bed_reader(bed_path);
        bed_reader.read_regions(regions);
        bed_reader.subset_by_regions_name(regions, names);
    }

    cerr << "Iterating alignments...\n" << std::flush;

    // Launch threads for parsing alignments and generating matrices
    get_coverage_labels<T>(bam_path,
            input_directory,
            output_directory,
            read_paths,
            ref_runlength_sequences,
            regions,
            max_threads);

    cerr << '\n';
}


void measure_runlength_distribution_from_marginpolish(path input_directory,
                                                       path reference_fasta_path,
                                                       path output_directory,
                                                       uint16_t max_runlength,
                                                       uint16_t max_threads,
                                                       path bed_path) {

    measure_runlength_distribution_from_coverage_data<MarginPolishReader>(
            input_directory,
            reference_fasta_path,
            output_directory,
            max_runlength,
            max_threads,
            bed_path);
}


void measure_runlength_distribution_from_shasta(
        path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads,
        path bed_path) {

    measure_runlength_distribution_from_coverage_data<ShastaReader>(
            input_directory,
            reference_fasta_path,
            output_directory,
            max_runlength,
            max_threads,
            bed_path);
}


void label_coverage_data_from_shasta(
        path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_threads,
        path bed_path) {

    label_coverage_data<ShastaReader>(
            input_directory,
            reference_fasta_path,
            output_directory,
            max_threads,
            bed_path);
}
