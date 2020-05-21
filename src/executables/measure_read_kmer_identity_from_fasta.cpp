#include "PileupGenerator.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "CigarKmer.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include "boost/program_options.hpp"
#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <thread>
#include <mutex>
#include <atomic>

using std::cerr;
using std::cout;
using std::string;
using std::ifstream;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::ref;
using std::move;
using std::atomic;
using std::atomic_fetch_add;

using std::experimental::filesystem::path;
using boost::program_options::value;
using boost::program_options::variables_map;
using boost::program_options::options_description;


void chunk_sequences_into_regions(vector<Region>& regions, unordered_map<string,SequenceElement>& sequences, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& [name, item]: sequences){
        chunk_sequence(regions, name, chunk_size, item.sequence.size());
    }
}


void iterate_read_kmers(
        path bam_path,
        path absolute_fasta_ref_path,
        path absolute_fasta_reads_path,
        uint8_t k,
        KmerIdentities& kmer_identities,
        vector <Region>& regions,
        atomic <uint64_t>& job_index){
    ///
    ///
    ///
    // Initialize readers and relevant containers
    FastaReader sequence_reader = FastaReader(absolute_fasta_reads_path);
    FastaReader ref_reader = FastaReader(absolute_fasta_ref_path);
    BamReader bam_reader = BamReader(bam_path);

    auto read_sequence = sequence_reader.generate_sequence_container();
    auto ref_sequence = ref_reader.generate_sequence_container();
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

    bool filter_secondary = true;
    bool filter_supplementary = false;
    uint16_t map_quality_cutoff = 5;

    // Volatiles
    bool in_left_bound;
    bool in_right_bound;
    string true_base;
    char observed_base;
    CigarKmer cigar_kmer(k);

    uint8_t match_code = Cigar::cigar_code_key.at("=");
    uint8_t mismatch_code = Cigar::cigar_code_key.at("X");
    uint8_t insert_code = Cigar::cigar_code_key.at("I");
    uint8_t delete_code = Cigar::cigar_code_key.at("D");

    // Allow all cigar operations //TODO: change this to an ordered set? or remove entirely?
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
        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary, filter_supplementary)) {
            // Reset kmer to empty deque
            cigar_kmer = CigarKmer(k);
            sequence_reader.get_sequence(read_sequence, aligned_segment.read_name);
            ref_reader.get_sequence(ref_sequence, aligned_segment.ref_name);

            // Iterate cigars that match the criteria
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
                // Skip any large indels and reset the kmer
                if (cigar.code != match_code and cigar.length > 30){
                    cigar_kmer = CigarKmer(k);
                    continue;
                }

                in_left_bound = (int64_t(region.start) <= coordinate.ref_index);
                in_right_bound = (coordinate.ref_index - 1 < int64_t(region.stop));

                // Subset alignment to portions of the read that are within the window/region
                if (in_left_bound and in_right_bound) {
                    observed_base = read_sequence.sequence[coordinate.read_true_index];
                    cigar_kmer.update(cigar, observed_base);
                    kmer_identities.update(cigar_kmer);
                }
            }

            i++;
        }

        cerr << "\33[2K\rParsed: " << region.to_string() << flush;
    }
}


void get_kmer_identity(path bam_path,
        path reference_fasta_path_rle,
        path reads_fasta_path_rle,
        uint8_t k,
        vector <Region>& regions,
        uint16_t max_threads,
        path output_directory){
    ///
    ///
    ///

    vector<KmerIdentities> kmer_identities_per_thread(max_threads, k);

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            threads.emplace_back(thread(iterate_read_kmers,
                                        ref(bam_path),
                                        ref(reference_fasta_path_rle),
                                        ref(reads_fasta_path_rle),
                                        k,
                                        ref(kmer_identities_per_thread[i]),
                                        ref(regions),
                                        ref(job_index)));

            // Call thread safe function to read and write to file
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

    // Sum together the kmer identity/cigar stats from all the threads
    KmerIdentities sum(k);
    for (auto& kmer_identities: kmer_identities_per_thread){
        sum += kmer_identities;
    }

    // Prepare output file
    path output_path = output_directory / "kmer_stats.csv";
    ofstream output_file(output_path);
    if (not output_file.is_open()){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    cerr << "WRITING FILE: " << absolute(output_path) << '\n';

    // Print kmer identities to file
    sum.write_to_output(output_file);
}


void load_fasta_sequences_from_file(path& fasta_path,
        vector <pair <string, FastaIndex> >& read_index_vector,
        unordered_map<string, SequenceElement>& sequences,
        mutex& map_mutex,
        atomic<uint64_t>& job_index){

    while (job_index < read_index_vector.size()) {
        // Fetch add
        uint64_t thread_job_index = job_index.fetch_add(1);

        // Initialize containers
        SequenceElement sequence;

        // Fetch Fasta sequence
        FastaReader fasta_reader = FastaReader(fasta_path);

        // Get read name and index
        string read_name = read_index_vector[thread_job_index].first;
        uint64_t read_index = read_index_vector[thread_job_index].second.byte_index;

        // First of each element is read_name, second is its index
        fasta_reader.get_sequence(sequence, read_name, read_index);

        // Append the sequence to a map of names:sequence
        map_mutex.lock();
        sequences[sequence.name] = move(sequence);
        map_mutex.unlock();

        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << read_name << flush;
    }
}


void load_fasta_file(path input_file_path,
        unordered_map <string, SequenceElement>& sequences,
        uint16_t max_threads) {

    // This reader is used to fetch an index
    FastaReader fasta_reader(input_file_path);

    // Get index
    unordered_map <string, FastaIndex> read_indexes;
    fasta_reader.get_indexes_mapped_by_name(read_indexes);

    // Convert the map object into an indexable object
    vector <pair <string, FastaIndex> > read_index_vector;
    get_vector_from_index_map(read_index_vector, read_indexes);

    mutex map_mutex;

    atomic<uint64_t> job_index = 0;
    vector<thread> threads;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        // Get data to send to threads (must not be sent by reference, or will lead to concurrency issues)
        try {
            // Call thread safe function to RL encode and write to file
            threads.emplace_back(thread(load_fasta_sequences_from_file,
                                        ref(input_file_path),
                                        ref(read_index_vector),
                                        ref(sequences),
                                        ref(map_mutex),
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


void measure_kmer_accuracy_from_fasta(
        path reference_fasta_path,
        path reads_fasta_path,
        path output_directory,
        string minimap_preset,
        uint16_t minimap_k,
        uint16_t max_threads,
        uint8_t k){

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 50*1000;

    // Index any input files before threading
    FastaReader reads_fasta_reader = FastaReader(reads_fasta_path);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);
    reads_fasta_reader.index();
    ref_fasta_reader.index();

    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    bool explicit_mismatch = true;

    // Align reads to the reference
    path bam_path;
    bam_path = align(
            reference_fasta_path,
            reads_fasta_path,
            output_directory,
            sort,
            index,
            delete_intermediates,
            minimap_k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    // Store ref sequences in memory
    unordered_map<string,SequenceElement> ref_sequences;
    load_fasta_file(reference_fasta_path, ref_sequences, max_threads);

    // Chunk alignment regions
    vector<Region> regions;
    chunk_sequences_into_regions(regions, ref_sequences, chunk_size);

    cerr << "Iterating alignments...\n";

    get_kmer_identity(bam_path,
        reference_fasta_path,
        reads_fasta_path,
        k,
        regions,
        max_threads,
        output_directory);
}



int main(int argc, char* argv[]){
    path ref_fasta_path;
    path reads_fasta_path;
    path output_dir;
    string minimap_preset;
    uint16_t max_threads;
    uint16_t minimap_k;
    uint16_t k;

    options_description options("Arguments");

    options.add_options()
        ("ref",
        value<path>(&ref_fasta_path),
        "File path of FASTA file containing REFERENCE sequences to be Run-length encoded")

        ("sequences",
        value<path>(&reads_fasta_path),
        "File path of FASTA file containing QUERY sequences to be Run-length encoded")

        ("minimap_preset",
        value<string>(&minimap_preset),
        "Minimap preset to be used, e.g. asm20, map-ont")

        ("minimap_k",
        value<uint16_t>(&minimap_k)->
        default_value(0),
        "Kmer seed size for minimap (overrides preset)")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. Files will be named based on input file name")

        ("max_threads",
        value<uint16_t>(&max_threads)->
        default_value(1),
        "Maximum number of threads to launch")

        ("k",
        value<uint16_t>(&k)->
        default_value(6),
        "kmer size to evaluate. Default = 6");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    cout << int(k) << '\n';

    measure_kmer_accuracy_from_fasta(
            ref_fasta_path,
            reads_fasta_path,
            output_dir,
            minimap_preset,
            minimap_k,
            max_threads,
            k);

    return 0;
}