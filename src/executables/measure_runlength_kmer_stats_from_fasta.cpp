#include "FastaWriter.hpp"
#include "PileupGenerator.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "Align.hpp"
#include "Kmer.hpp"
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



void chunk_sequences_into_regions(vector<Region>& regions, unordered_map<string,RunlengthSequenceElement>& sequences, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& [name, item]: sequences){
        chunk_sequence(regions, name, chunk_size, item.sequence.size());
    }
}


void iterate_pileup_kmers(path bam_path,
        path absolute_fasta_ref_path,
        path absolute_fasta_reads_path,
        vector <Region>& regions,
        size_t window_size,
        uint8_t k,
        KmerStats& kmer_stats,
        atomic <uint64_t>& job_index){

    path absolute_bam_path = absolute(bam_path);

    Pileup pileup;
    Region region;

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        PileupGenerator pileup_generator = PileupGenerator(absolute_bam_path, 20);

        FastaReader ref_reader = FastaReader(absolute_fasta_ref_path);
        FastaReader sequence_reader = FastaReader(absolute_fasta_reads_path);

        pileup_generator.fetch_region(region, sequence_reader, pileup);

        PileupKmerIterator pileup_iterator(pileup, window_size, k);

        for (size_t i = 0; i < pileup.pileup.size() - 1; i++) {
            pileup_iterator.step(pileup);
        }

        kmer_stats += pileup_iterator.kmer_stats;
    }
}


KmerStats get_kmer_stats(path bam_path,
        path ref_fasta_path,
        path reads_fasta_path,
        size_t window_size,
        uint8_t k,
        vector <Region>& regions,
        uint16_t max_threads){
    ///
    ///
    ///

    vector<KmerStats> kmer_stats_per_thread;

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(iterate_pileup_kmers,
                                        ref(bam_path),
                                        ref(ref_fasta_path),
                                        ref(reads_fasta_path),
                                        ref(regions),
                                        window_size,
                                        k,
                                        ref(kmer_stats_per_thread[i]),
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

    KmerStats sum_of_stats(k);

    for (auto& stat: kmer_stats_per_thread){
        sum_of_stats += stat;
    }

    return sum_of_stats;
}


void measure_kmer_stats_from_fasta(path reference_fasta_path,
        path reads_fasta_path,
        path output_directory,
        string minimap_preset,
        uint16_t max_threads,
        size_t window_size,
        uint8_t k){

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
    store_in_memory = false;
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
    uint16_t k_align = 19;
    bool explicit_mismatch = true;

    // Align reads to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle,
            reads_fasta_path_rle,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k_align,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    // Chunk alignment regions
    vector<Region> regions;
    chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);

    KmerStats stats = get_kmer_stats(bam_path,
        reference_fasta_path,
        reads_fasta_path,
        window_size,
        k,
        regions,
        max_threads);

    cout << stats.to_string() << '\n';
}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path reads_fasta_path;
    path output_dir;
    string minimap_preset;
    uint16_t max_threads;
    size_t window_size;
    uint8_t k;

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
        value<uint8_t>(&k)->
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

        ("window_size",
        value<size_t>(&window_size)->
        default_value(5),
        "Size of window in kmer space to consider when looking for matching kmers")

        ("k",
        value<uint8_t>(&k)->
        default_value(6),
        "kmer size to evaluate");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    measure_kmer_stats_from_fasta(reads_fasta_path,
                                ref_fasta_path,
                                output_dir,
                                minimap_preset,
                                max_threads,
                                window_size,
                                k);

    return 0;
}