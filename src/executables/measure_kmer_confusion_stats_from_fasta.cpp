#include "PileupGenerator.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "PileupKmer.hpp"
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



void chunk_sequences_into_regions(vector<Region>& regions, unordered_map <string, FastaIndex>& read_indexes, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& [name, item]: read_indexes){
        chunk_sequence(regions, name, chunk_size, item.length);
    }
}


void iterate_pileup_kmers(path bam_path,
        path absolute_fasta_ref_path,
        path absolute_fasta_reads_path,
        vector <Region>& regions,
        size_t window_size,
        uint8_t k,
        KmerConfusionStats& kmer_confusion_stats,
        bool reference_based,
        atomic <uint64_t>& job_index){

    path absolute_bam_path = absolute(bam_path);

    Pileup ref_pileup;
    Pileup read_pileup;
    Region region;

    FastaReader ref_reader = FastaReader(absolute_fasta_ref_path);
    FastaReader sequence_reader = FastaReader(absolute_fasta_reads_path);

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        PileupGenerator pileup_generator = PileupGenerator(absolute_bam_path, 60);

        pileup_generator.fetch_region(region, sequence_reader, read_pileup);
        pileup_generator.generate_reference_pileup(read_pileup, ref_pileup, region, ref_reader);

        PileupKmerIterator ref_pileup_iterator(ref_pileup, window_size, k);
        PileupKmerIterator read_pileup_iterator(read_pileup, window_size, k);

        for (size_t i = 0; i < read_pileup.pileup.size() - 1; i++) {
            ref_pileup_iterator.step(ref_pileup);
            read_pileup_iterator.step(read_pileup);

            if (reference_based) {
                read_pileup_iterator.update_ref_kmer_confusion_stats(ref_pileup_iterator, kmer_confusion_stats);
            }
            else{
                read_pileup_iterator.update_read_kmer_confusion_stats(ref_pileup_iterator, kmer_confusion_stats);
            }
        }

        cerr << "\33[2K\rParsed: " << region.to_string() << flush;
    }
}


void get_kmer_stats(path bam_path,
        path reference_fasta_path_rle,
        path reads_fasta_path_rle,
        size_t window_size,
        uint8_t k,
        vector <Region>& regions,
        uint16_t max_threads,
        bool reference_based,
        path output_directory){
    ///
    ///
    ///

    vector<KmerConfusionStats> kmer_stats_per_thread(max_threads, k);

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(iterate_pileup_kmers,
                                        ref(bam_path),
                                        ref(reference_fasta_path_rle),
                                        ref(reads_fasta_path_rle),
                                        ref(regions),
                                        window_size,
                                        k,
                                        ref(kmer_stats_per_thread[i]),
                                        reference_based,
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

    // Prepare output file
    path output_path = output_directory / "kmer_identity.csv";
    ofstream output_file(output_path);
    if (not output_file.is_open()){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    KmerConfusionStats sum_of_stats;
    for (auto &stat: kmer_stats_per_thread) {
        sum_of_stats += stat;
    }

    output_file << sum_of_stats.to_string();
}


void measure_kmer_stats_from_fasta(path reference_fasta_path,
        path reads_fasta_path,
        path output_directory,
        string minimap_preset,
        uint16_t minimap_k,
        uint16_t max_threads,
        size_t window_size,
        uint8_t k,
        bool reference_based){

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 250*1000;

    // Initialize readers
    FastaReader reads_fasta_reader = FastaReader(reads_fasta_path);
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Ensure that the FASTAs are indexed before starting threading
    FastaReader ref_reader = FastaReader(reference_fasta_path);
    FastaReader sequence_reader = FastaReader(reads_fasta_path);
    ref_reader.index();
    sequence_reader.index();

    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    bool explicit_mismatch = true;

    // Align reads to the reference
    path bam_path;
    bam_path = align(reference_fasta_path,
            reads_fasta_path,
            output_directory,
            sort,
            index,
            delete_intermediates,
            minimap_k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    // Chunk alignment regions
    vector<Region> regions;
    unordered_map <string, FastaIndex > read_indexes;
    ref_reader.get_indexes_mapped_by_name(read_indexes);

    chunk_sequences_into_regions(regions, read_indexes, chunk_size);

    get_kmer_stats(bam_path,
        reference_fasta_path,
        reads_fasta_path,
        window_size,
        k,
        regions,
        max_threads,
        reference_based,
        output_directory);
}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path reads_fasta_path;
    path output_dir;
    string minimap_preset;
    uint16_t max_threads;
    size_t window_size;
    uint16_t minimap_k;
    uint16_t k;
    bool reference_based;

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

        ("window_size",
        value<size_t>(&window_size)->
        default_value(5),
        "Size of window in kmer space to consider when looking for matching kmers")

        ("k",
        value<uint16_t>(&k)->
        default_value(6),
        "kmer size to evaluate. Default = 6")

        ("ref_kmer",
        value<bool>(&reference_based)->
        default_value(false),
        "Whether or not to compute proportion of kmer matches on a reference basis or read basis");

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

    measure_kmer_stats_from_fasta(ref_fasta_path,
                                reads_fasta_path,
                                output_dir,
                                minimap_preset,
                                minimap_k,
                                max_threads,
                                window_size,
                                k,
                                reference_based
    );

    return 0;
}