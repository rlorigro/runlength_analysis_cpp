#include "BinaryRunnieReader.hpp"
#include "BinaryRunnieWriter.hpp"
#include "PileupGenerator.hpp"
#include "RunlengthReader.hpp"
#include "RunlengthWriter.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Identity.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include <iostream>
#include <SimpleBayesianConsensusCaller.hpp>
#include <thread>
#include <mutex>
#include <atomic>
#include <exception>

using std::exception;
using std::thread;
using std::mutex;
using std::lock_guard;
using std::move;
using std::ref;
using std::atomic;
using std::atomic_fetch_add;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::min;
using std::max;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using std::experimental::filesystem::create_directories;


void chunk_sequences_into_regions(vector<Region>& regions, vector<RunlengthIndex> indexes, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& item: indexes){
        chunk_sequence(regions, item.name, chunk_size, item.sequence_length);
    }
}


path align_as_runlength(path fasta_ref_path,
        path fasta_reads_path,
        path runlength_fasta_ref_path,
        path runlength_fasta_reads_path,
        path runlength_ref_path,
        path runlength_reads_path,
        path output_directory,
        uint32_t max_threads){

    SequenceElement sequence;
    RunlengthSequenceElement runlength_sequence;

    FastaReader ref_fasta_reader(fasta_ref_path);
    ref_fasta_reader.index();

    FastaReader reads_fasta_reader(fasta_reads_path);
    reads_fasta_reader.index();

    FastaWriter ref_fasta_writer(runlength_fasta_ref_path);
    FastaWriter reads_fasta_writer(runlength_fasta_reads_path);

    RunlengthWriter ref_runlength_writer(runlength_ref_path);
    RunlengthWriter reads_runlength_writer(runlength_reads_path);

    cerr << "Writing FASTAs as runlength files...\n";

    while (not ref_fasta_reader.end_of_file) {
        ref_fasta_reader.next_element(sequence);
        runlength_encode(runlength_sequence, sequence);
        ref_fasta_writer.write(runlength_sequence);
        ref_runlength_writer.write_sequence(runlength_sequence);
    }
    ref_runlength_writer.write_indexes();

    while (not reads_fasta_reader.end_of_file) {
        reads_fasta_reader.next_element(sequence);
        runlength_encode(runlength_sequence, sequence);
        reads_fasta_writer.write(runlength_sequence);
        reads_runlength_writer.write_sequence(runlength_sequence);
    }
    reads_runlength_writer.write_indexes();

    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;
    uint16_t k = 19;
    string minimap_preset = "asm20";
    bool explicit_mismatch = true;

    // Align reads to the reference
    path bam_path;
    bam_path = align(
            runlength_fasta_ref_path,
            runlength_fasta_reads_path,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    path absolute_bam_path = absolute(bam_path);

    return absolute_bam_path;
}


void print_column(vector <vector <float> > column){
    string s;
    size_t i;
    for (auto& data_vector: column) {
        i = 0;
        for (auto& element: data_vector) {
            if (i==0) {
                s += float_to_base(element);
            }
            else{
                s += to_string(int(element));
            }
            i++;
        }
        s += " ";
        cout << s;
    }
}


void append_consensus_sequence(string& consensus_sequence, vector<float>& consensus){
    char base;
    size_t length;
    string consensus_string;

    base = float_to_base_char(consensus[0]);
    length = size_t(consensus[1]);

    if (length > 0) {
        consensus_sequence += string(length, base);
    }
}


void predict_consensus(PileupGenerator& pileup_generator,
        RunlengthReader& ref_runlength_reader,
        RunlengthReader& reads_runlength_reader,
        Pileup& pileup,
        Pileup& ref_pileup,
        SimpleBayesianConsensusCaller& consensus_caller,
        Region& region,
        vector<ofstream>& output_files,
        vector<mutex>& file_write_mutexes,
        uint32_t max_coverage){

    vector<float> consensus = {4,0};
    vector<string> consensus_sequences((max_coverage/5));

    for (size_t width_index = 0; width_index<pileup.pileup.size(); width_index++) {
        vector <vector <float> > column = pileup.pileup[width_index];

        for (int64_t i=(max_coverage/5)-1; i>=0; i--) {
            // Subset the coverage (iteratively shrinking the vector) IFF the column is bigger than the current value
            column.resize(min(int64_t(column.size()), (i+1) * 5));

            // Call consensus on ref columns
            consensus_caller(column, consensus);
            append_consensus_sequence(consensus_sequences[i], consensus);
        }

        // Call consensus on insert columns
        if (pileup.inserts.count(width_index) > 0) {
            for (auto &insert_column: pileup.inserts.at(width_index)) {
                column = insert_column;

                for (int64_t i = (max_coverage/5)-1; i >= 0; i--) {
                    column.resize(min(int64_t(column.size()), (i+1) * 5));

                    consensus_caller(column, consensus);
                    append_consensus_sequence(consensus_sequences[i], consensus);
                }
            }
        }
    }

    for (size_t i=0; i < (max_coverage/5); i++) {
        // Add fasta formatting to sequence string
        if (not consensus_sequences[i].empty()) {
            consensus_sequences[i] = ">" + region.to_string() + "_coverage-" + to_string((i+1)*5) + "\n" + consensus_sequences[i] + "\n";
        }

        file_write_mutexes[i].lock();
        output_files[i] << consensus_sequences[i];
        file_write_mutexes[i].unlock();
    }
}


void predict_chunk_consensus(path& bam_path,
        RunlengthReader& ref_runlength_reader,
        RunlengthReader& reads_runlength_reader,
        vector<Region>& regions,
        vector<ofstream>& output_files,
        vector<mutex>& file_write_mutex,
        uint16_t& max_coverage,
        atomic<size_t>& job_index){

    PileupGenerator pileup_generator = PileupGenerator(bam_path, max_coverage);

    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path config_path = project_directory / "config/SimpleBayesianConsensusCaller-5.csv";
    SimpleBayesianConsensusCaller consensus_caller(config_path);

    Pileup pileup;
    Pileup ref_pileup;

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);

        pileup_generator.fetch_region(regions[thread_job_index], reads_runlength_reader, pileup);
        pileup_generator.generate_reference_pileup(pileup, ref_pileup, regions[thread_job_index], ref_runlength_reader);
        predict_consensus(pileup_generator, ref_runlength_reader, reads_runlength_reader, pileup, ref_pileup, consensus_caller, regions[thread_job_index], output_files, file_write_mutex, max_coverage);
    }
}


void get_consensus(path bam_path,
        path runlength_ref_path,
        path runlength_reads_path,
        vector <ofstream>& output_files,
        uint16_t max_coverage,
        uint16_t max_threads){

    RunlengthReader ref_runlength_reader(runlength_ref_path);
    RunlengthReader reads_runlength_reader(runlength_reads_path);

    vector<Region> regions;
    uint64_t chunk_size = 100*1000;
    chunk_sequences_into_regions(regions, ref_runlength_reader.indexes, chunk_size);

//    vector<Region> regions = {Region("hg38_dna", 2*1000*1000, 2*1000*1000+400*1000)}; //TODO: UNDO THIS TEST
//    vector<Region> regions = {Region("CM010818.1", 74000000, 74000000+400*1000)}; //TODO: UNDO THIS TEST

    vector<thread> threads;
    atomic<size_t> job_index = 0;
    vector<mutex> file_write_mutexes(output_files.size());

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(predict_chunk_consensus,
                    ref(bam_path),
                    ref(ref_runlength_reader),
                    ref(reads_runlength_reader),
                    ref(regions),
                    ref(output_files),
                    ref(file_write_mutexes),
                    ref(max_coverage),
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
}


void test(path fasta_ref_path, path fasta_reads_path, path output_directory, uint16_t max_threads, uint16_t max_coverage) {
    create_directories(output_directory);

    // Create filenames for runlength files
    path runlength_ref_path = fasta_ref_path;
    runlength_ref_path.replace_extension(".rnq");
    path runlength_reads_path = fasta_reads_path;
    runlength_reads_path.replace_extension(".rnq");

    // Create filenames for runlength FASTA files.........
    path runlength_fasta_ref_path = fasta_ref_path;
    runlength_fasta_ref_path.replace_extension("rle.fasta");
    path runlength_fasta_reads_path = fasta_reads_path;
    runlength_fasta_reads_path.replace_extension("rle.fasta");

    path bam_path = align_as_runlength(
            fasta_ref_path,
            fasta_reads_path,
            runlength_fasta_ref_path,
            runlength_fasta_reads_path,
            runlength_ref_path,
            runlength_reads_path,
            output_directory,
            max_threads);

    path output_file_prefix = output_directory / "consensus";

    // Initialize a vector of fasta files to correspond to each coverage
    vector<ofstream> output_files;
    vector<path> output_file_paths;
    for (size_t i=1; i<size_t(max_coverage/5) + 1; i++) {
        path file_path = output_file_prefix.string() + "_coverage-" + to_string(i*5) + ".fasta";

        output_file_paths.emplace_back(file_path);
        output_files.emplace_back(file_path.string());

        if (not output_files.back().is_open()){
            throw runtime_error("ERROR: output file could not be written to or created: " + file_path.string());
        }
    }

    get_consensus(bam_path,
            runlength_ref_path,
            runlength_reads_path,
            output_files,
            max_coverage,
            max_threads);

    path results_path = output_directory / "results.txt";
    ofstream results_file(results_path);
    CigarStats stats;

    vector<Region> regions = {Region("chr1", 20000000, 24999999),
                              Region("chr1", 25000000, 29999999),
                              Region("chr1", 30000000, 34999999),
                              Region("chr1", 35000000, 40000000)
    };

    for (size_t i=0; i<output_file_paths.size(); i++) {
        stats = measure_identity_from_fasta(
                output_file_paths[i],
                fasta_ref_path,
                output_directory,
                "asm20",
                max_threads,
                false,
                0,
                regions);

        results_file << "Coverage " << (i+1)*5 << '\n';
        results_file << "identity (M/(M+X+I+D)):\t" << stats.calculate_identity() << '\n';
        results_file << stats.to_string() << '\n';
    }
}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path reads_fasta_path;
    path output_dir;
    uint16_t max_threads;
    uint16_t max_coverage;

    options_description options("Arguments");

    options.add_options()
        ("ref",
        value<path>(&ref_fasta_path),
        "File path of FASTA file containing REFERENCE sequences to be Run-length encoded")

        ("sequences",
        value<path>(&reads_fasta_path),
        "File path of FASTA file containing QUERY sequences to be Run-length encoded")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. File will be named based on input file name")

        ("max_threads",
        value<uint16_t>(&max_threads)->
        default_value(1),
        "Maximum number of threads to launch")

        ("max_coverage",
        value<uint16_t>(&max_coverage)->
        default_value(120),
        "Maximum coverage to test at. IF coverage is lower than this, it will still be used as is.");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    test(
            ref_fasta_path,
            reads_fasta_path,
            output_dir,
            max_threads,
            max_coverage);

    return 0;
}
