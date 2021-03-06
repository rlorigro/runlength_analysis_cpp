#include "BinaryRunnieReader.hpp"
#include "BinaryRunnieWriter.hpp"
#include "PileupGenerator.hpp"
#include "RunlengthReader.hpp"
#include "RunlengthWriter.hpp"
#include "SequenceElement.hpp"
#include "RunnieReader.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Identity.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include <iostream>
#include <SimpleBayesianRunnieConsensusCaller.hpp>
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
        path runnie_reads_directory,
        path runlength_fasta_ref_path,
        path runlength_fasta_reads_path,
        path runlength_ref_path,
        path runlength_reads_path,
        path output_directory,
        uint32_t max_threads){

    SequenceElement sequence;
    RunlengthSequenceElement runlength_sequence;
    RunnieSequenceElement runnie_sequence;

    FastaReader ref_fasta_reader(fasta_ref_path);
    ref_fasta_reader.index();

    RunnieReader runnie_reads_reader(runnie_reads_directory);
    runnie_reads_reader.index();

    FastaWriter ref_fasta_writer(runlength_fasta_ref_path);
    FastaWriter reads_fasta_writer(runlength_fasta_reads_path);

    RunlengthWriter ref_runlength_writer(runlength_ref_path);
    BinaryRunnieWriter reads_runnie_writer(runlength_reads_path);

    cerr << "Writing FASTAs as runlength files...\n";

    while (not ref_fasta_reader.end_of_file) {
        ref_fasta_reader.next_element(sequence);
        runlength_encode(runlength_sequence, sequence);
        ref_fasta_writer.write(runlength_sequence);
        ref_runlength_writer.write_sequence(runlength_sequence);
    }
    ref_runlength_writer.write_indexes();

    for (auto& item: runnie_reads_reader.read_indexes){
        string name = item.first;
        runnie_reads_reader.fetch_sequence(runnie_sequence, name);
        reads_fasta_writer.write(runnie_sequence);
        reads_runnie_writer.write_sequence(runnie_sequence);
    }
    reads_runnie_writer.write_indexes();

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


void predict_consensus(Pileup& pileup, Pileup& ref_pileup, SimpleBayesianRunnieConsensusCaller& consensus_caller, Region& region, ofstream& output_file, mutex& file_write_mutex){
    vector<float> consensus;
    string consensus_sequence;

    for (size_t width_index = 0; width_index<pileup.pileup.size(); width_index++) {
        // Call standard alignment columns
        consensus_caller(pileup.pileup[width_index], consensus);
        append_consensus_sequence(consensus_sequence, consensus);

        if (consensus[0] != ref_pileup.pileup[width_index][0][Pileup::BASE]){
            cout << "BASE ERROR\n";
            cout << consensus[0] << " " << ref_pileup.pileup[width_index][0][Pileup::BASE] << '\n';
        }
        if (consensus[1] != ref_pileup.pileup[width_index][0][Pileup::REVERSAL+1]){
            cout << "LENGTH ERROR\n";
            cout << consensus[1] << " " << ref_pileup.pileup[width_index][0][Pileup::REVERSAL+1] << '\n';
            cout << region.name << ": " << region.start + width_index << '\n';

            size_t min_index;
            size_t max_index;

            min_index = max(size_t(width_index-1), size_t(0));
            max_index = min(size_t(width_index+1), size_t(pileup.pileup.size())-1);
            PileupGenerator::print(ref_pileup, min_index, max_index);
            PileupGenerator::print(pileup, min_index, max_index);

            min_index = max(size_t(width_index-10), size_t(0));
            max_index = min(size_t(width_index+10), size_t(pileup.pileup.size())-1);
            PileupGenerator::print(ref_pileup, min_index, max_index);
            PileupGenerator::print(pileup, min_index, max_index);
        }

        // Call insert columns if they exist
        if (pileup.inserts.count(width_index) > 0) {
            for (auto &column: pileup.inserts.at(width_index)) {
                // Inserts are always wrong if called as anything other than (gap,0)
                consensus_caller(column, consensus);
                append_consensus_sequence(consensus_sequence, consensus);
            }
        }
    }

    // Add fasta formatting to sequence string
    if (not consensus_sequence.empty()) {
        consensus_sequence = ">" + region.to_string() + "\n" + consensus_sequence + "\n";
    }

    file_write_mutex.lock();
    output_file << consensus_sequence;
    file_write_mutex.unlock();
}


void predict_chunk_consensus(path& bam_path,
        RunlengthReader& ref_runlength_reader,
        BinaryRunnieReader& reads_runnie_reader,
        vector<Region>& regions,
        ofstream& output_file,
        mutex& file_write_mutex,
        atomic<size_t>& job_index){

    PileupGenerator pileup_generator = PileupGenerator(bam_path, 80);

    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path config_path = project_directory / "config/SimpleBayesianConsensusCaller-6-runnie-raw-reads-5mb-chr11.csv";
    SimpleBayesianRunnieConsensusCaller consensus_caller(config_path);

    Pileup pileup;
    Pileup ref_pileup;

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);

        pileup_generator.fetch_region(regions[thread_job_index], reads_runnie_reader, pileup);
        pileup_generator.generate_reference_pileup(pileup, ref_pileup, regions[thread_job_index], ref_runlength_reader);
        predict_consensus(pileup, ref_pileup, consensus_caller, regions[thread_job_index], output_file, file_write_mutex);
    }
}


void get_consensus(path bam_path,
        path runlength_ref_path,
        path runnie_reads_path,
        path output_file_path,
        uint16_t max_threads){

    ofstream output_file(output_file_path);

    if (not output_file.is_open()){
        throw runtime_error("ERROR: output file could not be written to or created: " + output_file_path.string());
    }

    RunlengthReader ref_runlength_reader(runlength_ref_path);
    BinaryRunnieReader reads_runnie_reader(runnie_reads_path);

//    vector<Region> regions;
//    uint64_t chunk_size = 100*1000;
//    chunk_sequences_into_regions(regions, ref_runlength_reader.indexes, chunk_size);

    vector<Region> regions = {Region("hg38_dna", 2*1000*1000, 2*1000*1000+50*1000)}; //TODO: UNDO THIS TEST

    vector<thread> threads;
    atomic<size_t> job_index = 0;
    mutex file_write_mutex;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(predict_chunk_consensus,
                    ref(bam_path),
                    ref(ref_runlength_reader),
                    ref(reads_runnie_reader),
                    ref(regions),
                    ref(output_file),
                    ref(file_write_mutex),
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


void test(path fasta_ref_path, path runnie_reads_directory, path output_directory, uint32_t max_threads) {
    create_directories(output_directory);

    path reads_output_prefix = runnie_reads_directory;
    if (reads_output_prefix.stem() == ".") {
        reads_output_prefix = reads_output_prefix.parent_path().stem();
    }
    else {
        reads_output_prefix = reads_output_prefix.stem();
    }

    // Create filenames for runlength files
    path runlength_ref_path = fasta_ref_path;
    runlength_ref_path.replace_extension(".rnq");
    path runnie_reads_path = runnie_reads_directory;
    runnie_reads_path.replace_extension(".rnq");

    // Create filenames for runlength FASTA files.........
    path runlength_fasta_ref_path = fasta_ref_path;
    runlength_fasta_ref_path.replace_extension("rle.fasta");
    path runlength_fasta_reads_path = runnie_reads_directory / path(reads_output_prefix.string() + ".rle.fasta");

    path bam_path = align_as_runlength(
            fasta_ref_path,
            runnie_reads_directory,
            runlength_fasta_ref_path,
            runlength_fasta_reads_path,
            runlength_ref_path,
            runnie_reads_path,
            output_directory,
            max_threads);

    path output_file_path = output_directory / "consensus.fasta";

    get_consensus(bam_path,
            runlength_ref_path,
            runnie_reads_path,
            output_file_path,
            max_threads);

    measure_identity_from_fasta(
            output_file_path,
            fasta_ref_path,
            output_directory,
            "asm20",
            max_threads);
}


int main(int argc, char* argv[]){
//    path script_path = __FILE__;
//    path project_directory = script_path.parent_path().parent_path().parent_path();
//
//    path ref_fasta_path = "/home/ryan/data/Nanopore/Human/paolo/5mb_test_region/reference.fa";
//    path reads_directory = "/home/ryan/data/Nanopore/Human/paolo/5mb_test_region/runnie/out";
//    uint32_t max_threads = 30;

    path ref_fasta_path;
    path reads_directory;
    path output_dir;
    uint16_t max_threads;

    options_description options("Arguments");

    options.add_options()
        ("ref",
        value<path>(&ref_fasta_path),
        "File path of FASTA file containing REFERENCE sequences to be Run-length encoded")

        ("sequences",
        value<path>(&reads_directory),
        "Path of DIRECTORY containing runnie .OUT QUERY sequences to be aligned")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. File will be named based on input file name")

        ("max_threads",
        value<uint16_t>(&max_threads)->
        default_value(1),
        "Maximum number of threads to launch");

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
            reads_directory,
            output_dir,
            max_threads);

    return 0;
}
