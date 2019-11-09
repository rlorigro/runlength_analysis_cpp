#include "boost/program_options.hpp"
#include "RunlengthReader.hpp"
#include "PileupGenerator.hpp"
#include "RunlengthWriter.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "Identity.hpp"
#include "FastaWriter.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include <iostream>
#include <SimpleBayesianConsensusCaller.hpp>

using std::cerr;
using std::ifstream;
using std::ofstream;
using std::cout;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using std::experimental::filesystem::create_directories;


void chunk_sequences_into_regions(vector<Region>& regions, vector<RunlengthSequenceElement>& sequences, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///

    // For every sequence
    for (auto& item: sequences){
        chunk_sequence(regions, item.name, chunk_size, item.sequence.size());
    }
}


void write_fasta_as_runlength(
        path& absolute_fasta_reads_path,
        path& absolute_fasta_ref_path,
        path& absolute_runlength_fasta_reads_path,
        path& absolute_runlength_fasta_ref_path,
        path& absolute_runlength_reads_path,
        path& absolute_runlength_ref_path,
        vector <RunlengthSequenceElement>& ref_runlength_sequences){

    SequenceElement sequence;
    RunlengthSequenceElement runlength_sequence;

    cerr << "READING: " << absolute_fasta_reads_path << '\n';
    cerr << "READING: " << absolute_fasta_ref_path << '\n';
    FastaReader reads_fasta_reader(absolute_fasta_reads_path);
    FastaReader ref_fasta_reader(absolute_fasta_ref_path);
    reads_fasta_reader.index();

    cerr << "WRITING: " << absolute_runlength_fasta_reads_path << '\n';
    cerr << "WRITING: " << absolute_runlength_fasta_ref_path << '\n';
    FastaWriter reads_fasta_writer(absolute_runlength_fasta_reads_path);
    FastaWriter ref_fasta_writer(absolute_runlength_fasta_ref_path);

    cerr << "WRITING: " << absolute_runlength_reads_path << '\n';
    cerr << "WRITING: " << absolute_runlength_ref_path << '\n';
    RunlengthWriter reads_writer(absolute_runlength_reads_path);
    RunlengthWriter ref_writer(absolute_runlength_ref_path);


    while (not reads_fasta_reader.end_of_file) {
        reads_fasta_reader.next_element(sequence);
        runlength_encode(runlength_sequence, sequence);
        reads_fasta_writer.write(runlength_sequence);
        reads_writer.write_sequence(runlength_sequence);
    }
    reads_writer.write_indexes();

    while (not ref_fasta_reader.end_of_file) {
        ref_fasta_reader.next_element(sequence);
        runlength_encode(runlength_sequence, sequence);
        ref_fasta_writer.write(runlength_sequence);
        ref_writer.write_sequence(runlength_sequence);

        ref_runlength_sequences.push_back(runlength_sequence);
    }
    ref_writer.write_indexes();
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


void predict_consensus(Pileup& pileup, SimpleBayesianConsensusCaller& consensus_caller, Region& region, ofstream& output_file){
    vector<float> consensus;
    string consensus_sequence;

    for (size_t width_index = 0; width_index<pileup.pileup.size(); width_index++) {
        // Call standard alignment columns
        print_column(pileup.pileup[width_index]);
        cout << '\n';
        consensus_caller(pileup.pileup[width_index], consensus);
        append_consensus_sequence(consensus_sequence, consensus);

        // Call insert columns if they exist
        if (pileup.inserts.count(width_index) > 0) {
            for (auto &column: pileup.inserts.at(width_index)) {
                consensus_caller(column, consensus);
                append_consensus_sequence(consensus_sequence, consensus);
            }
        }
    }

    // Add fasta formatting to sequence string
    if (not consensus_sequence.empty()) {
        consensus_sequence = ">" + region.to_string() + "\n" + consensus_sequence + "\n";
    }

    output_file << consensus_sequence;
}

void predict_runlength_from_minimap(path absolute_fasta_ref_path, path absolute_fasta_reads_path, path output_directory, uint16_t max_threads) {
    // Create filenames for runlength files
    path absolute_runlength_reads_path = absolute_fasta_reads_path;
    absolute_runlength_reads_path.replace_extension(".rlq");
    path absolute_runlength_ref_path = absolute_fasta_ref_path;
    absolute_runlength_ref_path.replace_extension(".rlq");

    // Create filenames for runlength FASTA files.........
    path absolute_runlength_fasta_reads_path = absolute_fasta_reads_path;
    absolute_runlength_fasta_reads_path.replace_extension("rle.fasta");
    path absolute_runlength_fasta_ref_path = absolute_fasta_ref_path;
    absolute_runlength_fasta_ref_path.replace_extension("rle.fasta");

    vector<RunlengthSequenceElement> ref_runlength_sequences;
    write_fasta_as_runlength(
        absolute_fasta_reads_path,
        absolute_fasta_ref_path,
        absolute_runlength_fasta_reads_path,
        absolute_runlength_fasta_ref_path,
        absolute_runlength_reads_path,
        absolute_runlength_ref_path,
        ref_runlength_sequences);

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
            absolute_runlength_fasta_ref_path,
            absolute_runlength_fasta_reads_path,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    path absolute_bam_path = absolute(bam_path);

    PileupGenerator pileup_generator = PileupGenerator(absolute_bam_path, 80);

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;  // TODO: un-hardcode this

    // Chunk alignment regions
    vector<Region> regions;
    chunk_sequences_into_regions(regions, ref_runlength_sequences, chunk_size);

    ifstream check(absolute_runlength_reads_path);
    cerr << check.good() << '\n';

    RunlengthReader reads_runlength_reader(absolute_runlength_reads_path);
    RunlengthReader ref_runlength_reader(absolute_runlength_ref_path);

    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path config_path = project_directory / "config/SimpleBayesianConsensusCaller-4.csv";
    SimpleBayesianConsensusCaller consensus_caller(config_path);

    Pileup pileup;

    path output_file_path = output_directory / "consensus.fasta";
    ofstream output_file(output_file_path);

    if (not output_file.is_open()){
        throw runtime_error("ERROR: output file could not be written to or created: " + output_file_path.string());
    }

    for (auto& region: regions) {
        output_file << ">" << region.to_string() << '\n';
        pileup_generator.fetch_region(region, ref_runlength_reader, reads_runlength_reader, pileup);
        predict_consensus(pileup, consensus_caller, region, output_file);
    }

    measure_identity_from_fasta(
            output_file_path,
            absolute_fasta_ref_path,
            output_directory,
            max_threads);

}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path reads_fasta_path;
    path output_dir;
    uint16_t max_threads;

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

    predict_runlength_from_minimap(
            ref_fasta_path,
            reads_fasta_path,
            output_dir,
            max_threads);

    return 0;
}
