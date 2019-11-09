#include "BinaryRunnieReader.hpp"
#include "PileupGenerator.hpp"
#include "RunlengthReader.hpp"
#include "RunlengthWriter.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "Identity.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include <iostream>
#include <SimpleBayesianRunnieConsensusCaller.hpp>

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
    ref_fasta_reader.index();

    cerr << "WRITING: " << absolute_runlength_reads_path << '\n';
    cerr << "WRITING: " << absolute_runlength_ref_path << '\n';
//    BinaryRunnieWriter reads_writer(absolute_runlength_reads_path);
    RunlengthWriter ref_writer(absolute_runlength_ref_path);


//    while (not reads_fasta_reader.end_of_file) {
//        reads_fasta_reader.next_element(sequence);
//        runlength_encode(runlength_sequence, sequence);
//        reads_writer.write_sequence(runlength_sequence);
//    }
//    reads_writer.write_indexes();

    while (not ref_fasta_reader.end_of_file) {
        ref_fasta_reader.next_element(sequence);
        runlength_encode(runlength_sequence, sequence);
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


void predict_consensus(Pileup& pileup, SimpleBayesianRunnieConsensusCaller& consensus_caller, Region& region, ofstream& output_file){
    vector<float> consensus;
    string consensus_sequence;

    for (size_t width_index = 0; width_index<pileup.pileup.size(); width_index++) {
        // Call standard alignment columns
        consensus_caller(pileup.pileup[width_index], consensus);
        append_consensus_sequence(consensus_sequence, consensus);
//        cout << float_to_base(consensus[0]) << consensus[1] << '\n';

        // Call insert columns if they exist
        if (pileup.inserts.count(width_index) > 0) {
            for (auto &column: pileup.inserts.at(width_index)) {
                consensus_caller(column, consensus);
                append_consensus_sequence(consensus_sequence, consensus);
            }
        }
    }

    cout << "CONSENSUS: " << consensus_sequence << '\n';

    // Add fasta formatting to sequence string
    if (not consensus_sequence.empty()) {
        consensus_sequence = ">" + region.to_string() + "\n" + consensus_sequence + "\n";
    }

    output_file << consensus_sequence;
}


void test(path absolute_fasta_ref_path, path absolute_runlength_ref_path, path absolute_runnie_reads_path, path bam_path, path output_directory) {
    create_directories(output_directory);

    path absolute_bam_path = absolute(bam_path);

    PileupGenerator pileup_generator = PileupGenerator(absolute_bam_path, 80);

    ifstream check(absolute_runnie_reads_path);

    BinaryRunnieReader reads_runnie_reader(absolute_runnie_reads_path);
    RunlengthReader ref_runlength_reader(absolute_runlength_ref_path);

    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path config_path = project_directory / "config/SimpleBayesianConsensusCaller-6-runnie-raw-reads-5mb-chr11.csv";
    SimpleBayesianRunnieConsensusCaller consensus_caller(config_path);

    Pileup pileup;

    path output_file_path = output_directory / "consensus.fasta";
    ofstream output_file(output_file_path);

    if (not output_file.is_open()){
        throw runtime_error("ERROR: output file could not be written to or created: " + output_file_path.string());
    }

    Region region("hg38_dna",2*1000*1000,2*1000*1000+2*1000);

    pileup_generator.fetch_region(region, ref_runlength_reader, reads_runnie_reader, pileup);
    pileup_generator.print(pileup);
    predict_consensus(pileup, consensus_caller, region, output_file);

    RunlengthSequenceElement ref_runlength_sequence;
    SequenceElement ref_sequence;
    string ref_name = "hg38_dna";

    ref_runlength_reader.get_sequence(ref_runlength_sequence, ref_name);
    runlength_decode(ref_runlength_sequence, ref_sequence, region.start, region.stop+1);

    cout << "REFERENCE: " << ref_sequence.sequence << '\n';

    measure_identity_from_fasta(
            output_file_path,
            absolute_fasta_ref_path,
            output_directory,
            30);

}


int main(int argc, char* argv[]){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    path bam_path = project_directory / "data/test/runnie_bayes/out_rle_VS_reference_rle.sorted.bam";
    path ref_fasta_path = project_directory / "data/test/runnie_bayes/reference.fa";
    path ref_runlength_path = project_directory / "data/test/runnie_bayes/reference.rlq";
    path reads_runnie_path = project_directory / "data/test/runnie_bayes/out.rnq";
    path output_dir = "output/minimap_consensus_runnie_small";

    test(
            ref_fasta_path,
            ref_runlength_path,
            reads_runnie_path,
            bam_path,
            output_dir);

    return 0;
}
