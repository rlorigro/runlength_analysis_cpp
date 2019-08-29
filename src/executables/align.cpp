#include <stdexcept>
#include "Miscellaneous.hpp"
#include "Align.hpp"
#include "boost/program_options.hpp"

using std::runtime_error;
using boost::program_options::options_description;


int main(int argc, char* argv[]){
    path ref_sequence_path;
    path read_sequence_path;
    path output_dir;
    string minimap_preset;
    bool explicit_mismatch;
    uint16_t max_threads;
    uint16_t k;
    bool sort;
    bool index;
    bool keep;

    options_description options("Arguments:");

    options.add_options()
        ("sequences",
        value<path>(&read_sequence_path),
        "File path of FASTA file containing query sequences to be aligned")

        ("ref",
        value<path>(&ref_sequence_path),
        "File path of FASTA file containing reference sequences to be aligned")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. File will be named based on input file name")

        ("minimap_preset",
        value<string>(&minimap_preset)->
        default_value("map-ont"),
        "Minimap preset configuration name. Affects multiple alignment/mapping parameters")

        ("explicit_mismatch",
        value<bool>(&explicit_mismatch)->
        default_value(true),
        "Whether to use the '='/'X' convention to explicitly label match/mismatch in the SAM cigar string")

        ("k",
        value<uint16_t>(&k)->
        default_value(15),
        "Minimap kmer/seed size")

        ("max_threads",
        value<uint16_t>(&max_threads)->
        default_value(1),
        "Maximum number of threads to launch")

        ("sort",
        value<bool>(&sort)->
        default_value(true),
        "Whether to sort the BAM output")

        ("index",
        value<bool>(&index)->
        default_value(true),
        "Whether to index the BAM output")

        ("keep",
        value<bool>(&keep)->
        default_value(false),
        "Whether to keep intermediate files");

    variables_map vm = parse_arguments(argc, argv, options);
    notify(vm);

    align(ref_sequence_path, read_sequence_path, output_dir, sort, index, !keep, k, minimap_preset, explicit_mismatch, max_threads);

    return 0;
}
