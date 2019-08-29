#include "Runlength.hpp"
#include "boost/program_options.hpp"
#include <iostream>

using std::cout;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path reads_fasta_path;
    path output_dir;
    uint16_t max_threads;
    uint16_t max_runlength;

    options_description options("Arguments");

    options.add_options()
        ("ref",
        value<path>(&ref_fasta_path),
        "File path of reference FASTA file containing REFERENCE sequences to be Run-length encoded")

        ("sequences",
        value<path>(&reads_fasta_path),
        "File path of reference FASTA file containing QUERY sequences to be Run-length encoded")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. File will be named based on input file name")

        ("max_threads",
        value<uint16_t>(&max_threads)->
        default_value(1),
        "Maximum number of threads to launch")

        ("max_runlength",
        value<uint16_t>(&max_runlength)->
        default_value(50),
        "Maximum length of a run to use in the model");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    measure_runlength_distribution_from_fasta(reads_fasta_path,
        ref_fasta_path,
        output_dir,
        max_runlength,
        max_threads);

    return 0;
}

