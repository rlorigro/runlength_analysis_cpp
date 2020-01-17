#include "Runlength.hpp"
#include "boost/program_options.hpp"
#include <iostream>

using std::cout;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path bed_path;
    path input_dir;
    path output_dir;
    uint16_t max_threads;

    options_description options("Arguments");

    options.add_options()
        ("ref",
        value<path>(&ref_fasta_path),
        "File path of reference FASTA file containing sequences to be Run-length encoded")

        ("input_dir",
        value<path>(&input_dir),
        "Path to directory containing Shasta CSVs")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. File will be named based on input file name")

        ("bed",
        value<path>(&bed_path)->
        default_value(path()),
        "File path of BED file contain regions to subset ")

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

    label_coverage_data_from_shasta(input_dir,
            ref_fasta_path,
            output_dir,
            max_threads,
            bed_path);

    return 0;
}

