#include <iostream>
#include <experimental/filesystem>
#include "FastqReader.hpp"
#include "boost/program_options.hpp"

using std::cout;
using std::flush;
using std::experimental::filesystem::path;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void filter_fastq_by_quality(path fastq_path, uint64_t minimum_average_quality){
    FastqReader reader(fastq_path);
    reader.filter_by_quality(minimum_average_quality);
}


int main(int argc, char* argv[]){
    path input_file_path;
    path output_dir;
    uint16_t max_threads;

    options_description options("Required options");

    options.add_options()
        ("fastq",
        value<path>(&input_file_path),
        "File path of FASTA file containing sequences to be Run-length encoded")

        ("q",
        value<uint16_t>(&max_threads)->
        default_value(0),
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

    cout << "READING FILE: " << string(input_file_path) << "\n";

    filter_fastq_by_quality(input_file_path, max_threads);

    return 0;
}

