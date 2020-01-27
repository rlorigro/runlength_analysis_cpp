#include <iostream>
#include <experimental/filesystem>
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Runlength.hpp"
#include "boost/program_options.hpp"
#include "Runlength.hpp"

using std::cout;
using std::flush;
using std::experimental::filesystem::path;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void fasta_to_RLE_fasta(path input_file_path, path output_dir, uint max_threads) {
    // Generate parent directories if necessary
    create_directories(output_dir);

    // Runlength encode the READ SEQUENCES, rewrite to another FASTA, and DON'T store in memory
    unordered_map<string, RunlengthSequenceElement> _;
    path reads_fasta_path_rle;
    bool store_in_memory = false;
    reads_fasta_path_rle = runlength_encode_fasta_file(
            input_file_path,
            _,
            output_dir,
            store_in_memory,
            max_threads);

    cout << '\n';
}


int main(int argc, char* argv[]){
    path input_file_path;
    path output_dir;
    uint16_t max_threads;

    options_description options("Required options");

    options.add_options()
        ("fasta",
        value<path>(&input_file_path),
        "File path of FASTA file containing sequences to be Run-length encoded")

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

    cout << "READING FILE: " << string(input_file_path) << "\n";

    fasta_to_RLE_fasta(input_file_path, output_dir, max_threads);

    return 0;
}

