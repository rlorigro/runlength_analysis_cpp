#include <iostream>
#include <experimental/filesystem>
#include <assert.h>
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Runlength.hpp"
#include "boost/program_options.hpp"

using std::cout;
using std::flush;
using std::experimental::filesystem::path;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void fasta_to_RLE_fasta(path input_file_path, path output_dir) {
    // Generate parent directories if necessary
    create_directories(output_dir);

    string output_filename;
    output_filename = string(input_file_path.filename());
    output_filename = output_filename.substr(0, output_filename.find_last_of(".")) + "_RLE.fasta";

    path output_file_path = output_dir / output_filename;

    FastaReader fasta_reader(input_file_path);
    FastaWriter fasta_writer(output_file_path);

    string line;
    SequenceElement element;
    RunlengthSequenceElement runlength_element;

    while (!fasta_reader.end_of_file) {
        // Initialize empty containers
        element = {};
        runlength_element = {};

        // Parse next sequence element
        fasta_reader.next_element(element);

        // Print status update to stdout
        cout << "\33[2K\r" << element.name << flush;

        // Convert to Run-length Encoded sequence element
        runlength_encode(runlength_element, element);

        // Write RLE sequence to file (no lengths written)
        fasta_writer.write(runlength_element);
    }
}


int main(int argc, char* argv[]){
    path input_file_path;
    path output_dir;

    options_description options("Required options");

    options.add_options()
        ("fasta",
        value<path>(&input_file_path),
        "File path of FASTA file containing sequences to be Run-length encoded")

        ("output_dir",
        value<path>(&output_dir)->
        default_value("output/"),
        "Destination directory. File will be named based on input file name");

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

    fasta_to_RLE_fasta(input_file_path, output_dir);
}

