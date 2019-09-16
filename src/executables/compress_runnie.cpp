#include <iostream>
#include <experimental/filesystem>
#include "RunnieReader.hpp"
#include "CompressedRunnieWriter.hpp"
#include "boost/program_options.hpp"

using std::cout;
using std::flush;
using std::experimental::filesystem::path;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void compress_runnie(path input_dir, path output_dir){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    path relative_config_path = "/data/test/runnie/config/compression_parameters.tsv";     //TODO: move to `/config/`
    path absolute_config_path = project_directory / relative_config_path;

    path output_filename;

    if (input_dir.stem() == ".") {
        output_filename = input_dir.parent_path().stem().string() + ".rq";
    }
    else{
        output_filename = input_dir.stem().string() + ".rq";
    }

    path output_path = output_dir / output_filename;

    cout << "\nWriting to: " << output_path << '\n';

    RunnieReader reader = RunnieReader(input_dir);
    CompressedRunnieWriter writer = CompressedRunnieWriter(output_path, absolute_config_path);

    RunnieSequence sequence;
    string read_name;

    for (auto& index: reader.read_indexes){
        sequence = {};
        read_name = index.first;
        reader.fetch_sequence(sequence, read_name);
        writer.write_sequence(sequence);
    }

    writer.write_indexes();
}


int main(int argc, char* argv[]){
    path input_dir;
    path output_dir;

    options_description options("Required options");

    options.add_options()
        ("input_dir",
        value<path>(&input_dir),
        "File path of directory containing Runnie .out file containing sequences to be compressed")

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

    cout << "READING DIR: " << string(input_dir) << "\n";

    compress_runnie(input_dir, output_dir);

    return 0;
}


