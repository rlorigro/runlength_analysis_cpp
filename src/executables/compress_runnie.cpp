#include <iostream>
#include <experimental/filesystem>
#include "RunnieReader.hpp"
#include "CompressedRunnieWriter.hpp"
#include "boost/program_options.hpp"

using std::cout;
using std::cerr;
using std::flush;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void compress_runnie(path config_path, path input_dir, path output_dir){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    path output_filename;

    config_path = absolute(config_path);
    input_dir = absolute(input_dir);

    if (input_dir.stem() == ".") {
        output_filename = input_dir.parent_path().stem().string() + ".rq";
    }
    else{
        output_filename = input_dir.stem().string() + ".rq";
    }

    path output_path = output_dir / output_filename;

    cerr << "\nWriting to: " << output_path << '\n';

    RunnieReader reader = RunnieReader(input_dir);
    reader.index();

    CompressedRunnieWriter writer = CompressedRunnieWriter(output_path, config_path);

    RunnieSequence sequence;
    string read_name;

    for (auto& index: reader.read_indexes){
        read_name = index.first;
        reader.fetch_sequence(sequence, read_name);
        writer.write_sequence(sequence);
    }

    writer.write_indexes();
}


int main(int argc, char* argv[]){
    path config_path;
    path input_dir;
    path output_dir;

    options_description options("Required options");

    options.add_options()
        ("config_path",
        value<path>(&config_path)->required(),
        "File path of TSV containing bounds and centroids for discretization")

        ("input_dir",
        value<path>(&input_dir)->required(),
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

    compress_runnie(config_path, input_dir, output_dir);

    return 0;
}


