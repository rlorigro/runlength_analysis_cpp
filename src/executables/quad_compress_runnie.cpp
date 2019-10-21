#include <iostream>
#include <experimental/filesystem>
#include "RunnieReader.hpp"
#include "QuadLoss.hpp"
#include "QuadCompressor.hpp"
#include "boost/program_options.hpp"

using std::cout;
using std::cerr;
using std::flush;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void compress_runnie(path input_dir, path output_dir){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path output_filename;
    input_dir = absolute(input_dir);

    RunnieReader reader = RunnieReader(input_dir);
    reader.index();

    float x = 25;
    float y = 25;

    float size = 50;

    QuadCoordinate center = QuadCoordinate(x, y);
    BoundingBox bounds = BoundingBox(center, size/2);
    QuadCompressor tree = QuadCompressor(bounds);

    RunnieSequence sequence;
    string read_name;
    double scale;
    double shape;

    size_t cutoff = 10*1000*1000;
    size_t n_bases = 0;
    for (auto& index: reader.read_indexes){
        read_name = index.first;
        reader.fetch_sequence(sequence, read_name);

        for (size_t i=0; i<sequence.scales.size(); i++){
            scale = sequence.scales[i];
            shape = sequence.shapes[i];

            tree.insert(QuadCoordinate(scale, shape));

            n_bases++;
        }

        if (n_bases > cutoff){
            break;
        }
    }

    DiscreteWeibullLoss loss = DiscreteWeibullLoss(50);

    vector<QuadCoordinate> results;
    tree.query_range(results, bounds);

    tree.compress(255, loss, output_dir);
}


int main(int argc, char* argv[]){
    path input_dir;
    path output_dir;

    options_description options("Required options");

    options.add_options()
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

    compress_runnie(input_dir, output_dir);

    return 0;
}
