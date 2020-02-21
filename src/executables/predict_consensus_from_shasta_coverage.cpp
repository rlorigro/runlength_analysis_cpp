#include "PileupGenerator.hpp"
#include "RunlengthReader.hpp"
#include "RunlengthWriter.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Identity.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include <iostream>
#include <SimpleBayesianConsensusCaller.hpp>
#include <thread>
#include <mutex>
#include <atomic>
#include <exception>


// TODO FINISH THIS STUB

using std::exception;
using std::thread;
using std::mutex;
using std::lock_guard;
using std::move;
using std::ref;
using std::atomic;
using std::atomic_fetch_add;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::min;
using std::max;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using std::experimental::filesystem::create_directories;


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path input_dir;
    path output_dir;

    options_description options("Arguments");

    options.add_options()
        ("ref",
        value<path>(&ref_fasta_path),
        "File path of FASTA file containing REFERENCE sequences to be Run-length encoded")

        ("input_dir",
        value<path>(&input_dir),
        "File path of FASTA file containing QUERY sequences to be Run-length encoded")

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

    return 0;
}
