
#include "ReferenceRunlength.hpp"
#include "boost/program_options.hpp"
#include <iostream>
#include <experimental/filesystem>

using std::cout;
using std::pair;
using std::experimental::filesystem::path;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


int main(int argc, char* argv[]){
    path input_path;
    uint16_t max_length;

    options_description options("Arguments");

    options.add_options()
        ("fasta",
        value<path>(&input_path),
        "Path to FASTA file containing (reference) sequences")

        ("max_length",
        value<uint16_t>(&max_length)->
        default_value(50),
        "Maximum runlength to count when collecting length observations");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    measure_runlength_priors_from_reference(input_path, max_length);

    return 0;
}
