#include "MarginPolishReader.hpp"
#include "boost/program_options.hpp"
#include <iostream>
#include <utility>

using std::cout;
using std::pair;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void count_marginpolish_coverage(path input_directory){
    MarginPolishReader marginpolish_reader = MarginPolishReader(input_directory);
    CoverageSegment segment;
    marginpolish_reader.index();

    double total_coverage = 0;

    // Iterate file names
    for (pair<string,path> element: marginpolish_reader.file_paths){
        marginpolish_reader.fetch_read(segment, element.first);

        for (auto& pileup: segment.coverage_data){
            for (auto& element: pileup){
                total_coverage += element.weight;
            }
        }
    }

    cout << "total_coverage: " << int(total_coverage) << '\n';
}


int main(int argc, char* argv[]){
    path input_dir;

    options_description options("Arguments");

    options.add_options()
        ("input_dir",
        value<path>(&input_dir),
        "Path to directory containing MarginPolish TSVs");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    count_marginpolish_coverage(input_dir);

    return 0;
}
