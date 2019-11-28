
#include "FastaReader.hpp"
#include "Miscellaneous.hpp"
#include <algorithm>
#include <experimental/filesystem>
#include <boost/program_options.hpp>
#include "boost/algorithm/string.hpp"

using std::make_heap;
using std::push_heap;
using std::pop_heap;
using std::sort_heap;
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::trim_left_if;
using boost::trim_right;
using boost::split;


bool is_tab(char c){
    return (c == '\t');
}


void read_lengths_from_fasta_index(path index_path, vector<uint32_t>& lengths){
    ifstream index_file = ifstream(index_path);
    vector <string> elements;
    string line;

    // Check if file is readable or exists
    if (!index_file.good()){
        throw runtime_error("ERROR: file read error: " + string(index_path));
    }

    // Iterate .fai to collect byte start positions of sequences for each fasta element
    while(getline(index_file, line)){
        split(elements, line, is_tab);

        // Each index element is a pair of sequence name (column 0), byte position (column 2), and sequence length (column 1)
        lengths.emplace_back(stoul(elements[1]));
    }
}


void measure_read_length_stats_from_fasta(string comma_separated_paths, path output_dir, uint16_t max_threads) {
    vector<string> paths;
    vector<uint32_t> lengths;
    string separators = ",";
    split_as_string(paths, comma_separated_paths, separators);

    create_directories(output_dir);
    path output_path = output_dir / "read_lengths.txt";
    ofstream out_file(output_path);

    for (auto& path: paths) {
        lengths = {};

        out_file << '>' << path << '\n';

        FastaReader reader(path);
        reader.build_fasta_index();

        read_lengths_from_fasta_index(reader.index_path, lengths);
        make_heap(lengths.begin(), lengths.end());
        sort_heap(lengths.begin(), lengths.end());

        for (size_t i=0; i<lengths.size(); i++){
            out_file << lengths[i];
            if (i < lengths.size()-1){
                out_file << ',';
            }
        }
        out_file << '\n';
    }
}


int main(int argc, char* argv[]){
    string fasta_paths;
    path output_dir;
    uint16_t max_threads;

    options_description options("Arguments");

    options.add_options()
            ("fastas",
             value<string>(&fasta_paths),
             "File path of reference FASTA file containing QUERY sequences to be Run-length encoded")

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

    measure_read_length_stats_from_fasta(fasta_paths,
            output_dir,
            max_threads);

    return 0;
}

