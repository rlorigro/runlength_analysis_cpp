
#include "FastaReader.hpp"
#include "Miscellaneous.hpp"
#include <algorithm>
#include <random>
#include <experimental/filesystem>
#include <boost/program_options.hpp>
#include "boost/algorithm/string.hpp"

using std::cerr;
using std::make_heap;
using std::push_heap;
using std::pop_heap;
using std::sort_heap;
using std::random_device;
using std::mt19937;
using std::shuffle;
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


void read_lengths_from_fasta_index(path index_path, vector<uint32_t>& lengths, uint32_t min_length=0){
    ifstream index_file = ifstream(index_path);
    vector <string> elements;
    string line;
    uint32_t length;

    // Check if file is readable or exists
    if (!index_file.good()){
        throw runtime_error("ERROR: file read error: " + string(index_path));
    }

    // Iterate .fai to collect byte start positions of sequences for each fasta element
    while(getline(index_file, line)){
        split(elements, line, is_tab);

        // Each index element is a pair of sequence name (column 0), byte position (column 2), and sequence length (column 1)
        length = stoul(elements[1]);

        if (length > min_length){
            lengths.emplace_back(length);
        }
    }
}


void measure_read_length_stats_from_fasta(string comma_separated_paths, path output_dir, uint32_t min_length, uint64_t max_cumulative_length, uint16_t max_threads) {
    vector<string> paths;
    vector<uint32_t> lengths;
    vector<uint32_t> lengths_subset;
    string separators = ",";
    split_as_string(paths, comma_separated_paths, separators);
    uint64_t cumulative_length;

    create_directories(output_dir);
    path output_path = output_dir / "read_lengths.txt";
    output_path = absolute(output_path);
    ofstream out_file(output_path);

    cerr << "WRITING: " << output_path << '\n';

    for (auto& path: paths) {
        cumulative_length = 0;
        lengths_subset = {};
        lengths = {};

        out_file << '>' << path << '\n';

        FastaReader reader(path);
        reader.build_fasta_index();

        read_lengths_from_fasta_index(reader.index_path, lengths);

        random_device device;
        mt19937 random_generator(device());
        shuffle(lengths.begin(), lengths.end(), random_generator);

        for (auto &length: lengths) {
            if (cumulative_length < max_cumulative_length) {
                lengths_subset.emplace_back(length);
                cumulative_length += length;
            } else {
                break;
            }
        }

        make_heap(lengths_subset.begin(), lengths_subset.end());
        sort_heap(lengths_subset.begin(), lengths_subset.end());

        for (size_t i=0; i<lengths_subset.size(); i++) {
            out_file << lengths_subset[i];
            if (i < lengths_subset.size() - 1) {
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
    uint32_t min_length;
    uint64_t max_cumulative_length;

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
             "Maximum number of threads to launch")

            ("min_length",
             value<uint32_t>(&min_length)->
                     default_value(1),
             "Minimum length of read to record")

             ("max_cumulative_length",
             value<uint64_t>(&max_cumulative_length)->
                     default_value(std::numeric_limits<uint64_t>::max()),
             "Minimum length of read to record");

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
            min_length,
            max_cumulative_length,
            max_threads);

    return 0;
}

