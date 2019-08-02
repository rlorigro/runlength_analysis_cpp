#include "Runlength.hpp"
#include "FastaReader.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <experimental/filesystem>
#include <assert.h>

using std::cout;
using std::map;
using std::vector;
using std::experimental::filesystem::path;

int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/test_sequences.fasta";
    path absolute_data_path = project_directory / relative_data_path;

    cout << absolute_data_path << "\n";

    FastaReader fasta_reader(absolute_data_path);

    vector<string> sequence_names = {"test1", "test2", "test3", "test4"};
    map<string,string> sequence_truth_map;
    map<string,vector<uint16_t>> length_truth_map;

    for (auto& name: sequence_names){
        if (name == "test1"){
            sequence_truth_map[name] = "ACAC";
            vector<uint16_t> lengths = {1,2,3,4};
            length_truth_map[name] = lengths;
        }
        if (name == "test2"){
            sequence_truth_map[name] = "GTGT";
            length_truth_map[name] = vector<uint16_t> {4,3,2,1};
        }
        if (name == "test3"){
            sequence_truth_map[name] = "GTGTGTGTGTGT";
            length_truth_map[name] = vector<uint16_t> {4,3,2,1,4,3,2,1,4,3,2,1};
        }
        if (name == "test4"){
            sequence_truth_map[name] = "ACAC";
            length_truth_map[name] = vector<uint16_t> {1,2,3,4};
        }
    }

    string line;
    SequenceElement element;
    runlength_sequence_element runlength_element;

    while (!fasta_reader.end_of_file) {
        // Initialize empty containers
        element = {};
        runlength_element = {};

        // Parse next sequence element
        fasta_reader.next_element(element);

        // Convert to Run-length Encoded sequence element
        runlength_encode(runlength_element, element);

        assert (runlength_element.sequence == sequence_truth_map[runlength_element.name]);
        cout << "PASS\n";

        assert (runlength_element.lengths == length_truth_map[runlength_element.name]);
        cout << "PASS\n";

        cout << "\n";
    }
}
