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

    unordered_map <string, RunlengthSequenceElement> runlength_sequences;

//    runlength_encode_fasta_file(absolute_data_path, runlength_sequences, "output/", 4);
    runlength_encode_fasta_file("/home/ryan/data/GRCh38/hg38_no_alts.fasta", runlength_sequences, "output/", false, 32);

    cout << runlength_sequences.size() << "\n";
    for (auto& element: runlength_sequences){
        cout << element.first << "\n";
    }
}
