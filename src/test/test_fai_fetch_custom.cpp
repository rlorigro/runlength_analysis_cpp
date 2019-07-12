#include "FastaWriter.hpp"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "boost/algorithm/string/split.hpp"
#include <string>
#include <map>
#include <iostream>
#include <experimental/filesystem>
#include <assert.h>

using std::cout;
using std::cerr;
using std::map;
using std::experimental::filesystem::path;
using boost::split;


bool is_tab(char query){
    return (query == '\t');
}


void read_fasta_index(map<string, uint64_t>& index_map, path fasta_index_path){
    ifstream index_file = ifstream(fasta_index_path);
    vector <string> elements;
    string line;

    // Check if file is readable or exists
    if (!index_file.good()){
        throw runtime_error("ERROR: file read error: " + string(fasta_index_path));
    }

    while(getline(index_file, line)){
        split(elements, line, is_tab);

        index_map[elements[0]] = stoull(elements[2]);

        for (auto& element: elements){
            cout << element << " ";
        }
        cout << "\n";
    }
}


int main(){
    // Find absolute path to test file within repo
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_input_path = "/data/test/";
    path filename = "test_sequences.fasta";
    path fasta_path = project_directory / relative_input_path / filename;

    // Build faidx
    int fai_exit_code = fai_build(fasta_path.c_str());
    if(fai_exit_code != 0) {
        cerr << "Error running faidx_build on " << fasta_path << "\n";
        exit(fai_exit_code);
    }

    // Start reading index file (assuming it was created in above step)
    path fasta_index_path = fasta_path.string() + ".fai";

    cout << "Reading: " << fasta_index_path << "\n";

    map <string,uint64_t> index_map;
    read_fasta_index(index_map, fasta_index_path);

    ifstream fasta_file = ifstream(fasta_path);
    string line;
    for (auto& pair: index_map){
        cout << pair.first << "," << pair.second << "\n";
        fasta_file.seekg(pair.second);
        getline(fasta_file, line);
        cout << line << "\n";
    }
}


/*
    ifstream fasta_file = ifstream(fasta_path);
    ifstream index_file = ifstream(fasta_index_path);

    // Check if file is readable or exists
    if (!fasta_file.good()){
        throw runtime_error("ERROR: file read error: " + string(fasta_index_path));
    }

    // Check if file is readable or exists
    if (!index_file.good()){
        throw runtime_error("ERROR: file read error: " + string(fasta_index_path));
    }

    map <string,uint64_t> index_map;
    vector <string> elements;
    string line;

    // Read index file
    while(getline(index_file, line)){
        split(elements, line, is_tab);

        index_map[elements[0]] = stoull(elements[2]);

        for (auto& element: elements){
            cout << element << " ";
        }
        cout << "\n";
    }

    // Read fasta file (using indexes)
    for (auto& pair: index_map){
        cout << pair.first << "," << pair.second << "\n";
        fasta_file.seekg(pair.second);
        getline(fasta_file, line);
        cout << line << "\n";
    }

*/