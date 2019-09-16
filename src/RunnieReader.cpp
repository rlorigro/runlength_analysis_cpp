
#include "RunnieReader.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::experimental::filesystem::directory_iterator;
using std::experimental::filesystem::is_regular_file;
using std::replace;
using std::ostream;
using std::exception;
using std::cout;
using std::cerr;


RunnieIndex::RunnieIndex(path file_path, uint64_t byte_index, uint64_t length){
    this->file_path = file_path;
    this->byte_index = byte_index;
    this->length = length;
}


ostream& RunnieIndex::operator<<(ostream& s){
    cout << this->file_path << " " << this->byte_index << " " << this->length;

    return s;
}


RunnieReader::RunnieReader(path directory_path){
    this->directory_path = directory_path;
    this->index();
}


void RunnieReader::update_index(path& file_path,
        ifstream& file,
        string& line,
        uint64_t& l,
        string& read_name,
        uint64_t& byte_index,
        uint64_t& read_length,
        uint64_t& current_header_line_index,
        uint64_t& current_header_byte_index){

    // If a full read has been parsed, i.e. a header other than the first header has been reached, do some stuff
    if (l != 0) {
        read_length = l - current_header_line_index - 1;
        byte_index = current_header_byte_index;

        bool no_conflict = this->read_indexes.try_emplace(read_name, RunnieIndex(file_path, byte_index, read_length)).second;
        if (not no_conflict) {
            throw runtime_error(
                    "ERROR: duplicate or nonexistent read detected in Runnie: " + file_path.string() + " " + read_name);
        }
    }

    // If it's not the end of the file
    if (not file.eof()) {
        // Current header saved until next header is found
        read_name = line.substr(2, line.size() - 1);
        current_header_line_index = l;
        current_header_byte_index = file.tellg();
    }
}


void RunnieReader::index_file(path file_path){
    // Open file
    ifstream file = ifstream(file_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not open file " + file_path.string());
    }

    // File iteration variables
    string line;
    uint64_t l = 0;

    // Read iteration variables
    string read_name;
    uint64_t byte_index = 0;
    uint64_t read_length = 0;
    uint64_t current_header_line_index = 0;
    uint64_t current_header_byte_index = 0;

    while (getline(file, line)){
        if (line[0] == '#'){
            this->update_index(file_path,file,line,l,read_name,byte_index,read_length,current_header_line_index,current_header_byte_index);
        }
        l++;
    }
    this->update_index(file_path,file,line,l,read_name,byte_index,read_length,current_header_line_index,current_header_byte_index);
}


unordered_map <string, RunnieIndex> RunnieReader::get_index(){
    return this->read_indexes;
}


void RunnieReader::set_index(unordered_map <string, RunnieIndex> index){
    this->read_indexes = index;
}



void RunnieReader::index(){
    ///
    /// Load all the filenames and iterate their contents to find read names and byte indexes
    ///
    path filename_prefix;
    string read_name;

    for (const path& file_path: directory_iterator(this->directory_path)){
        if (is_regular_file(file_path)) {
            this->index_file(file_path);
        }
    }
}


void RunnieReader::parse_line(RunnieSequence& sequence, string& line){
    sequence.sequence += line[0];

    string scale_string;
    string shape_string;
    double scale;
    double shape;

    uint8_t column_index = 0;
    for (auto& c: line){
        if (c == '\t'){
            // Count the current column
            column_index++;
        }
        else{
            // Shape param (beta) is in column 1
            if (column_index == 1){
                shape_string += c;
            }
            // Scale param (alpha) is in column 2
            else if (column_index == 2){
                scale_string += c;
            }
        }
    }

    // Convert string representations of hex floating points to doubles
    scale = stod(scale_string);
    shape = stod(shape_string);

    sequence.scales.emplace_back(scale);
    sequence.shapes.emplace_back(shape);
}


void RunnieReader::fetch_sequence_bases(RunnieSequence& sequence, string& read_name){
    ///
    /// Dont read the whole file, just fetch the nucleotide bases for the sequence
    ///

    // Clear the container
    sequence = {};

    // Name is known from prior indexing
    sequence.name = read_name;

    // Fetch index
    RunnieIndex read_index = this->read_indexes.at(read_name);

    if (this->read_indexes.empty()){
        cerr << "Index not loaded for runnie directory, generating now...\n";
        this->index();
    }

    // Open file
    ifstream file = ifstream(read_index.file_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not open file " + read_index.file_path.string());
    }

    // Skip to sequence start position in file
    file.seekg(read_index.byte_index);

    // File iteration variables
    string line;
    uint64_t l = 0;

    // The sequence length is known from the first pass during indexing
    while (getline(file, line) and (l < read_index.length)){
        sequence.sequence += line[0];
        l++;
    }
}



void RunnieReader::fetch_sequence(RunnieSequence& sequence, string& read_name){
    // Clear the container
    sequence = {};

    // Name is known from prior indexing
    sequence.name = read_name;

    if (this->read_indexes.empty()){
        cerr << "Index not loaded for runnie directory, generating now...\n";
        this->index();
    }

    // Fetch index
    RunnieIndex read_index = this->read_indexes.at(read_name);

    // Open file
    ifstream file = ifstream(read_index.file_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not open file " + read_index.file_path.string());
    }

    // Skip to sequence start position in file
    file.seekg(read_index.byte_index);

    // File iteration variables
    string line;
    uint64_t l = 0;

    // The sequence length is known from the first pass during indexing
    while (getline(file, line) and (l < read_index.length)){
        try {
            this->parse_line(sequence, line);
        }
        catch (const exception& e){
            cout << e.what() << "\n";
            cout << "runnie parser failed at line " << l << " in file: " << read_index.file_path << "\n";
        }
        l++;
    }
}


void RunnieReader::fetch_all_sequences(vector<RunnieSequence>& sequences){
    RunnieSequence sequence;
    string read_name;

    for (auto& item: this->read_indexes){
        read_name = item.first;
        sequences.push_back(sequence);
        this->fetch_sequence(sequences.back(), read_name);
    }
}

