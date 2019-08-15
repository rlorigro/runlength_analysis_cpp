
#include "RunnieReader.hpp"
#include <algorithm>
#include <iostream>

using std::experimental::filesystem::directory_iterator;
using std::replace;
using std::ostream;
using std::cout;


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
                    "ERROR: duplicate reads detected in Runnie: " + file_path.string() + " " + read_name);
        }
        cout << read_name << " " << file_path.string() << " " << byte_index << " " << read_length << "\n";
    }

    // If its not the end of the file, update as if it
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
//        cout << line << "\n";
        if (line[0] == '#'){
            this->update_index(file_path,file,line,l,read_name,byte_index,read_length,current_header_line_index,current_header_byte_index);
        }
        l++;
    }
    this->update_index(file_path,file,line,l,read_name,byte_index,read_length,current_header_line_index,current_header_byte_index);
}


void RunnieReader::index(){
    ///
    /// Load all the filenames and iterate their contents to find read names and byte indexes
    ///
    path filename_prefix;
    string read_name;

    for (const path& file_path: directory_iterator(this->directory_path)){
        this->index_file(file_path);
    }

}