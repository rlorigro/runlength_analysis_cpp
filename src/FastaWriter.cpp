#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <experimental/filesystem>
#include "FastaWriter.hpp"
#include "FastaReader.hpp"

using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;


FastaWriter::FastaWriter(path file_path) {
    this->file_path = file_path;
    this->fasta_file = ofstream(file_path);
    this->header_symbol = '>';
    this->index = 0;

    // Check if file is readable or exists
    if (!this->fasta_file.good()){
        throw runtime_error("ERROR: file write error: " + string(file_path));
    }
}


