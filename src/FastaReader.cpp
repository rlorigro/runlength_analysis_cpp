#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <experimental/filesystem>
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "boost/algorithm/string.hpp"
#include "FastaReader.hpp"

using std::map;
using std::string;
using std::vector;
using std::to_string;
using std::cout;
using std::cerr;
using std::ifstream;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::exists;
using boost::trim_left_if;
using boost::trim_right;
using boost::split;


bool is_caret(char c){
    return (c == '>');
}


bool is_tab(char c){
    return (c == '\t');
}


FastaReader::FastaReader(path file_path){
    this->file_path = file_path;
    this->index_path = file_path.string() + ".fai";
    this->fasta_file = ifstream(file_path);
    this->header_symbol = '>';
    this->end_of_file = false;
    this->line_index = 0;

    // Check if file is readable or exists
    if (!this->fasta_file.good()){
        throw runtime_error("ERROR: file read error: " + string(file_path));
    }
}


map<string,uint64_t> FastaReader::get_index(){
    if (this->read_indexes.empty()){
        this->index();
    }

    return this->read_indexes;
}


void FastaReader::set_index(map <string,uint64_t>& index){
    this->read_indexes = index;
}


void FastaReader::fetch_sequence(SequenceElement& element, string& sequence_name){
    if (this->read_indexes.empty()){
        this->index();
    }

    // Fill in the "name" field of the sequence element
    element.name = sequence_name;

    if (this->fasta_file.tellg() == -1){
        this->end_of_file = false;
        this->fasta_file.clear();
    }

    // Set ifstream cursor to the start of this read's sequence
    this->fasta_file.seekg(this->read_indexes[sequence_name]);

    // Fill in the "sequence" field of the sequence element
    this->read_next_sequence(element);
}


void FastaReader::fetch_sequence(SequenceElement& element, string& sequence_name, uint64_t fasta_byte_index){
    // Fill in the "name" field of the sequence element
    element.name = sequence_name;

    if (this->fasta_file.tellg() == -1){
        this->end_of_file = false;
        this->fasta_file.clear();
    }

    // Set ifstream cursor to the start of this read's sequence
    this->fasta_file.seekg(fasta_byte_index);

    // Fill in the "sequence" field of the sequence element
    this->read_next_sequence(element);
}


void FastaReader::index(){
    // If not fai exists build one
    this->build_fasta_index();

    // Read fai and store in this->read_indexes
    this->read_fasta_index();
}


void FastaReader::build_fasta_index(){
    if (!exists(this->index_path)){
        cout << "No index found, generating .fai for " << this->file_path << "\n";

        // Build faidx using htslib
        int fai_exit_code = fai_build(this->file_path.c_str());
        if(fai_exit_code != 0) {
            throw runtime_error("Error running faidx_build on " + this->file_path.string() + "\n");
        }
    }
}


void FastaReader::read_fasta_index(){
    ifstream index_file = ifstream(this->index_path);
    vector <string> elements;
    string line;

    // Check if file is readable or exists
    if (!index_file.good()){
        throw runtime_error("ERROR: file read error: " + string(this->index_path));
    }

    // Iterate .fai to collect byte start positions of sequences for each fasta element
    while(getline(index_file, line)){
        split(elements, line, is_tab);

        // Each index element is a pair of sequence name (column 0) and byte position (column 2)
        read_indexes[elements[0]] = stoull(elements[2]);
    }
}


// Fetch header which starts at current ifstream position
void FastaReader::read_next_header(SequenceElement& element){
    // Try fetching header line
    getline(this->fasta_file, element.name);
    this->line_index++;

    // Update class status
    if (this->fasta_file.peek() == EOF){
        this->end_of_file = true;
    }

    // Verify that header starts with header symbol '>', and trim the symbol from the line. Also trim trailing space.
    if (element.name[0] == this->header_symbol){
        trim_left_if(element.name, &is_caret);
        trim_right(element.name);
    }
    else{
        throw runtime_error("Unrecognized FASTA header character on line: " + to_string(this->line_index));
    }
}


// Fetch sequence which starts at current ifstream position
void FastaReader::read_next_sequence(SequenceElement& element){
    // Make sure the element is empty before starting
    if (element.sequence.size() != 0){
        element = {};
    }

    // Try fetching sequence line (of which there may be multiple)
    while (this->fasta_file.peek() != std::char_traits<char>::to_int_type(this->header_symbol) && !this->end_of_file) {
        // If this is a single line, just place it in the sequence object
        if (element.sequence.length() == 0){
            getline(this->fasta_file, element.sequence);
        }
        // If this is a multi-line FASTA, append the sequence object using a copy of the line (sadly)
        else{
            string buffer;
            getline(this->fasta_file, buffer);
            element.sequence += buffer;
        }
        this->line_index++;

        // Trim any trailing whitespace on the sequence line
        trim_right(element.sequence);

        // Update class status
        if (this->fasta_file.peek() == EOF){
            this->end_of_file = true;
        }
    }
}


// Read the next lines of the file and update the SequenceElement object. If no lines exists, toggle this->file_end
void FastaReader::next_element(SequenceElement& element){
    this->read_next_header(element);
    this->read_next_sequence(element);
}
