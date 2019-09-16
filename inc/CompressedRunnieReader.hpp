
#ifndef RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
#define RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP

#include "CompressedRunnieWriter.hpp"
#include "BinaryIO.hpp"
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <unordered_map>

using std::pair;
using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using std::unordered_map;


//ostream& operator<<(ostream& s, CompressedRunnieIndex& index);


class CompressedRunnieSequence {
public:
    /// Attributes ///
    string sequence;
    vector <uint8_t> encoding;

    /// Methods ///
    void print_encoding();
};


class NamedCompressedRunnieSequence: public CompressedRunnieSequence{
public:
    /// Attributes ///
    string name;
};


class CompressedRunnieReader{
public:

    /// Methods ///

    // Initialize the class with a file path
    CompressedRunnieReader(string file_path);

    // Fetch the name of a read based on its number (ordering in file, 0-based)
    const string& get_read_name(uint64_t read_number);

    // Fetch the length of a read based on its number (ordering in file, 0-based)
    uint64_t get_length(uint64_t read_number);

    // Fetch the sequence of a read based on its number (ordering in file, 0-based)
    void get_sequence_data(CompressedRunnieSequence& sequence, uint64_t read_number);

    // Fetch sequence data, and the 'name' field is also filled in.
    void get_sequence_data(NamedCompressedRunnieSequence& sequence, uint64_t read_number);

    // Fetch the number of reads in the file
    size_t get_read_count();

    // Return the file path that this reader is reading from
    const string& get_file_name();

private:

    /// Attributes ///
    string sequence_file_path;
    int sequence_file_descriptor;

    uint64_t indexes_start_position;
    uint64_t channel_metadata_start_position;
    off_t file_length;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    uint64_t n_channels;

    // What is the unit size of each channel
    vector<uint64_t> channel_sizes;

    /// Methods ///
    void read_footer();
    void read_channel_metadata();
    void read_indexes();
    void read_index_entry(CompressedRunnieIndex& index_element, off_t& byte_index);

    vector<CompressedRunnieIndex> indexes;
    unordered_map<string,size_t> index_map;

};

#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
