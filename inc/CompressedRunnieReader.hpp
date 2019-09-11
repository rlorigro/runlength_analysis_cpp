
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
#include <experimental/filesystem>

using std::pair;
using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;


class CompressedRunnieSequence {
public:
    /// Attributes ///
    string name;
    string sequence;
    vector <uint8_t> encoding;

    /// Methods ///
    void print_encoding();
};


class CompressedRunnieReader{
public:
    /// Attributes ///
    path sequence_file_path;
    int sequence_file_descriptor;

    uint64_t indexes_start_position;
    uint64_t channel_metadata_start_position;
    off_t file_length;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    uint64_t n_channels;

    // What is the unit size of each channel
    vector<uint64_t> channel_sizes;

    vector<CompressedRunnieIndex> indexes;
    unordered_map<string,size_t> index_map;

    /// Methods ///
    CompressedRunnieReader(path file_path);
    void read_footer();
    void read_channel_metadata();
    void read_indexes();
    void read_index_entry(CompressedRunnieIndex& index_element, off_t& byte_index);
    void read_sequence(CompressedRunnieSequence& sequence, CompressedRunnieIndex& index_element);
};

#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
