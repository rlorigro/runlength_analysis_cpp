
#ifndef RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
#define RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP

#include "CompressedRunnieWriter.hpp"
#include "RunnieReader.hpp"
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


class CompressedRunnieReader{
public:
    /// Attributes ///
    path sequence_file_path;
    path index_file_path;
    path params_path;
    ifstream sequence_file;

    uint64_t indexes_start_position;
    uint64_t channel_metadata_start_position;
    uint64_t file_length;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    uint64_t n_channels;

    // What is the unit size of each channel
    uint64_t channel_size_1;

    vector<CompressedRunnieIndex> indexes;
    unordered_map<string,size_t> index_map;

    /// Methods ///
    CompressedRunnieReader(path file_path);
    void read_footer();
    void read_indexes();
    void read_index_entry(CompressedRunnieIndex& index_element);

};

#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
