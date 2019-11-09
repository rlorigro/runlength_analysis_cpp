
#ifndef RUNLENGTH_ANALYSIS_RUNLENGTHREADER_HPP
#define RUNLENGTH_ANALYSIS_RUNLENGTHREADER_HPP

#include "RunlengthSequenceElement.hpp"
#include "RunlengthIndex.hpp"
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


class RunlengthReader{
public:

    /// Methods ///

    // Initialize the class with a file path
    RunlengthReader(string file_path);

    // Fetch the name of a read based on its number (ordering in file, 0-based)
    const string& get_read_name(uint64_t read_number);

    // Fetch the length of a read based on its number (ordering in file, 0-based)
    uint64_t get_length(uint64_t read_number);

    // Fetch the sequence of a read based on its number (ordering in file, 0-based)
    void get_sequence(RunlengthSequenceElement& sequence, uint64_t read_number);

    // Fetch the sequence of a read based on its number (ordering in file, 0-based)
    void get_sequence(RunlengthSequenceElement& sequence, string& read_name);

    // Fetch the number of reads in the file
    size_t get_read_count();

    // Return the file path that this reader is reading from
    const string& get_file_name();

    // Sometimes you just need to take a shortcut
    RunlengthSequenceElement generate_sequence_container();

    vector<RunlengthIndex> indexes;

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
    void read_index_entry(RunlengthIndex& index_element, off_t& byte_index);

    unordered_map<string,size_t> index_map;
};

#endif //RUNLENGTH_ANALYSIS_RUNLENGTHREADER_HPP
