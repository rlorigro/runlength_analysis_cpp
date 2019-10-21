
#ifndef RUNLENGTH_ANALYSIS_RUNLENGTHWRITER_HPP
#define RUNLENGTH_ANALYSIS_RUNLENGTHWRITER_HPP

#include "RunlengthSequenceElement.hpp"
#include "RunlengthIndex.hpp"
#include "Miscellaneous.hpp"
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


ostream& operator<<(ostream& s, RunlengthIndex& index);


class RunlengthWriter {
public:
    /// Attributes ///
    path sequence_file_path;
    ofstream sequence_file;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    static const uint64_t n_channels = 1;

    // What is the unit size of that channel
    static const vector<uint64_t> channel_sizes;

    // Named constant for accessing channels
    static const size_t LENGTH = 0;

    // When writing the binary file, this vector is appended, so the position of each sequence is stored
    vector<RunlengthIndex> indexes;

    /// Methods ///
    RunlengthWriter(path file_path);

    void write_sequence(RunlengthSequenceElement& sequence);
    void write_sequence_block(RunlengthSequenceElement& sequence);
    void write_length_block(RunlengthSequenceElement& sequence);
    void write_index(RunlengthIndex& index);
    void write_indexes();
};


#endif //RUNLENGTH_ANALYSIS_RUNLENGTHWRITER_HPP
