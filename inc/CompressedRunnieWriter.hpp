
#ifndef RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEWRITER_HPP
#define RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEWRITER_HPP

#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"
#include "RunnieReader.hpp"
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <experimental/filesystem>

using boost::icl::interval_map;
using boost::icl::interval;
using boost::icl::total_enricher;
using std::pair;
using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;

using lower_interval_map = interval_map<double,uint8_t,total_enricher>;


class CompressedRunnieIndex {
public:
    /// Attributes ///
    string name;
    uint64_t sequence_byte_index;
    uint64_t sequence_length;

    /// Methods ///
};

class CompressedRunnieWriter {
public:
    /// Attributes ///
    path sequence_file_path;
    path index_file_path;
    path params_path;
    ofstream sequence_file;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    const uint64_t n_channels = 1;

    // What is the unit size of that channel
    const uint64_t channel_size_1 = sizeof(uint8_t);

    // Interval data structures used to initialize the interval map
    vector <pair <double,double> > scale_intervals;
    vector < vector <pair <double,double> > > shape_intervals;

    // The interval map encodes any given pair of scale/shape.
    // The first tree is for scale, which then refers to a second tree for shape, which finds the encoding
    interval_map < double,lower_interval_map,total_enricher > recursive_interval_map = interval_map < double,lower_interval_map,total_enricher>();

    // When writing the binary file, this vector is appended, so the position of each sequence is stored
    vector<CompressedRunnieIndex> indexes;

    /// Methods ///
    CompressedRunnieWriter(path file_path, path params_path);
    void load_parameters();
    void build_recursive_interval_tree();
    uint8_t fetch_encoding(double scale, double shape);

    void write_sequence(RunnieSequence& sequence);
    void write_sequence_block(RunnieSequence& sequence);
    void write_encoding_block(RunnieSequence& sequence);
    void write_index(CompressedRunnieIndex& index);
    void write_indexes();
};


#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEWRITER_HPP
