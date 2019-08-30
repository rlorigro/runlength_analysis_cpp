
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


class CompressedRunnieWriter {
public:
    /// Attributes ///
    path sequence_file_path;
    path index_file_path;
    path params_path;

    ofstream sequence_file;
    ofstream index_file;

    vector <pair <double,double> > scale_intervals;
    vector < vector <pair <double,double> > > shape_intervals;
    interval_map < double,lower_interval_map,total_enricher > recursive_interval_map = interval_map < double,lower_interval_map,total_enricher>();

    /// Methods ///
    CompressedRunnieWriter(path file_path, path params_path);
    void load_parameters();
    void build_recursive_interval_tree();
    uint8_t fetch_encoding(double scale, double shape);

    void write_sequence(RunnieSequence& sequence);
};


#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEWRITER_HPP
