
#ifndef RUNLENGTH_ANALYSIS_SHASTAREADER_HPP
#define RUNLENGTH_ANALYSIS_SHASTAREADER_HPP

#include "CoverageSegment.hpp"
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <experimental/filesystem>

using std::move;
using std::cout;
using std::cerr;
using std::flush;
using std::ifstream;
using std::to_string;
using std::getline;
using std::replace;
using std::runtime_error;
using std::unordered_map;
using std::experimental::filesystem::directory_iterator;
using std::experimental::filesystem::path;


class ShastaReader{
public:
    /// Attributes ///
    path directory_path;
    unordered_map<string,path> file_paths;
    bool store_length_consensus;
    bool store_coverage_data;

    /// Methods ///
    ShastaReader(path directory_path, bool store_length_consensus=false, bool store_coverage_data=true);

    // Read all the file paths in directory and store in a map of names:path where `name` is derived from the filename
    void index();

    // Return a copy of the read indexes
    unordered_map<string,path> get_index();
    void set_index(unordered_map<string, path>& file_paths);

    void read_file(CoverageSegment& segment, path& file_path);
    void fetch_read(CoverageSegment& segment, string& read_name);
    void parse_coverage_string(CoverageSegment& segment, string& line);
    size_t parse_consensus(CoverageSegment& segment, string& line);
    bool parse_reversal_string(string reversal_string);
    void read_consensus_sequence_from_file(CoverageSegment& segment, path& file_path);
    void fetch_consensus_sequence(CoverageSegment& segment, string& read_name);

private:
    /// Attributes ///

    /// Methods ///
};

#endif //RUNLENGTH_ANALYSIS_SHASTAREADER_HPP
