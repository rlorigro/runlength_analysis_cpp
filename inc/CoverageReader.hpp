
#ifndef RUNLENGTH_ANALYSIS_COVERAGEREADER_H
#define RUNLENGTH_ANALYSIS_COVERAGEREADER_H

#include "CoverageSegment.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

using std::unordered_map;
using std::vector;
using std::string;
using std::ifstream;
using std::experimental::filesystem::path;


class CoverageReader {
public:
    /// Attributes ///
    path directory_path;
    unordered_map<string,path> file_paths;

    /// Methods ///
    void index();
    unordered_map<string,path> get_index();
    void set_index(unordered_map<string, path>& file_paths);

    void fetch_read(CoverageSegment& mp_segment, string& read_name);
    void fetch_consensus_sequence(CoverageSegment& mp_segment, string& read_name);
};


#endif //RUNLENGTH_ANALYSIS_COVERAGEREADER_H
