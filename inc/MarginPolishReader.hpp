#ifndef RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP
#define RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP

#include "CoverageSegment.hpp"
//#include "Base.hpp"
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


class MarginPolishReader{
public:
    /// Attributes ///
    path directory_path;
    unordered_map<string,path> file_paths;

    /// Methods ///
    MarginPolishReader(path directory_path);

    // Read all the file paths in directory and store in a map of names:path where `name` is derived from the filename
    void index();

    // Return a copy of the read indexes
    unordered_map<string,path> get_index();
    void set_index(unordered_map<string, path>& file_paths);

    void read_file(CoverageSegment& mp_segment, path& file_path);
    void fetch_read(CoverageSegment& mp_segment, string& read_name);
    void parse_coverage_string(CoverageSegment& mp_segment, string& line);
    bool parse_reversal_string(string reversal_string);
    void read_consensus_sequence_from_file(CoverageSegment& mp_segment, path& file_path);
    void fetch_consensus_sequence(CoverageSegment& mp_segment, string& read_name);

private:
    /// Attributes ///

    /// Methods ///
};


#endif //RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP
