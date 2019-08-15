
#ifndef RUNLENGTH_ANALYSIS_RUNNIE_HPP
#define RUNLENGTH_ANALYSIS_RUNNIE_HPP

#include "Base.hpp"
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
using std::ostream;
using std::experimental::filesystem::path;


class RunnieIndex {
public:
    path file_path;
    uint64_t byte_index;
    uint64_t length;

    RunnieIndex(path file_path, uint64_t byte_index, uint64_t length);
    ostream& operator<<(ostream& s);
};


class RunnieSequence{
public:
    /// Attributes ///
    string sequence;
    vector <double> scales;
    vector <double> shapes;

    /// Methods ///
};


class RunnieReader{
public:
    /// Attributes ///
    path directory_path;
    unordered_map <string, RunnieIndex> read_indexes;

    /// Methods ///
    RunnieReader(path directory_path);

    void index();
    void index_file(path file_path);
    void fetch_read(string& read_name);
    void update_index(path& file_path,
        ifstream& file,
        string& line,
        uint64_t& l,
        string& read_name,
        uint64_t& byte_index,
        uint64_t& read_length,
        uint64_t& current_header_line_index,
        uint64_t& current_header_byte_index);
};


#endif //RUNLENGTH_ANALYSIS_RUNNIE_HPP
