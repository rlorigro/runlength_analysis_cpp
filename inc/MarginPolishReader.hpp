#ifndef RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP
#define RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP

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


class CoverageElement{
public:
    /// Attributes ///
    string base;
    uint16_t length;
    bool reversal;
    double weight;

    CoverageElement(string base, uint16_t length, bool reversal, double weight);
    static const string reversal_string;

    /// Methods ///
    string to_string();
};


class MarginPolishSegment{
public:
    /// Attributes ///
    string name;
    vector <vector <CoverageElement> > coverage_data;
    string sequence;

    /// Methods ///
    void print();
};


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

    void read_file(MarginPolishSegment& mp_segment, path& file_path);
    void fetch_read(MarginPolishSegment& mp_segment, string& read_name);
    void parse_coverage_string(MarginPolishSegment& mp_segment, string& line);
    bool parse_reversal_string(string reversal_string);
    void read_consensus_sequence_from_file(MarginPolishSegment& mp_segment, path& file_path);
    void fetch_consensus_sequence(MarginPolishSegment& mp_segment, string& read_name);

private:
    /// Attributes ///

    /// Methods ///
};


#endif //RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP
