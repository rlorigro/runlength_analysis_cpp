#ifndef RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP
#define RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

using std::map;
using std::string;
using std::ifstream;
using std::experimental::filesystem::path;


class MarginPolishReader{
public:
    /// Attributes ///
    path directory_path;
    map<string,path> read_paths;

    /// Methods ///
    MarginPolishReader(path directory_path);

    // Read all the file paths in directory and store in a map of names:path where `name` is derived from the filename
    void index();

    // Return a copy of the read indexes
    map<string,path> get_index();

private:
    /// Attributes ///

    /// Methods ///
};

#endif //RUNLENGTH_ANALYSIS_MARGINPOLISHREADER_HPP
