
#ifndef RUNLENGTH_ANALYSIS_RUNNIEINDEX_HPP
#define RUNLENGTH_ANALYSIS_RUNNIEINDEX_HPP

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

#endif //RUNLENGTH_ANALYSIS_RUNNIEINDEX_HPP
