
#ifndef RUNLENGTH_ANALYSIS_FASTQREADER_HPP
#define RUNLENGTH_ANALYSIS_FASTQREADER_HPP

#include "IterativeSummaryStats.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::experimental::filesystem::path;


class FastqReader {
public:
    /// Attributes ///
    path file_path;

    /// Methods ///
    FastqReader(path file_path);

    void filter_by_quality(uint64_t minimum_average_quality);
    void parse_filter_candidate(
        string& name,
        string& sequence,
        string& qualities,
        IterativeSummaryStats<uint64_t>& stats,
        uint64_t minimum_average_quality,
        ofstream& output_file);
};

#endif //RUNLENGTH_ANALYSIS_FASTQREADER_HPP
