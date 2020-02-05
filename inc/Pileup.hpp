#ifndef RUNLENGTH_ANALYSIS_PILEUP_HPP
#define RUNLENGTH_ANALYSIS_PILEUP_HPP

#include <utility>
#include <vector>
#include <iostream>
#include <cassert>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <experimental/filesystem>

using std::cout;
using std::string;
using std::vector;
using std::ostream;
using std::deque;
using std::unordered_map;
using std::unordered_set;
using std::experimental::filesystem::path;


class Pileup{
public:
    /// Attributes ///
    static const size_t BASE = 0;
    static const size_t REVERSAL = 1;
    constexpr static const float EMPTY = 6.0;
    constexpr static const float INSERT_CODE = 5.0;
    constexpr static const float DELETE_CODE = 4.0;

    size_t n_alignments = 0;
    size_t n_channels;
    size_t max_observed_depth = 0;
    unordered_map <int64_t, vector <vector <vector <float> > > > inserts;
    vector <vector <vector <float> > > pileup;
    vector<uint16_t> coverage_per_position;

    /// Methods ///
    Pileup();
    Pileup(size_t n_channels, size_t region_size, size_t maximum_depth, vector<float>& default_read_data);
};


#endif //RUNLENGTH_ANALYSIS_PILEUP_HPP
