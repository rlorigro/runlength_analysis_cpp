#ifndef RUNLENGTH_ANALYSIS_PILEUP_HPP
#define RUNLENGTH_ANALYSIS_PILEUP_HPP

#include <utility>
#include <vector>
#include <iostream>
#include <cassert>
#include <deque>
#include <unordered_map>
#include <experimental/filesystem>

using std::cout;
using std::vector;
using std::ostream;
using std::deque;
using std::unordered_map;
using std::experimental::filesystem::path;


class Pileup{
public:
    /// Attributes ///
    size_t n_channels;
    unordered_map <int64_t, vector <vector <vector <float> > > > inserts;
    vector <vector <vector <float> > > pileup;

    /// Methods ///
    Pileup(size_t n_channels);
};

#endif //RUNLENGTH_ANALYSIS_PILEUP_HPP
