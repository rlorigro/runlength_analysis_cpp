
#ifndef RUNLENGTH_ANALYSIS_KMER_HPP
#define RUNLENGTH_ANALYSIS_KMER_HPP

#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>

#include "IterativeSummaryStats.hpp"

using std::vector;
using std::array;
using std::deque;
using std::string;
using std::unordered_map;
using std::ofstream;


class KmerStats{
public:
    /// Attributes
    vector <IterativeSummaryStats <double> > qualities;
    uint8_t k;

    /// Methods
    KmerStats();
    KmerStats(uint8_t k);
    void update_quality(uint64_t kmer_index, double quality);
    string to_string();
};


class KmerConfusionStats{
public:
    /// Attributes
    // Initialize a 2d vector (confusion matrix) with zeros
    unordered_map <uint64_t, unordered_map<uint64_t, uint32_t> > confusion;
    uint8_t k;

    /// Methods
    KmerConfusionStats();
    KmerConfusionStats(uint8_t k);
    string to_string();
};


void operator+=(KmerStats& a, KmerStats& b);

void operator+=(KmerConfusionStats& a, KmerConfusionStats& b);

uint64_t kmer_to_index(deque<uint8_t>& kmer);

uint64_t kmer_to_index(deque<char>& kmer);

void index_to_kmer(deque<uint8_t>& kmer, uint64_t index, uint8_t k);

string kmer_index_to_string(uint64_t kmer_index, uint8_t k);


#endif //RUNLENGTH_ANALYSIS_KMER_HPP
