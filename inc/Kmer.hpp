
#ifndef RUNLENGTH_ANALYSIS_KMER_HPP
#define RUNLENGTH_ANALYSIS_KMER_HPP

//#include <Pileup.hpp>
#include <vector>
#include <deque>
#include <string>
#include "IterativeSummaryStats.hpp"

using std::vector;
using std::deque;
using std::string;


class KmerStats{
public:
    /// Attributes
    vector <IterativeSummaryStats <double> > qualities;
    vector <vector <uint64_t> > confusion;
    uint8_t k;

    /// Methods
    KmerStats();
    KmerStats(uint8_t k);
    void update_quality(uint64_t kmer_index, double quality);
    string to_string();
};

void operator+=(KmerStats& a, KmerStats& b);

uint64_t kmer_to_index(deque<uint8_t>& kmer);

void index_to_kmer(deque<uint8_t>& kmer, uint64_t index, uint8_t k);

string kmer_index_to_string(uint64_t kmer_index, uint8_t k);

#endif //RUNLENGTH_ANALYSIS_KMER_HPP
