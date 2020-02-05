
#ifndef RUNLENGTH_ANALYSIS_CIGARKMER_HPP
#define RUNLENGTH_ANALYSIS_CIGARKMER_HPP

#include <map>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>

#include "Base.hpp"
#include "Kmer.hpp"
#include "AlignedSegment.hpp"
#include "IterativeSummaryStats.hpp"

using std::map;
using std::vector;
using std::array;
using std::deque;
using std::string;
using std::unordered_map;
using std::ofstream;


class CigarKmer {
public:
    /// Attributes
    deque <char> kmer;
    deque <map <uint8_t, uint64_t> > cigar_counts;
    map <uint8_t, uint64_t> total_cigar_counts;
    uint8_t k;

    /// Methods
    CigarKmer(uint8_t k);
    void update(Cigar& cigar, uint8_t base);
    uint64_t to_index();
};


class KmerIdentities{
public:
    /// Attributes
    vector <map <uint8_t, uint64_t> > cigar_counts_per_kmer;
    uint8_t k;

    /// Methods
    KmerIdentities();
    KmerIdentities(uint8_t k);
    void update(uint64_t kmer_index, array<uint64_t, 9>& cigar_counts);
    void update(CigarKmer& cigar_kmer);
    void write_to_output(ofstream& output);
};


void operator+=(KmerIdentities& a, KmerIdentities& b);



#endif //RUNLENGTH_ANALYSIS_CIGARKMER_HPP
