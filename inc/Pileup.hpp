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
#include "Kmer.hpp"

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


class PileupReadKmerIterator {
public:
    vector <unordered_set <size_t> > kmer_indexes;

    deque <uint8_t> current_kmer;   // for efficiency of updating the kmer
    deque <uint64_t> window_kmer_indexes;  // This is all the kmers in the current window, possibly including inserts
    deque <uint64_t> n_operations;  // This is how we know how many kmers to cycle in the deque, it will always be
                                    // exactly the window size, where each element is n kmers per ref index
    size_t depth_index;
    size_t width_index;
    size_t window_size;
    size_t middle_index;
    uint8_t k;

    PileupReadKmerIterator(size_t window_size, size_t depth_index, uint8_t k);
    void step(Pileup& pileup, unordered_map<uint64_t, unordered_set <uint64_t> >& middle_kmers);
    void pop_left_kmers(Pileup& pileup);
    void push_right_kmer(uint8_t base);
    void update_middle_kmers(unordered_map<uint64_t, unordered_set <uint64_t> >& middle_kmers);
};


class PileupKmerIterator {
public:
    vector <PileupReadKmerIterator> read_iterators;     // Iterating the pileup is simplified by abstracting each row
                                                        // as their own iterators

    unordered_map <uint64_t, unordered_set <uint64_t> > middle_kmers;   // The kmers and their coverages in the middle
                                                                        // of the window

    vector <unordered_set <size_t> > kmer_indexes;      // This is all the kmers in the current window, possibly
                                                        // including inserts, and the indexes of their supporting reads

    vector <vector <uint64_t> > kmer_confusion_matrix;  // This represents the frequency of overlap between kmers, for
                                                        // all possible pairs of kmers

    size_t width_index = 0;
    size_t window_size;
    uint8_t k;

    PileupKmerIterator(Pileup& pileup, size_t window_size, uint8_t k);
    PileupKmerIterator(Pileup& pileup, size_t window_size, uint8_t k, KmerStats& kmer_stats);
    string to_string();
    void step(Pileup& pileup);
    void update_coverage_stats(Pileup& pileup, KmerStats& kmer_stats);
    void update_confusion_stats(PileupKmerIterator& ref_pileup_iterator, KmerConfusionStats& kmer_confusion_stats);
};


#endif //RUNLENGTH_ANALYSIS_PILEUP_HPP
