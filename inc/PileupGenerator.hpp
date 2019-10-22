#ifndef RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
#define RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP

#include "FastaReader.hpp"
#include "BamReader.hpp"
#include "Runlength.hpp"
#include "Pileup.hpp"
#include <utility>
#include <iostream>
#include <cassert>
#include <deque>
#include <experimental/filesystem>

using std::cout;
using std::ostream;
using std::deque;
using std::unordered_map;
using std::experimental::filesystem::path;

class PileupGenerator{
public:
    /// Attributes ///
    path bam_path;
    BamReader bam_reader;
    uint16_t maximum_depth;

    string null_character = "_";
    float null_value = -1;

    /// Methods ///
    PileupGenerator(path bam_path, uint16_t maximum_depth=80);
    void print_lowest_free_indexes();
    void print(Pileup& pileup);

    template <class T> void fetch_region(Region& region, FastaReader& ref_fasta_reader, T& sequence_reader);
    int64_t find_depth_index(int64_t start_index);
//    void get_base(string& read_base, Cigar& cigar, Coordinate& coordinate, SequenceElement& read_sequence);
    void parse_insert(Pileup& pileup, int64_t pileup_width_index, int64_t pileup_depth_index, AlignedSegment& aligned_segment, vector<float>& read_data);

private:
    deque <pair <int64_t, int64_t> > lowest_free_index_per_depth;
    vector <float> default_data_vector;
    vector <vector <float> > default_insert_column;
    vector <vector <vector <float> > > default_insert_pileup;
};

#endif //RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
