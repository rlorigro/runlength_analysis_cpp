#ifndef RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
#define RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP

#include "AlignedSegment.hpp"
#include "FastaReader.hpp"
#include "BamReader.hpp"
#include "Runlength.hpp"
#include <utility>
#include <iostream>
#include <cassert>
#include <deque>
#include <experimental/filesystem>

using std::cout;
using std::deque;
using std::unordered_map;
using std::experimental::filesystem::path;

class PileupGenerator{
public:
    /// Attributes ///
    path bam_path;
    path ref_fasta_path;
    path reads_fasta_path;
    BamReader bam_reader;
    FastaReader ref_fasta_reader;
    FastaReader reads_fasta_reader;

    unordered_map <int64_t, vector <vector <string> > > insert_columns;
    vector <vector <string> > pileup;

    /// Methods ///
    PileupGenerator(path bam_path, path ref_fasta_path, path reads_fasta_path);
    void fetch_region(Region region);
    int64_t find_depth_index(int64_t start_index);

private:
    deque <pair <int64_t, int64_t> > lowest_free_index_per_depth;

};

#endif //RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
