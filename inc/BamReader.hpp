#ifndef RUNLENGTH_ANALYSIS_CIGARPARSER_H
#define RUNLENGTH_ANALYSIS_CIGARPARSER_H

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "RunlengthSequenceElement.hpp"
#include "AlignedSegment.hpp"
#include "FastaReader.hpp"
#include <experimental/filesystem>
#include <unordered_map>
#include <map>
#include <string>
#include <vector>
#include <array>

using std::unordered_map;
using std::map;
using std::array;
using std::string;
using std::vector;
using std::experimental::filesystem::path;


/**
 * ****************************
 * *** CIGAR related macros ***
 * ****************************
 *
 * #define BAM_CMATCH      0
 * #define BAM_CINS        1
 * #define BAM_CDEL        2
 * #define BAM_CREF_SKIP   3
 * #define BAM_CSOFT_CLIP  4
 * #define BAM_CHARD_CLIP  5
 * #define BAM_CPAD        6
 * #define BAM_CEQUAL      7
 * #define BAM_CDIFF       8
 * #define BAM_CBACK       9
 *
**/




class CigarStats {
public:
    uint64_t n_matches;
    uint64_t n_mismatches;
    uint64_t n_inserts;
    uint64_t n_deletes;

    map <uint8_t, map <uint64_t, uint64_t> > cigar_lengths;

    CigarStats();
    string to_string();
    double calculate_identity();
    void update_lengths(Cigar& cigar);
};


void operator+=(CigarStats& cigar_stats_a, CigarStats& cigar_stats_b);


// For occasional convenience
class Region {
public:
    string name;
    uint64_t start;
    uint64_t stop;

    Region(string name, uint64_t start, uint64_t stop);
    Region();
    string to_string();
};


class BamReader{
public:
    /// Attributes ///
    path bam_path;

    // Samtools objects
    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    hts_itr_t* bam_iterator;
    bam1_t* alignment;

    // Bit operations
    static const uint16_t secondary_mask = 256;
    static const uint16_t supplementary_mask = 2048;

    // Volatile metadata
    bool valid_region;
    string ref_name;
    int ref_id;
    uint64_t region_start;
    uint64_t region_stop;
    AlignedSegment aligned_segment;

    /// Methods ///
    BamReader(path bam_path);
    BamReader();
    ~BamReader();
    void initialize_region(string& ref_name, uint64_t start, uint64_t stop);
    void load_alignment(AlignedSegment& aligned_segment, bam1_t* alignment, bam_hdr_t* bam_header);
    bool next_alignment(AlignedSegment& aligned_segment,
                        uint16_t map_quality_cutoff=0,
                        bool filter_secondary=true,
                        bool filter_supplementary=false);

private:
    /// Attributes ///

    /// Methods ///

};


void chunk_sequence(vector<Region>& regions, string read_name, uint64_t chunk_size, uint64_t length);

#endif //RUNLENGTH_ANALYSIS_CIGARPARSER_H
