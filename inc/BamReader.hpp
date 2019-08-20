#ifndef RUNLENGTH_ANALYSIS_CIGARPARSER_H
#define RUNLENGTH_ANALYSIS_CIGARPARSER_H

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "AlignedSegment.hpp"
#include "FastaReader.hpp"
#include <experimental/filesystem>
#include <unordered_map>
#include <string>
#include <vector>
#include <array>

using std::unordered_map;
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
    uint16_t secondary_mask;

    // Volatile metadata
    bool valid_region;
    string ref_name;
    uint32_t ref_id;
    uint64_t region_start;
    uint64_t region_stop;
    AlignedSegment aligned_segment;

    /// Methods ///
    BamReader(path bam_path);
    BamReader();
    void initialize_region(string& ref_name, uint64_t start, uint64_t stop);
    void load_alignment(AlignedSegment& aligned_segment, bam1_t* alignment, bam_hdr_t* bam_header);
    bool next_alignment(AlignedSegment& aligned_segment, uint16_t map_quality_cutoff=0, bool filter_secondary=false);

private:
    /// Attributes ///

    /// Methods ///

};


void chunk_sequence(vector<Region>& regions, string read_name, uint64_t chunk_size, uint64_t length);

vector<Region> chunk_sequences_from_fasta_index_into_regions(unordered_map<string,FastaIndex> index_map, uint64_t chunk_size);

#endif //RUNLENGTH_ANALYSIS_CIGARPARSER_H
