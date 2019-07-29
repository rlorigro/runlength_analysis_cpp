
#ifndef RUNLENGTH_ANALYSIS_CIGARPARSER_H
#define RUNLENGTH_ANALYSIS_CIGARPARSER_H

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "AlignedSegment.hpp"
#include <experimental/filesystem>
#include <string>
#include <array>

using std::array;
using std::string;
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

    // Volatile metadata
    bool valid_region;
    string ref_name;
    uint32_t ref_id;
    uint64_t region_start;
    uint64_t region_stop;
    AlignedSegment aligned_segment;

    /// Methods ///
    explicit BamReader(path bam_path);
    void initialize_region(string& ref_name, uint64_t start, uint64_t stop);
    void load_alignment(AlignedSegment& aligned_segment, bam1_t* alignment, bam_hdr_t* bam_header);
//    string get_read_base(AlignedSegment& aligned_segment, int64_t read_start_index, int64_t i);
//    void get_cigar(Cigar& cigar, AlignedSegment& aligned_segment, int64_t i);
//    int64_t infer_reference_stop_position_from_alignment(AlignedSegment& aligned_segment);
//    uint64_t get_ref_index_increment(Cigar& cigar);
//    uint64_t get_read_index_increment(Cigar& cigar);
    bool next_alignment(AlignedSegment& aligned_segment);
//    void print_region(string ref_name, uint64_t ref_start, uint64_t ref_stop);

    // Iterate cigars to return ref/read index for every match
//    pair<uint64_t,uint64_t> get_next_match_index();

private:
    /// Attributes ///

    /// Methods ///

};


#endif //RUNLENGTH_ANALYSIS_CIGARPARSER_H
