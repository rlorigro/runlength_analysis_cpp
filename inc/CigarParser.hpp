
#ifndef RUNLENGTH_ANALYSIS_CIGARPARSER_H
#define RUNLENGTH_ANALYSIS_CIGARPARSER_H

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "AlignedSegment.hpp"
#include <experimental/filesystem>
#include <utility>
#include <array>

using std::array;
using std::pair;
using std::experimental::filesystem::path;


class CigarParser{
public:
    // ---- Attributes ----
    path bam_path;

    // Samtools objects
    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    hts_itr_t* bam_iterator;
    bam1_t* alignment;

    // Class constants
    const string bases_forward = "=ACMGRSVTWYHKDBN";
    const string bases_reverse = "=TGKCYSBAWRDKHVN";
    const array <string, 2> bases = {bases_forward, bases_reverse};

    // Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
    const uint8_t BAM_SEQUENCE_SHIFT = 4;
    const uint8_t BAM_SEQUENCE_MASK = 15;     // 0b1111

    // Volatile metadata
    string ref_name;
    uint32_t ref_id;
    uint64_t region_start;
    uint64_t region_stop;
    AlignedSegment aligned_segment;

    // ---- Methods ----
    CigarParser(path bam_path);
    void initialize_region(string ref_name, uint64_t start, uint64_t stop);
    void print_region(string ref_name, uint64_t start, uint64_t stop);

    // Iterate cigars to return ref/read index for every match
    pair<uint64_t,uint64_t> get_next_match_index();

private:
    // ---- Attributes ----


    // ---- Methods ----


};


#endif //RUNLENGTH_ANALYSIS_CIGARPARSER_H
