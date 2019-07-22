
#ifndef RUNLENGTH_ANALYSIS_READ_H
#define RUNLENGTH_ANALYSIS_READ_H

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <string>

using std::string;


// Used primarily as a container inside CigarParser, which needs metadata to iterate the cigar/alignment data
class AlignedSegment{
public:
    // ---- Attributes ----
    int64_t ref_start_index;     // Left most position of alignment
    string ref_name;              // Reference contig name (usually chromosome)
    int64_t read_length;         // Length of the read.
    uint8_t* read_sequence;      // DNA sequence
    string read_name;
    uint32_t* cigar;
    uint32_t n_cigar;
    uint64_t cigar_code;
    uint64_t cigar_length;
    bool reversal;

    // ---- Methods ----
    string to_string();
};

#endif //RUNLENGTH_ANALYSIS_READ_H
