#include "CigarParser.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <iostream>
#include <stdexcept>
#include <experimental/filesystem>

using std::string;
using std::to_string;
using std::cout;
using std::vector;
using std::array;
using std::runtime_error;
using std::abs;
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
 * #define BAM_CIGAR_STR   "MIDNSHP=XB"
 * #define BAM_CIGAR_SHIFT 4
 * #define BAM_CIGAR_MASK  0xf
 * #define BAM_CIGAR_TYPE  0x3C1A7
 **/


CigarParser::CigarParser(path bam_path){
    this->bam_path = bam_path.string();
    this->bam_file = NULL;
    this->bam_index = NULL;
    this->bam_iterator = NULL;

    // bam file
    if ((this->bam_file = hts_open(this->bam_path.string().c_str(), "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file" + string(this->bam_path));
    }

    // bam index
    if ((this->bam_index = sam_index_load(this->bam_file, this->bam_path.string().c_str())) == 0) {
        throw runtime_error("ERROR: Cannot open index for bam file " + string(this->bam_path) + "\n");
    }

    this->bam_header = sam_hdr_read(this->bam_file);
    this->alignment = bam_init1();
}


void CigarParser::initialize_region(string ref_name, uint64_t start, uint64_t stop){
    this->ref_name = ref_name;

    // Find the number ID for this contig/chromosome/region
    this->ref_id = bam_name2id(bam_header, ref_name.c_str());

    // sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
    this->region_start = start;
    this->region_stop = stop;
    this->bam_iterator = sam_itr_queryi(this->bam_index, this->ref_id, this->region_start, this->region_stop);

    if (this->bam_iterator == 0) {
        throw runtime_error("ERROR: Cannot open iterator for region "
                            + this->ref_name + ":" + to_string(this->region_start) + ":" + to_string(this->region_stop)
                            + " for bam file " + string(this->bam_path) + "\n");
    }
}


void CigarParser::print_region(string ref_name, uint64_t start, uint64_t stop){
    // Initialize indexes/pointers
    this->initialize_region(ref_name, start, stop);

    // For tracking the position of the read
    uint64_t sequence_index;
    int64_t sequence_start;
    int64_t cigar_start;
    int8_t increment;

    // Iterate/fetch alignments
    int64_t result;
    while ((result = sam_itr_next(this->bam_index, this->bam_iterator, alignment)) >= 0) {
        this->aligned_segment = {};
        this->aligned_segment.ref_start_index = alignment->core.pos + 1;
        this->aligned_segment.read_length = alignment->core.l_qseq;
        this->aligned_segment.read_sequence = bam_get_seq(alignment);
        this->aligned_segment.read_name = bam_get_qname(alignment);
        this->aligned_segment.cigar = bam_get_cigar(alignment);
        this->aligned_segment.n_cigar = alignment->core.n_cigar;
        this->aligned_segment.reversal = bam_is_rev(alignment);

        sequence_start = 0;
        cigar_start = 0;
        increment = 1;

        cout << this->aligned_segment.to_string();
    }
}

//CigarParser::load_alignment(){
//
//}
