#include "hts.h"
#include "sam.h"
#include "bgzf.h"
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <stdexcept>

using std::string;
using std::cout;
using std::vector;
using std::runtime_error;

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


//TODO: Make a new class for handling alignments
//TODO: Practice threading, and see if bam data can be cast as atomic

void read_bam_file(char* bam_path) {
    samFile *in = hts_open(bam_path, "r");

    if (in == NULL) {
        throw runtime_error("ERROR: Cannot open bam file" + string(bam_path));
    }

    bam_hdr_t* bam_header = sam_hdr_read(in);
    bam1_t* alignment = bam_init1();

    int64_t pos;            // Left most position of alignment
    char* chr;              // Contig name (chromosome)
    int64_t len;            // Length of the read.
    uint8_t* seq;           // DNA sequence
    char* read_name;
    uint32_t* cigar;
    uint32_t n_cigar;
    uint64_t cigar_code;
    uint64_t cigar_length;

    while (sam_read1(in, bam_header, alignment) > 0) {
        pos = alignment->core.pos + 1;
        chr = bam_header->target_name[alignment->core.tid];
        len = alignment->core.l_qseq;
        seq = bam_get_seq(alignment);
        read_name = bam_get_qname(alignment);
        cigar = bam_get_cigar(alignment);
        n_cigar = alignment->core.n_cigar;

        cout << chr << " " << pos << " " << len << " " << read_name << " " << " " << n_cigar << "\n";

        for (uint32_t i=0; i < n_cigar; i++){
            cigar_code = cigar[i] & BAM_CIGAR_MASK;
            cigar_length = cigar[i] >> BAM_CIGAR_SHIFT;

            if (i < 10){
                cout << "(" << cigar_code << ", " << cigar_length << "), ";
            }
        }

        cout << "\n";

    }
}


int main() {
    string bam_path = "/home/ryan/data/nanopore/human/NA12878.np.chr3.100kb.0.bam";
    char* bam_path_chars = const_cast<char*>(bam_path.c_str());

    read_bam_file(bam_path_chars);
}
