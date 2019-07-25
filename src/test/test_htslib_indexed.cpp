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


void read_bam_file(char* bam_path) {
    string ref_name = "synthetic_ref_0";

    samFile *in = NULL;

    // bam file
    if ((in = hts_open(bam_path, "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file" + string(bam_path));
    }

    bam_hdr_t* bam_header = sam_hdr_read(in);
    bam1_t* alignment = bam_init1();


    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;

    // bam index
    if ((idx = sam_index_load(in, bam_path)) == 0) {
        throw runtime_error("ERROR: Cannot open index for bam file " + string(bam_path) + "\n");
    }

    //    sam_hdr_name2tid(sam_hdr_t *h, const char *ref);
    int32_t tid = bam_name2id(bam_header, ref_name.c_str());

    // sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
    uint64_t start = 0;
    uint64_t end = 500;
    iter = sam_itr_queryi(idx, tid, start, end);
    if (iter == 0) {
        throw runtime_error("ERROR: Cannot open iterator for region " + to_string(tid) + ":" + to_string(start) + ":" + to_string(end) + " for bam file " + string(bam_path) + "\n");
    }

    int64_t ref_start_index;     // Left most position of alignment
//    char* ref_name;              // Reference contig name (usually chromosome)
    int64_t read_length;         // Length of the read.
    uint8_t* read_sequence;      // DNA sequence
    char* read_name;
    uint32_t* cigar;
    uint32_t n_cigar;
    uint64_t cigar_code;
    uint64_t cigar_length;
    bool reversal;
    string bases_forward = "=ACMGRSVTWYHKDBN";
    string bases_reverse = "=TGKCYSBAWRDKHVN";
    array <string, 2> bases = {bases_forward, bases_reverse};

    // Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
    uint8_t BAM_SEQUENCE_SHIFT = 4;
    uint8_t BAM_SEQUENCE_MASK = 15;     // 0b1111

    // For tracking the position of the read
    uint64_t sequence_index;
    int64_t sequence_start;
    int64_t cigar_start;
    int8_t increment;

    // Fetch alignments
    int64_t result;
    while ((result = sam_itr_next(in, iter, alignment)) >= 0) {
        ref_start_index = alignment->core.pos + 1;
        ref_name = bam_header->target_name[alignment->core.tid];
        read_length = alignment->core.l_qseq;
        read_sequence = bam_get_seq(alignment);
        read_name = bam_get_qname(alignment);
        cigar = bam_get_cigar(alignment);
        n_cigar = alignment->core.n_cigar;
        reversal = bam_is_rev(alignment);
        sequence_start = 0;
        cigar_start = 0;
        increment = 1;

        string reversal_string = reversal ? "R" : "F";

        cout << ref_name << " " << ref_start_index << " " << reversal_string << " " << read_length << " " << read_name << " " << " " << n_cigar << "\n";

        // Index the sequence in its F direction (not alignment direction)
        if (reversal){
            sequence_start = read_length - 1;
            cigar_start = n_cigar - 1;
            increment = -1;
        }

        // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
        for (int64_t i=sequence_start; llabs(sequence_start-i)<read_length; i+=increment){
            sequence_index = i/2;
            if (i%2 == 0){
                // Perform bit shift and decode using the standard or complemented base map
                cout << bases[reversal][read_sequence[sequence_index] >> BAM_SEQUENCE_SHIFT];
            }
            else {
                // Perform bit mask and decode using the standard or complemented base map
                cout << bases[reversal][read_sequence[sequence_index] & BAM_SEQUENCE_MASK];
            }
        }
        cout << "\n";

        // Fetch the 4 bit cigar operation and 4 bit cigar length from the
        for (int64_t i=cigar_start; llabs(cigar_start-i)<n_cigar; i+=increment){
            cigar_code = cigar[i] & BAM_CIGAR_MASK;
            cigar_length = cigar[i] >> BAM_CIGAR_SHIFT;

            cout << "(" << cigar_code << ", " << cigar_length << "), ";
        }
        cout << "\n";
    }
}


int main() {
    // Find absolute path to test file within repo
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_input_path = "/data/test/";
    path filename = "test_alignable_sequences_non_RLE_VS_test_alignable_reference_non_RLE.sorted.bam";
    path absolute_input_path = project_directory / relative_input_path / filename;

    char* bam_path_chars = const_cast<char*>(absolute_input_path.c_str());

    read_bam_file(bam_path_chars);
}
