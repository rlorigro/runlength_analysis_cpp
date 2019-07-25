#include "BamReader.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
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


const array <string, 2> BamReader::bases = {"=ACMGRSVTWYHKDBN", "=TGKCYSBAWRDKHVN"};

const array <bool, 10> BamReader::cigar_ref_move = {true,      //BAM_CMATCH      0
                                                      false,     //BAM_CINS        1
                                                      true,      //BAM_CDEL        2
                                                      true,      //BAM_CREF_SKIP   3
                                                      false,     //BAM_CSOFT_CLIP  4
                                                      false,     //BAM_CHARD_CLIP  5
                                                      false,     //BAM_CPAD        6
                                                      true,      //BAM_CEQUAL      7
                                                      true,      //BAM_CDIFF       8
                                                      false};    //BAM_CBACK       9

const array <bool, 10> BamReader::cigar_read_move = {true,     //BAM_CMATCH      0
                                                       true,     //BAM_CINS        1
                                                       false,    //BAM_CDEL        2
                                                       false,    //BAM_CREF_SKIP   3
                                                       true,     //BAM_CSOFT_CLIP  4
                                                       false,    //BAM_CHARD_CLIP  5
                                                       false,    //BAM_CPAD        6
                                                       true,     //BAM_CEQUAL      7
                                                       true,     //BAM_CDIFF       8
                                                       false};   //BAM_CBACK       9


BamReader::BamReader(path bam_path){
    this->bam_path = bam_path.string();
    this->bam_file = nullptr;
    this->bam_index = nullptr;
    this->bam_iterator = nullptr;
    this->alignment = bam_init1();

    // bam file
    if ((this->bam_file = hts_open(this->bam_path.string().c_str(), "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file" + string(this->bam_path));
    }

    // bam index
    if ((this->bam_index = sam_index_load(this->bam_file, this->bam_path.string().c_str())) == 0) {
        throw runtime_error("ERROR: Cannot open index for bam file " + string(this->bam_path) + "\n");
    }

    this->bam_header = sam_hdr_read(this->bam_file);
}


void BamReader::initialize_region(string& ref_name, uint64_t start, uint64_t stop){
    this->ref_name = ref_name;

    // Find the ID for this contig/chromosome/region
    this->ref_id = bam_name2id(bam_header, ref_name.c_str());

    // sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
    this->region_start = start;
    this->region_stop = stop;
    this->bam_iterator = sam_itr_queryi(this->bam_index, this->ref_id, this->region_start, this->region_stop);

    if (this->bam_iterator == nullptr) {
        throw runtime_error("ERROR: Cannot open iterator for region "
                            + this->ref_name + ":" + to_string(this->region_start) + ":" + to_string(this->region_stop)
                            + " for bam file " + string(this->bam_path) + "\n");
    }
}


string BamReader::get_read_base(AlignedSegment& aligned_segment, int64_t read_start_index, int64_t i){
    ///
    /// Convert the compressed representation of an aligned sequence into a string
    ///

    uint8_t byte;
    uint8_t base_code;
    string base;

    // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
    if (llabs(read_start_index - i) < aligned_segment.read_length){
        uint64_t sequence_index = i/2;
        byte = aligned_segment.read_sequence[sequence_index];

        if (i%2 == 0){
            // Perform bit shift and decode using the standard or complemented base map
            base_code = byte >> BamReader::bam_sequence_shift;
            base = BamReader::bases[aligned_segment.reversal][base_code];
        }
        else {
            // Perform bit mask and decode using the standard or complemented base map
            base_code = byte & BamReader::bam_sequence_mask;
            base = BamReader::bases[aligned_segment.reversal][base_code];
        }
    }
    else{
        throw runtime_error("ERROR: read index out of bounds: " + aligned_segment.read_name + " " + to_string(i));
    }
    return base;
}


void BamReader::get_cigar(Cigar& cigar, AlignedSegment& aligned_segment, int64_t i){
    ///
    /// Decompress byte from BAM into the 2 components of a cigar operation: (code, length)
    ///

    cigar.code = aligned_segment.cigar[i] & BamReader::bam_cigar_mask;
    cigar.length = aligned_segment.cigar[i] >> BamReader::bam_cigar_shift;
}


uint64_t BamReader::get_read_index_increment(Cigar& cigar){
    ///
    /// Determine whether a cigar operation should result in a step along the READ sequence
    ///

    uint64_t increment = 0;
    if (BamReader::cigar_read_move[cigar.code]){
        increment = cigar.length;
    }
    return increment;
}


uint64_t BamReader::get_ref_index_increment(Cigar& cigar){
    ///
    /// Determine whether a cigar operation should result in a step along the REF sequence
    ///

    uint64_t increment = 0;
    if (BamReader::cigar_ref_move[cigar.code]){
        increment = cigar.length;
    }
    return increment;
}


void BamReader::load_alignment(AlignedSegment& aligned_segment, bam1_t* alignment, bam_hdr_t* bam_header){
    ///
    /// Load data from shitty samtools structs into a cpp object
    ///

    aligned_segment.ref_start_index = alignment->core.pos + 1;
    aligned_segment.ref_name = bam_header->target_name[alignment->core.tid];
    aligned_segment.read_length = alignment->core.l_qseq;
    aligned_segment.read_sequence = bam_get_seq(alignment);
    aligned_segment.read_name = bam_get_qname(alignment);
    aligned_segment.cigar = bam_get_cigar(alignment);
    aligned_segment.n_cigar = alignment->core.n_cigar;
    aligned_segment.reversal = bam_is_rev(alignment);
}


int64_t BamReader::infer_reference_stop_position_from_alignment(AlignedSegment& aligned_segment){
    ///
    /// Iterate the cigar of an alignment and find the reference stop position
    ///

    Cigar cigar;
    int64_t ref_stop_position = aligned_segment.ref_start_index;
    for (int64_t i=0; i<aligned_segment.n_cigar; i++) {
        cigar = {};
        get_cigar(cigar, aligned_segment, i);
        ref_stop_position += get_ref_index_increment(cigar);
    }

    return ref_stop_position - 1;
}

//TODO: convert this to a next() style function
void BamReader::print_region(string ref_name, uint64_t ref_start, uint64_t ref_stop){
    ///
    /// Iterate all the alignments within a region and print their aligned cigar subsegments
    ///

    // Initialize indexes/pointers
    this->initialize_region(ref_name, ref_start, ref_stop);

    // For tracking the position of the read/ref
    int64_t ref_start_index;
    int64_t read_start_index;
    int64_t cigar_start;
    int8_t increment;
    int64_t ref_increment;
    int64_t read_increment;
    int64_t ref_index;
    int64_t read_index;

    string ref_sequence = "ACCTTGCGATGCTAGCATAGCATGGCTCATGAATGCGATCCGATTGCAGTCCGAATGCATTGCTACGCAGTGCATAGGCTAGCTCAGACTGCTAGCTAGGCTACTAGCATGCCTAGTCAGTGACTAGCCTAGCGTTAGTCGATCATATCAGCGTACTCATCGATGCAGCACATGCATGCTATGTCTAGTACTACGGTACGATTATTCGATTCGCCGATGATCTAGCATGCGTACTGCTAGATGCTATGCATGCGTATGATATCTGATGCATGTCAGTTATGCATATTCGATATGTACTAGTTGCAGTCATGTGCATTATGCAGCTATTATTACGCTGAGTGCATAGCATGTCTGTCGCTAGCTAGATCGTAGCATGATCAGCATCATTCATGCATTCGTAGCATGCTAATGTCATTATAGCTATCGGCATTATCAGTAGCATCTAGCATAGATCTCAGTACGTAGTATCTATCAGTCAGTAGCTAGTCGATACGTATCGTTGCATATGCTACGATGATGCTATGCATGGCAATGCATTGCACTAGCTACGTATACTATGCATGTCAGTAGATGCTATACTGATCGAATCTTGATCTAGTAGCCTAGTAGCACTGACTGCATGTCAGTACGTAGCTACGTTCGTATACGATCATCTAGATCATGCTAGCATGCATGCATATATGTGACTGATGCTGATGCCGGATTATCGCGTATCGATCGATCATCATGATCATGATGATGCTGCATCAGAATACGCTGACTGACATCGACGACTGCATCGCGACTGCATCGGCAGCTAGCATGCGCGATGCATGCATCGTACTGCATGCAGTCGATCGATGCACGATGCATGCATGCATCGAGTATAGCCGGATTAGCTACTGAGCGATTTATCTCTGAGGAGATCTCGATCGTAGCATGCTGCGCATCTGCTAATGTCGGATGCTAGCGCTAGCTGCTTAGCTCATATTACGTATCTGATCTGATTCGATGCATGCATTATCGATTCGTATTAGCATCGTACGTAGCTATGCATTCGTAGCTAGCATCGTAGCTGAGCGATGCTATGCGCTAGCTTAGCGATGCTGCCGATCGTAGCGTATCAGAGTCGATCGTAGCTAGCTACGCGTACTAGCTAGCTACGTTAGCGCTAGATGATCTAGGCGCTATTATCGAGAGTCTCTAGGCTACTGATATCTGAGCAGGAAGAGTCGATCGTATGCTGCTGCTAGTCGTACGTATCGTATCGATGCATGTCATGCATAGTATGCGATCGCATGCTACTGTGCTGATGCTAGCTAGCTAGTCGATGTCGTAGCGGCATGTAGCGTACGCGG";
    string ref_subsegment;
    string read_subsegment;

    // Iterate/fetch alignments
    int64_t result;
    while ((result = sam_itr_next(this->bam_file, this->bam_iterator, alignment)) >= 0) {
        this->aligned_segment = {};
        load_alignment(aligned_segment, this->alignment, this->bam_header);

        Cigar cigar;
        read_start_index = 0;
        ref_start_index = aligned_segment.ref_start_index;
        cigar_start = 0;
        increment = 1;

        cout << this->aligned_segment.to_string();

        // Index the sequence in its F direction (not alignment direction)
        if (aligned_segment.reversal){
            read_start_index = aligned_segment.read_length - 1;
            ref_start_index = infer_reference_stop_position_from_alignment(aligned_segment);
            cigar_start = aligned_segment.n_cigar - 1;
            increment = -1;
        }

        // Iterate cigar
        for (int64_t i=cigar_start; llabs(cigar_start-i)<aligned_segment.n_cigar; i+=increment) {
            cigar = {};
            get_cigar(cigar, aligned_segment, i);
            cout << cigar.to_string();
        }
        cout << "\n";

        // Iterate cigar and print aligned ref/read subsegments
        ref_index = ref_start_index - 1;
        read_index = read_start_index;

        for (int64_t i=cigar_start; llabs(cigar_start-i)<aligned_segment.n_cigar; i+=increment) {
            cigar = {};
            ref_subsegment = {};
            read_subsegment = {};

            get_cigar(cigar, aligned_segment, i);
            read_increment = get_read_index_increment(cigar)*increment;      // Reverse if necessary

            cout << read_index << " " << read_index + read_increment << "\n";
            for (int64_t i_read=read_index; llabs(read_index-i_read)<llabs(read_increment); i_read+=increment) {
                ref_subsegment += ref_sequence[ref_index];
                ref_index += BamReader::cigar_ref_move[cigar.code]*increment;   // Reverse if necessary
                read_subsegment += get_read_base(aligned_segment, read_index, i_read);
            }
            cout << ref_subsegment << "\n";
            cout << read_subsegment << "\n";

            read_index += read_increment;
        }
        cout << "\n";
    }
}
