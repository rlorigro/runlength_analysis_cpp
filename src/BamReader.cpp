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


BamReader::BamReader(path bam_path){
    this->bam_path = bam_path.string();
    this->bam_file = nullptr;
    this->bam_index = nullptr;
    this->bam_iterator = nullptr;
    this->alignment = bam_init1();

    this->valid_region = false;
    this->ref_name = "";
    this->region_start = -1;
    this->region_stop = -1;

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

    this-> valid_region = true;
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
    aligned_segment.cigars = bam_get_cigar(alignment);
    aligned_segment.n_cigar = alignment->core.n_cigar;
    aligned_segment.reversal = bam_is_rev(alignment);
}


bool BamReader::next_alignment(AlignedSegment& aligned_segment){
    ///
    /// Iterate the alignments within a region, returning references
    ///

    if (not this->valid_region) {
        throw runtime_error("ERROR: BAM reader must be initialized with `initialize_region()` before iterating");
    }

    // Iterate/fetch alignments
    int64_t result;
    if ((result = sam_itr_next(this->bam_file, this->bam_iterator, alignment)) >= 0) {
        aligned_segment = {};
        load_alignment(aligned_segment, this->alignment, this->bam_header);
    }
    else{
        this->valid_region = false;
    }

    return this->valid_region;
}


//void BamReader::print_region(string ref_name, uint64_t ref_start, uint64_t ref_stop){
//    ///
//    /// Iterate all the alignments within a region and print their aligned cigar subsegments
//    ///
//
//    // Initialize indexes/pointers
//    this->initialize_region(ref_name, ref_start, ref_stop);
//
//    // For tracking the position of the read/ref
//    int64_t ref_start_index;
//    int64_t read_start_index;
//    int64_t cigar_start;
//    int8_t increment;
//    int64_t ref_increment;
//    int64_t read_increment;
//    int64_t ref_index;
//    int64_t read_index;
//
//    string ref_sequence = "ACCTTGCGATGCTAGCATAGCATGGCTCATGAATGCGATCCGATTGCAGTCCGAATGCATTGCTACGCAGTGCATAGGCTAGCTCAGACTGCTAGCTAGGCTACTAGCATGCCTAGTCAGTGACTAGCCTAGCGTTAGTCGATCATATCAGCGTACTCATCGATGCAGCACATGCATGCTATGTCTAGTACTACGGTACGATTATTCGATTCGCCGATGATCTAGCATGCGTACTGCTAGATGCTATGCATGCGTATGATATCTGATGCATGTCAGTTATGCATATTCGATATGTACTAGTTGCAGTCATGTGCATTATGCAGCTATTATTACGCTGAGTGCATAGCATGTCTGTCGCTAGCTAGATCGTAGCATGATCAGCATCATTCATGCATTCGTAGCATGCTAATGTCATTATAGCTATCGGCATTATCAGTAGCATCTAGCATAGATCTCAGTACGTAGTATCTATCAGTCAGTAGCTAGTCGATACGTATCGTTGCATATGCTACGATGATGCTATGCATGGCAATGCATTGCACTAGCTACGTATACTATGCATGTCAGTAGATGCTATACTGATCGAATCTTGATCTAGTAGCCTAGTAGCACTGACTGCATGTCAGTACGTAGCTACGTTCGTATACGATCATCTAGATCATGCTAGCATGCATGCATATATGTGACTGATGCTGATGCCGGATTATCGCGTATCGATCGATCATCATGATCATGATGATGCTGCATCAGAATACGCTGACTGACATCGACGACTGCATCGCGACTGCATCGGCAGCTAGCATGCGCGATGCATGCATCGTACTGCATGCAGTCGATCGATGCACGATGCATGCATGCATCGAGTATAGCCGGATTAGCTACTGAGCGATTTATCTCTGAGGAGATCTCGATCGTAGCATGCTGCGCATCTGCTAATGTCGGATGCTAGCGCTAGCTGCTTAGCTCATATTACGTATCTGATCTGATTCGATGCATGCATTATCGATTCGTATTAGCATCGTACGTAGCTATGCATTCGTAGCTAGCATCGTAGCTGAGCGATGCTATGCGCTAGCTTAGCGATGCTGCCGATCGTAGCGTATCAGAGTCGATCGTAGCTAGCTACGCGTACTAGCTAGCTACGTTAGCGCTAGATGATCTAGGCGCTATTATCGAGAGTCTCTAGGCTACTGATATCTGAGCAGGAAGAGTCGATCGTATGCTGCTGCTAGTCGTACGTATCGTATCGATGCATGTCATGCATAGTATGCGATCGCATGCTACTGTGCTGATGCTAGCTAGCTAGTCGATGTCGTAGCGGCATGTAGCGTACGCGG";
//    string ref_subsegment;
//    string read_subsegment;
//
//    // Iterate/fetch alignments
//    int64_t result;
//    while ((result = sam_itr_next(this->bam_file, this->bam_iterator, alignment)) >= 0) {
//        cout << "\nRESULT: " << result << "\n";
//        this->aligned_segment = {};
//        load_alignment(aligned_segment, this->alignment, this->bam_header);
//
//        Cigar cigar;
//        read_start_index = 0;
//        ref_start_index = aligned_segment.ref_start_index;
//        cigar_start = 0;
//        increment = 1;
//
//        cout << this->aligned_segment.to_string();
//
//        // Index the sequence in its F direction (not alignment direction)
//        if (aligned_segment.reversal){
//            read_start_index = aligned_segment.read_length - 1;
//            ref_start_index = infer_reference_stop_position_from_alignment();
//            cigar_start = aligned_segment.n_cigar - 1;
//            increment = -1;
//        }
//
//        // Iterate cigar
//        for (int64_t i=cigar_start; llabs(cigar_start-i)<aligned_segment.n_cigar; i+=increment) {
//            cigar = {};
//            get_cigar(cigar, aligned_segment, i);
//            cout << cigar.to_string();
//        }
//        cout << "\n";
//
//        // Iterate cigar and print aligned ref/read subsegments
//        ref_index = ref_start_index - 1;
//        read_index = read_start_index;
//
//        for (int64_t i=cigar_start; llabs(cigar_start-i)<aligned_segment.n_cigar; i+=increment) {
//            cigar = {};
//            ref_subsegment = {};
//            read_subsegment = {};
//
//            get_cigar(cigar, aligned_segment, i);
//            read_increment = get_read_index_increment(cigar)*increment;      // Reverse if necessary
//
//            cout << read_index << " " << read_index + read_increment << "\n";
//            for (int64_t i_subcigar=read_index; llabs(read_index - i_subcigar) < llabs(read_increment); i_subcigar+=increment) {
//                ref_subsegment += ref_sequence[ref_index];
//                ref_index += BamReader::cigar_ref_move[cigar.code]*increment;   // Reverse if necessary
//                read_subsegment += get_read_base(aligned_segment, read_index, i_subcigar);
//            }
//            cout << ref_subsegment << "\n";
//            cout << read_subsegment << "\n";
//
//            read_index += read_increment;
//        }
//        cout << "\n";
//    }
//}
