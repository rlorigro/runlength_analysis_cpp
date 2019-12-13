#include "Region.hpp"
#include "BamReader.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <experimental/filesystem>

using std::string;
using std::to_string;
using std::cout;
using std::vector;
using std::array;
using std::runtime_error;
using std::abs;
using std::free;
using std::experimental::filesystem::path;


CigarStats::CigarStats(){
    this->n_matches = 0;
    this->n_mismatches = 0;
    this->n_inserts = 0;
    this->n_deletes = 0;
}


//void CigarStats::update_lengths(Cigar& cigar) {
//    if (this->cigar_lengths.count() > 0){
//        cigar_stats.cigar_lengths[cigar.code][cigar.length]++;
//    }
//}


void operator+=(CigarStats& cigar_stats_a, CigarStats& cigar_stats_b){
    cigar_stats_a.n_matches += cigar_stats_b.n_matches;
    cigar_stats_a.n_mismatches += cigar_stats_b.n_mismatches;
    cigar_stats_a.n_inserts += cigar_stats_b.n_inserts;
    cigar_stats_a.n_deletes += cigar_stats_b.n_deletes;

    for (auto& [cigar_code, lengths]: cigar_stats_b.cigar_lengths){
        for (auto& [length, count]: lengths){
            cigar_stats_a.cigar_lengths[cigar_code][length] += count;
        }
    }
}


double CigarStats::calculate_identity() {
    return double(this->n_matches)/double(this->n_matches + this->n_mismatches + this->n_inserts + this->n_deletes);
}


string CigarStats::to_string(bool verbose){
    string s;

    s += "n_matches:\t" + std::to_string(this->n_matches) + "\n";
    s += "n_mismatches:\t" + std::to_string(this->n_mismatches) + "\n";
    s += "n_inserts:\t" + std::to_string(this->n_inserts) + "\n";
    s += "n_deletes:\t" + std::to_string(this->n_deletes) + "\n";

    double sum = this->n_matches + this->n_mismatches + this->n_inserts + this->n_deletes;
    double mismatch_rate = this->n_mismatches / sum;
    double insert_rate = this->n_inserts / sum;
    double delete_rate = this->n_deletes / sum;
    double indel_rate = double(this->n_inserts + this->n_deletes) / sum;

    s += "mismatch_rate:\t" + std::to_string(mismatch_rate) + "\n";
    s += "insert_rate:\t" + std::to_string(insert_rate) + "\n";
    s += "delete_rate:\t" + std::to_string(delete_rate) + "\n";
    s += "indel_rate:\t" + std::to_string(indel_rate) + "\n";
    s += "mismatches per 100kb:\t" + std::to_string(mismatch_rate*100000) + "\n";
    s += "inserts per 100kb:\t" + std::to_string(insert_rate*100000) + "\n";
    s += "deletes per 100kb:\t" + std::to_string(delete_rate*100000) + "\n";
    s += "indels per 100kb:\t" + std::to_string(indel_rate*100000) + "\n";

    if (verbose) {
        for (auto&[cigar_code, lengths]: this->cigar_lengths) {
            s += ">\"" + Cigar::cigar_name_key.at(cigar_code) + "\"\n";
            for (auto&[length, count]: lengths) {
                s += std::to_string(length) + '\t' + std::to_string(count) + '\n';
            }
        }
    }

    return s;
}


// Helper function for iterating BAMs
void chunk_sequence(vector<Region>& regions, string read_name, uint64_t chunk_size, uint64_t length){
    uint64_t l = 0;
    uint64_t start;
    uint64_t stop;

    // Chunk the sequence by "chunk_size"
    while (l < length - 1){
        start = l;
        stop = l + chunk_size - 1;

        // If the current chunk would exceed the length of the sequence
        if (stop >= length){
            // Set the stop point to the last index. A chunk can be 0-width.. e.g.: [9,9]
            stop = length - 1;
        }

        regions.emplace_back(read_name, start, stop);
        l += chunk_size;
    }
}


BamReader::BamReader() = default;


BamReader::~BamReader() {
//    hts_close(this->bam_file);
//    bam_hdr_destroy(this->bam_header);
//    bam_destroy1(this->alignment);
//    free(this->bam_index);
}


BamReader::BamReader(path bam_path){
    this->bam_path = bam_path.string();
    this->bam_file = nullptr;
    this->bam_index = nullptr;
    this->bam_iterator = nullptr;
    this->alignment = bam_init1();

//    this->secondary_mask = 256;
//    this->supplementary_mask = 2048;

    this->valid_region = false;
    this->ref_name = "";
    this->region_start = -1;
    this->region_stop = -1;

    // bam file
    if ((this->bam_file = hts_open(this->bam_path.string().c_str(), "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file: " + string(this->bam_path));
    }

    // bam index
    if ((this->bam_index = sam_index_load(this->bam_file, this->bam_path.string().c_str())) == 0) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + string(this->bam_path) + "\n");
    }

    // bam header
    if ((this->bam_header = sam_hdr_read(this->bam_file)) == 0){
        throw runtime_error("ERROR: Cannot open header for bam file: " + string(this->bam_path) + "\n");
    }
}


void BamReader::initialize_region(string& reference_name, uint64_t start, uint64_t stop){
    this->ref_name = reference_name;

    // Find the ID for this contig/chromosome/region
    this->ref_id = bam_name2id(this->bam_header, reference_name.c_str());

    // sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
    this->region_start = start;
    this->region_stop = stop;
    this->bam_iterator = sam_itr_queryi(this->bam_index, this->ref_id, start, stop);

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

    aligned_segment = {};

    aligned_segment.ref_start_index = alignment->core.pos + 1;
    aligned_segment.ref_name = bam_header->target_name[alignment->core.tid];
    aligned_segment.read_length = alignment->core.l_qseq;
    aligned_segment.read_sequence = bam_get_seq(alignment);
    aligned_segment.read_name = bam_get_qname(alignment);
    aligned_segment.cigars = bam_get_cigar(alignment);
    aligned_segment.n_cigar = alignment->core.n_cigar;
    aligned_segment.reversal = bam_is_rev(alignment);
    aligned_segment.is_secondary = ((alignment->core.flag & BamReader::secondary_mask) == 0);
    aligned_segment.is_supplementary = ((alignment->core.flag & BamReader::supplementary_mask) == 0);
    aligned_segment.map_quality = alignment->core.qual;

    aligned_segment.initialize_cigar_iterator();
}


bool BamReader::next_alignment(AlignedSegment& aligned_segment,
        uint16_t map_quality_cutoff,
        bool filter_secondary,
        bool filter_supplementary){
    ///
    /// Iterate the alignments within a region, returning references
    ///

    if (not this->valid_region) {
        throw runtime_error("ERROR: BAM reader must be initialized with `initialize_region()` before iterating");
    }

    bool found_valid_alignment = false;
    while ((not found_valid_alignment) and this->valid_region) {
        // Free alignment member pointers
        bam_destroy1(this->alignment);
        this->alignment = bam_init1();

        // Call next() on samtools
        int64_t result;
        if ((result = sam_itr_next(this->bam_file, this->bam_iterator, this->alignment)) >= 0) {

            // Load alignment into container
            load_alignment(aligned_segment, this->alignment, this->bam_header);
            found_valid_alignment = true;

            // Check secondary filter
            if (filter_secondary and (not aligned_segment.is_secondary)){
                found_valid_alignment = false;
            }

            // Check secondary filter
            if (filter_supplementary and (not aligned_segment.is_supplementary)){
                found_valid_alignment = false;
            }

            // Check map quality filter
            if (uint16_t(aligned_segment.map_quality) <= map_quality_cutoff){
                found_valid_alignment = false;
            }

        } else {
            // No more alignments left
            this->valid_region = false;
        }

    }
    return this->valid_region;
}

