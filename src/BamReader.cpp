#include "BamReader.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <experimental/filesystem>

using std::string;
using std::to_string;
using std::vector;
using std::array;
using std::runtime_error;
using std::abs;
using std::free;
using std::experimental::filesystem::path;


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


vector<Region> chunk_sequences_from_fasta_index_into_regions(unordered_map<string,FastaIndex> index_map, uint64_t chunk_size){
    ///
    /// Take all the sequences in some iterable object and chunk their lengths
    ///
    vector <Region> regions;

    // For every sequence
    for (auto& [read_name, fasta_index]: index_map){
        chunk_sequence(regions, read_name, chunk_size, fasta_index.length);
    }

    return regions;
}


Region::Region(string name, uint64_t start, uint64_t stop){
    this->name = name;
    this->start = start;
    this->stop = stop;
}

Region::Region() = default;

string Region::to_string(){
    string s = this->name + ":" + std::to_string(this->start) + "-" + std::to_string(this->stop);
    return s;
}

BamReader::BamReader() = default;

BamReader::~BamReader() {
    free(this->bam_file);
    free(this->bam_index);
    free(this->bam_iterator);
    free(this->alignment);
}


BamReader::BamReader(path bam_path){
    this->bam_path = bam_path.string();
    this->bam_file = nullptr;
    this->bam_index = nullptr;
    this->bam_iterator = nullptr;
    this->alignment = bam_init1();

    this->secondary_mask = 256;

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

    aligned_segment = {};

    aligned_segment.ref_start_index = alignment->core.pos + 1;
    aligned_segment.ref_name = bam_header->target_name[alignment->core.tid];
    aligned_segment.read_length = alignment->core.l_qseq;
    aligned_segment.read_sequence = bam_get_seq(alignment);
    aligned_segment.read_name = bam_get_qname(alignment);
    aligned_segment.cigars = bam_get_cigar(alignment);
    aligned_segment.n_cigar = alignment->core.n_cigar;
    aligned_segment.reversal = bam_is_rev(alignment);
    aligned_segment.is_secondary = ((alignment->core.flag & this->secondary_mask) == 0);
    aligned_segment.map_quality = alignment->core.qual;

    aligned_segment.initialize_cigar_iterator();
}


bool BamReader::next_alignment(AlignedSegment& aligned_segment, uint16_t map_quality_cutoff, bool filter_secondary){
    ///
    /// Iterate the alignments within a region, returning references
    ///

    if (not this->valid_region) {
        throw runtime_error("ERROR: BAM reader must be initialized with `initialize_region()` before iterating");
    }

    bool found_valid_alignment = false;
    while ((not found_valid_alignment) and this->valid_region) {

        // Call next() on samtools
        int64_t result;
        if ((result = sam_itr_next(this->bam_file, this->bam_iterator, alignment)) >= 0) {

            // Load alignment into container
            load_alignment(aligned_segment, this->alignment, this->bam_header);
            found_valid_alignment = true;

            // Check secondary filter
            if (filter_secondary and (not aligned_segment.is_secondary)){
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

