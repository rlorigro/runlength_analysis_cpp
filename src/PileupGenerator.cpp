
#include "PileupGenerator.hpp"

#include "AlignedSegment.hpp"
#include "BamReader.hpp"
#include "FastaReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <cassert>
#include <Runlength.hpp>

using std::cout;
using std::experimental::filesystem::path;


PileupGenerator::PileupGenerator(path bam_path, path ref_fasta_path, path reads_fasta_path){
    this->bam_path = bam_path;
    this->ref_fasta_path = ref_fasta_path;
    this->reads_fasta_path = reads_fasta_path;

    this->bam_reader = BamReader(bam_path);
    this->ref_fasta_reader = FastaReader(ref_fasta_path);
    this->reads_fasta_reader = FastaReader(reads_fasta_path);
}


int64_t PileupGenerator::find_depth_index(int64_t start_index){
    ///
    /// Decide where to insert a read in the pileup, depending on what space is available.
    ///
    int64_t depth_index = -1;
    int64_t lowest_width_index;

    // Sort a vector of pairs
    sort(this->lowest_free_index_per_depth.begin(), this->lowest_free_index_per_depth.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });

    lowest_width_index = this->lowest_free_index_per_depth[0].second;
    // If this row is not empty
    if (lowest_width_index != 0){
        // Check if there is at least a space between the lowest free index and the start index for this read
        if (start_index > lowest_width_index + 1){
            depth_index = lowest_width_index;
        }
        // If there is not, then just add another row
        else{
            this->lowest_free_index_per_depth.emplace_front(this->lowest_free_index_per_depth.size(), start_index);
        }
    }

    return depth_index;
}


void PileupGenerator::fetch_region(Region region) {
    // Initialize BAM reader and relevant containers
    bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    // Initialize FASTA reader containers
    SequenceElement read_sequence;
    SequenceElement ref_sequence;

    ref_fasta_reader.fetch_sequence(ref_sequence, region.name);

    cout << cigar.to_string() << "\n";

    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;
    bool found_valid_cigar;
    bool in_left_bound;
    bool in_right_bound;
    int64_t pileup_width_index;
    int64_t pileup_depth_index;
    string read_base;

    this->lowest_free_index = {{0,0}};

    while (bam_reader.next_alignment(aligned_segment)) {
        reads_fasta_reader.fetch_sequence(read_sequence, aligned_segment.read_name);
        found_valid_cigar = false;

        cout << aligned_segment.to_string() << "\n";

        while (aligned_segment.next_coordinate(coordinate, cigar)) {
            in_left_bound = (coordinate.ref_index >= region.start);
            in_right_bound = (coordinate.ref_index <= region.stop);

            if (in_left_bound and in_right_bound){
                // Update the current width index
                pileup_width_index = coordinate.ref_index - region.start;

                // If this read is new, find the depth it should be inserted at
                if (not found_valid_cigar){
                    pileup_depth_index = find_depth_index(pileup_width_index);
                    found_valid_cigar = true;
                }

                read_base = read_sequence.sequence[coordinate.read_true_index];

                pileup[]
            }
        }
    }
}
