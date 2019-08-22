
#include "PileupGenerator.hpp"

#include "AlignedSegment.hpp"
#include "BamReader.hpp"
#include "FastaReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <cassert>
#include <Runlength.hpp>

using std::cout;
using std::ostream;
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
        // Check if there is at least 1 space between the lowest free index and the start index for this read
        if (start_index > lowest_width_index + 1){
            depth_index = this->lowest_free_index_per_depth[0].first;
            cout << "old row\n" << depth_index << " " << lowest_width_index << " " << start_index << "\n";
        }
        // If there is not, then just add another row, and set the depth index to that row
        else{
            depth_index = this->lowest_free_index_per_depth.size();
            cout << "new row\n" << depth_index;
            this->lowest_free_index_per_depth.emplace_front(this->lowest_free_index_per_depth.size(), start_index);
        }
    }
    else{
        depth_index = 0;
    }

    return depth_index;
}


void PileupGenerator::print_lowest_free_indexes(){
    for (auto& item: this->lowest_free_index_per_depth){
        cout << item.first << " " << item.second << "\n";
    }
}


void PileupGenerator::print_matrix(){
    vector<string> pileup_strings;

    cout << "Matrix of size " << this->pileup.size() << " " << this->pileup[0].size() << "\n";

    for (size_t width_index = 0; width_index<this->pileup.size(); width_index++){
        for (size_t depth_index = 0; depth_index < this->pileup[width_index].size(); depth_index++){
            if (depth_index>=pileup_strings.size()){
                pileup_strings.push_back("");
            }

            pileup_strings[depth_index] += this->pileup[width_index][depth_index];
        }
    }

    for (auto& s: pileup_strings){
        cout << s << "\n";
    }
}


void PileupGenerator::fetch_region(Region region, uint16_t maximum_depth) {
    // Initialize BAM reader and relevant containers
    bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    // Initialize FASTA reader containers
    SequenceElement read_sequence;
    SequenceElement ref_sequence;

    ref_fasta_reader.fetch_sequence(ref_sequence, region.name);

    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;
    bool in_left_bound;
    bool in_right_bound;
    int64_t pileup_width_index;
    int64_t pileup_depth_index;
    uint64_t insert_index;
    string read_base;

    this->lowest_free_index_per_depth = {{0,0}};

    size_t region_size = region.stop - region.start;

    this->pileup = vector <vector <string> >(region_size, vector <string>(maximum_depth, this->null_character));

    while (bam_reader.next_alignment(aligned_segment)) {
        reads_fasta_reader.fetch_sequence(read_sequence, aligned_segment.read_name);
        pileup_depth_index = find_depth_index(aligned_segment.ref_start_index);

        // Update the occupancy of this row
        this->lowest_free_index_per_depth[0].second = aligned_segment.infer_reference_stop_position_from_alignment() - region.start;

        while (aligned_segment.next_coordinate(coordinate, cigar) and (pileup_depth_index < maximum_depth)) {
            in_left_bound = (coordinate.ref_index >= int64_t(region.start));
            in_right_bound = (coordinate.ref_index <= int64_t(region.stop));

            if (in_left_bound and in_right_bound){
                // Update the current width index
                pileup_width_index = coordinate.ref_index - region.start;

                if (cigar.is_ref_move()) {
                    // Update the pileup base
                    read_base = read_sequence.sequence[coordinate.read_true_index];
                    this->pileup[pileup_width_index][pileup_depth_index] = read_base;
                }
                else{
                    // Do insert-related stuff
                    if (cigar.get_cigar_code_as_string() == "I"){
                        insert_index = aligned_segment.subcigar_index - 1;

                        // If there is already another insert anchored at this ref position:
                        if (this->insert_columns.count(pileup_width_index) > 0){

                            // If this insert will fit within the width of the insert columns already present
                            if (size_t(insert_index) < this->insert_columns.at(pileup_width_index).size()){
                                // Simply fill in the base
                                this->insert_columns.at(pileup_width_index)[insert_index][pileup_depth_index] = read_base;
                            }
                            else{
                                // Add another column and then fill in the base
                                this->insert_columns.at(pileup_width_index).emplace_back(vector<string>(maximum_depth,"_"));
                                this->insert_columns.at(pileup_width_index)[insert_index][pileup_depth_index] = read_base;
                            }
                        }
                        else{
                            // Add a new entry at this position and then fill in the base
                            this->insert_columns.emplace(pileup_width_index, vector <vector <string> >(1, vector<string>(maximum_depth,"_")));
                            this->insert_columns.at(pileup_width_index)[insert_index][pileup_depth_index] = read_base;
                        }
                    }
                }
            }
        }
        this->print_matrix();
    }
}
