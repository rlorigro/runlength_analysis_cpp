
#include "PileupGenerator.hpp"
#include "Base.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <cmath>

using std::cout;
using std::ostream;
using std::min;
using std::experimental::filesystem::path;


PileupGenerator::PileupGenerator(path bam_path, uint16_t maximum_depth){
    this->bam_path = bam_path;
    this->bam_reader = BamReader(bam_path);
    this->maximum_depth = maximum_depth;
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
        }
        // If there is not, then just add another row, and set the depth index to that row
        else{
            depth_index = this->lowest_free_index_per_depth.size();
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


void PileupGenerator::extract_runlength_sequences(vector<RunlengthSequenceElement>& pileup_sequences, Pileup& pileup, size_t min_index, size_t max_index) {
    pileup_sequences = vector<RunlengthSequenceElement>(pileup.pileup[0].size());
    vector<float> value;
    uint16_t length;
    string base;

    for (size_t width_index = min_index; width_index<max_index; width_index++) {
        for (size_t depth_index = 0; depth_index < pileup.pileup[width_index].size(); depth_index++) {
            value = pileup.pileup[width_index][depth_index];

            // WARNING: this is not a template safe method, not all pileups have length on channel 2
            length = value[2];
            base = index_to_base(value[0]);

            if (value[2] > 0) {
                pileup_sequences[depth_index].lengths.emplace_back(length);
                pileup_sequences[depth_index].sequence += base;
            }

            if (pileup.inserts.count(width_index) > 0) {
                for (auto& column: pileup.inserts.at(width_index)) {
                    value = column[depth_index];

                    // WARNING: this is not a template safe method, not all pileups have length on channel 2
                    length = value[2];
                    base = index_to_base(value[0]);

                    if (value[2] > 0) {
                        pileup_sequences[depth_index].lengths.emplace_back(length);
                        pileup_sequences[depth_index].sequence += base;
                    }
                }
            }
        }
    }
}


void PileupGenerator::fetch_sequence_indexes_from_region(Region& region,
        vector <tuple <string,int64_t,int64_t> >& read_indexes) {

    // Initialize BAM reader and relevant containers
    this->bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;
    bool in_left_bound;
    bool in_right_bound;

    size_t read_true_index_start;
    size_t read_true_index_stop;

    while (this->bam_reader.next_alignment(aligned_segment)) {
        cerr << "\33[2K\rParsed: "<< aligned_segment.to_string() << "\n";

        int i = 0;
        while (aligned_segment.next_coordinate(coordinate, cigar)) {
            in_left_bound = (coordinate.ref_index >= int64_t(region.start));
            in_right_bound = (coordinate.ref_index <= int64_t(region.stop));

            if (in_left_bound and in_right_bound){
                if (i==0){
                    read_true_index_start = coordinate.read_true_index;
                }
                else{
                    read_true_index_stop = coordinate.read_true_index;
                }
                i++;
            }
        }

        read_indexes.emplace_back(aligned_segment.read_name, read_true_index_start, read_true_index_stop);
    }
}


void PileupGenerator::to_strings(vector<vector<string>>& pileup_strings_per_channel, Pileup& pileup, size_t min_index, size_t max_index) {
    if (max_index == 0) {
        max_index = pileup.pileup.size();
    }

    size_t i = 0;
    string s_value;
    float value;

    cout << "n_columns: " << pileup.pileup.size() << '\n';
    cout << "n_alignments: " << pileup.n_alignments << '\n';

    for (size_t width_index = min_index; width_index < max_index; width_index++) {
        for (size_t depth_index = 0; depth_index < pileup.pileup[width_index].size(); depth_index++) {

            i = 0;
            for (auto &pileup_strings: pileup_strings_per_channel) {
                if (depth_index >= pileup_strings.size()) {
                    pileup_strings.push_back("");
                }

                value = pileup.pileup[width_index][depth_index][i];
                if (i == 0) {
                    s_value = float_to_base(value);
                } else {
                    s_value = to_string(min(int(9), int(value)));
                }

                pileup_strings[depth_index] += s_value;

                // If there are inserts in this column, append them to the strings
                if (pileup.inserts.count(width_index) > 0) {
                    for (auto &column: pileup.inserts.at(width_index)) {
                        value = column[depth_index][i];
                        if (i == 0) {
                            s_value = float_to_base(value);
                        } else {
                            s_value = to_string(min(int(9), int(value)));
                        }

                        pileup_strings[depth_index] += s_value;
                    }
                }
                i++;
            }
        }
    }
}


void PileupGenerator::print(Pileup& pileup, size_t min_index, size_t max_index){
    vector<vector<string>> pileup_strings_per_channel(pileup.pileup[0][0].size());
    PileupGenerator::to_strings(pileup_strings_per_channel, pileup, min_index, max_index);

    for (auto& pileup_strings: pileup_strings_per_channel){
        for (auto& s: pileup_strings) {
            cout << s << "\n";
        }
        cout << '\n';
    }
}


void PileupGenerator::update_insert_column(
        Pileup& pileup,
        int64_t pileup_depth_index,
        uint64_t insert_anchor_index,
        uint64_t insert_index,
        vector<float>& read_data){

    // If there is already another insert anchored at this ref position:
    if (pileup.inserts.count(insert_anchor_index) > 0) {

        // If this insert will fit within the width of the insert columns already present
        while (size_t(insert_index) > pileup.inserts.at(insert_anchor_index).size() - 1) {
            // Add another column and then fill in the base
            pileup.inserts.at(insert_anchor_index).push_back(this->default_insert_column);      // copy value because value is stored as class member
        }

        // Fill in the base
        pileup.inserts.at(insert_anchor_index)[insert_index][pileup_depth_index] = read_data;

    } else {
        // Add a new entry at this position, initialize with a vector, and then fill in the base
        pileup.inserts.insert({insert_anchor_index, this->default_insert_pileup});        // copy value because value is stored as class member

        // If this insert will fit within the width of the insert columns already present
        while (size_t(insert_index) > pileup.inserts.at(insert_anchor_index).size() - 1) {
            // Add another column and then fill in the base
            pileup.inserts.at(insert_anchor_index).push_back(this->default_insert_column);      // copy value because value is stored as class member
        }

        pileup.inserts.at(insert_anchor_index)[insert_index][pileup_depth_index] = read_data;
    }

}

void PileupGenerator::parse_insert(Pileup& pileup, int64_t pileup_width_index, int64_t pileup_depth_index, uint64_t cigar_length, AlignedSegment& aligned_segment, vector<float>& read_data){
    uint64_t insert_anchor_index = pileup_width_index + aligned_segment.reversal - 1;

    uint64_t insert_index;
    if (aligned_segment.reversal) {
        insert_index = cigar_length - (aligned_segment.subcigar_index-1) - 1;
    }
    else{
        insert_index = aligned_segment.subcigar_index - 1;
    }

    this->update_insert_column(pileup,
            pileup_depth_index,
            insert_anchor_index,
            insert_index,
            read_data);
}


void PileupGenerator::backfill_insert_columns(Pileup& pileup){
    float left_base;
    float right_base;
    bool left_not_empty;
    bool right_not_empty;

    for (size_t width_index = 0; width_index<pileup.pileup.size(); width_index++) {
        if (pileup.inserts.count(width_index) > 0) {
            for (auto &column: pileup.inserts.at(width_index)) {
                for (size_t depth_index = 0; depth_index < column.size(); depth_index++) {
                    // Don't overwrite existing insert data
                    if (is_valid_base_index(column[depth_index][Pileup::BASE])){
                        continue;
                    }

                    left_base = pileup.pileup[width_index][depth_index][Pileup::BASE];
                    left_not_empty = (left_base != Pileup::EMPTY);

                    right_not_empty = true;
                    if (width_index + 1 < pileup.pileup.size()) {
                        right_base = pileup.pileup[width_index+1][depth_index][Pileup::BASE];
                        right_not_empty = (right_base != Pileup::EMPTY);
                    }

                    // Fill in the insert values with insert codes and reversal codes if it is flanked by read data
                    if (left_not_empty and right_not_empty){
                        column[depth_index][Pileup::BASE] = Pileup::INSERT_CODE;
                        column[depth_index][Pileup::REVERSAL] = pileup.pileup[width_index][depth_index][Pileup::REVERSAL];
                    }
                }
            }
        }
    }
}
