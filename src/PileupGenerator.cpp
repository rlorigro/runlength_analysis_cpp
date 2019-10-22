
#include "PileupGenerator.hpp"
#include <iostream>
#include <experimental/filesystem>

using std::cout;
using std::ostream;
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


void PileupGenerator::print(Pileup& pileup){
    vector<string> pileup_strings;

    cout << "Matrix of size " << pileup.pileup.size() << " " << pileup.pileup[0].size() << "\n";

    for (size_t width_index = 0; width_index<pileup.pileup.size(); width_index++){
        for (size_t depth_index = 0; depth_index < pileup.pileup[width_index].size(); depth_index++){
            if (depth_index>=pileup_strings.size()){
                pileup_strings.push_back("");
            }

            pileup_strings[depth_index] += to_string(pileup.pileup[width_index][depth_index][0]);

            // If there are inserts in this column, append them to the strings
            if (pileup.inserts.count(width_index) > 0){
                for (auto& column: pileup.inserts.at(width_index)){
                    pileup_strings[depth_index] += to_string(column[depth_index][0]);
                }
            }
        }
    }

    for (auto& s: pileup_strings){
        cout << s << "\n";
    }
}


void PileupGenerator::parse_insert(Pileup& pileup, int64_t pileup_width_index, int64_t pileup_depth_index, AlignedSegment& aligned_segment, vector<float>& read_data){
    uint64_t insert_anchor_index = pileup_width_index + aligned_segment.reversal;
    uint64_t insert_index = aligned_segment.subcigar_index - 1;

    // If there is already another insert anchored at this ref position:
    if (pileup.inserts.count(insert_anchor_index) > 0) {

        // If this insert will fit within the width of the insert columns already present
        if (size_t(insert_index) < pileup.inserts.at(insert_anchor_index).size()) {
            // Simply fill in the base
            pileup.inserts.at(insert_anchor_index)[insert_index][pileup_depth_index] = read_data;
        } else {
            // Add another column and then fill in the base
            pileup.inserts.at(insert_anchor_index).push_back(this->default_insert_column);      // copy value because value is stored as class member
            pileup.inserts.at(insert_anchor_index)[insert_index][pileup_depth_index] = read_data;
        }
    } else {
        // Add a new entry at this position, initialize with a vector, and then fill in the base
        pileup.inserts.insert({insert_anchor_index, this->default_insert_pileup});        // copy value because value is stored as class member
        pileup.inserts.at(insert_anchor_index)[insert_index][pileup_depth_index] = read_data;
    }
}


//void PileupGenerator::get_base(string& read_base, Cigar& cigar, Coordinate& coordinate, SequenceElement& read_sequence){
//    if (cigar.is_read_move()) {
//        read_base = read_sequence.sequence[coordinate.read_true_index];
//    }
//    else{
//        read_base = this->null_character;
//    }
//}


template <class T> void PileupGenerator::fetch_region(Region& region,
        FastaReader& ref_reader,
        T& sequence_reader) {

    // Initialize BAM reader and relevant containers
    bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    // Initialize reader containers
    auto read_sequence = sequence_reader.generate_sequence_container();
    Pileup pileup = sequence_reader.generate_pileup_container();
    SequenceElement ref_sequence;

    this->default_data_vector = read_sequence.generate_default_data_vector();
    this->default_insert_column = vector <vector <float> >(this->maximum_depth, this->default_data_vector);
    this->default_insert_pileup = vector <vector <vector <float> > > (1, this->default_insert_column);

    ref_reader.get_sequence(ref_sequence, region.name);

    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;
    bool in_left_bound;
    bool in_right_bound;
    int64_t pileup_width_index;
    int64_t pileup_depth_index;
    vector<float> read_data;

    this->lowest_free_index_per_depth = {{0,0}};
    size_t region_size = region.stop - region.start;

    while (bam_reader.next_alignment(aligned_segment)) {
        cout << aligned_segment.to_string() << "\n";

        sequence_reader.get_sequence(read_sequence, aligned_segment.read_name);
        pileup_depth_index = find_depth_index(aligned_segment.ref_start_index);

        // Update the occupancy status of this row
        this->lowest_free_index_per_depth[0].second = aligned_segment.infer_reference_stop_position_from_alignment() - region.start;

        while (aligned_segment.next_coordinate(coordinate, cigar) and (pileup_depth_index < maximum_depth)) {
            in_left_bound = (coordinate.ref_index >= int64_t(region.start));
            in_right_bound = (coordinate.ref_index <= int64_t(region.stop));

            if (in_left_bound and in_right_bound){
                // Update the current width index
                pileup_width_index = coordinate.ref_index - region.start;

                //TODO: Replace with method from the SequenceElement that pushes values into a vector of length n_channels
                read_sequence.get_read_data(read_data, cigar, coordinate, read_sequence);

                if (cigar.is_ref_move()) {
                    // Update the pileup base
                    pileup.pileup[pileup_width_index][pileup_depth_index] = read_data;
                }
                else{
                    // Do insert-related stuff
                    if (cigar.get_cigar_code_as_string() == "I"){
                        this->parse_insert(pileup, pileup_width_index, pileup_depth_index, aligned_segment, read_data);
                    }
                }
            }
        }
    }
    this->print(pileup);
}
