#ifndef RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
#define RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP

#include "FastaReader.hpp"
#include "BamReader.hpp"
#include "Runlength.hpp"
#include "Pileup.hpp"
#include <utility>
#include <iostream>
#include <cassert>
#include <deque>
#include <experimental/filesystem>

using std::cout;
using std::ostream;
using std::deque;
using std::unordered_map;
using std::experimental::filesystem::path;


class PileupGenerator{
public:
    /// Attributes ///
    path bam_path;
    BamReader bam_reader;
    uint16_t maximum_depth;

    string null_character = "_";
    float null_value = -1;

    /// Methods ///
    PileupGenerator(path bam_path, uint16_t maximum_depth=80);
    void print_lowest_free_indexes();
    void print(Pileup& pileup);

    template <class T> void fetch_region(Region& region, FastaReader& ref_fasta_reader, T& sequence_reader);
    int64_t find_depth_index(int64_t start_index);
    void ensure_pileup_column_index_not_empty(Pileup& pileup, int64_t depth_index, int64_t width_index);
    void parse_insert(Pileup& pileup, int64_t pileup_width_index, int64_t pileup_depth_index, AlignedSegment& aligned_segment, vector<float>& read_data);

private:
    deque <pair <int64_t, int64_t> > lowest_free_index_per_depth;
    vector <float> default_data_vector;
    vector <vector <float> > default_insert_column;
    vector <vector <vector <float> > > default_insert_pileup;
};


template <class T> void PileupGenerator::fetch_region(Region& region,
        FastaReader& ref_reader,
        T& sequence_reader) {

    // Initialize BAM reader and relevant containers
    bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    size_t region_size = region.stop - region.start;

    // Initialize reader containers
    auto read_sequence = sequence_reader.generate_sequence_container();
    SequenceElement ref_sequence;

    read_sequence.generate_default_data_vector(this->default_data_vector);
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
    Pileup pileup = Pileup(read_sequence.n_channels, region_size, this->maximum_depth, this->default_data_vector);

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
                read_sequence.get_read_data(read_data, cigar, coordinate);

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

#endif //RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
