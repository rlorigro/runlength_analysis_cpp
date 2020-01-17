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
using std::cerr;
using std::ostream;
using std::deque;
using std::tuple;
using std::unordered_map;
using std::experimental::filesystem::path;


class PileupGenerator{
public:
    /// Attributes ///
    path bam_path;
    BamReader bam_reader;
    uint16_t maximum_depth;

    /// Methods ///
    PileupGenerator(path bam_path, uint16_t maximum_depth=80);
    void print_lowest_free_indexes();

    static void to_strings(vector<vector<string>>& pileup_strings_per_channel,
            Pileup& pileup,
            size_t min_index,
            size_t max_index);

    static void print(Pileup& pileup, size_t min_index=0, size_t max_index=0);

    static void extract_runlength_sequences(vector<RunlengthSequenceElement>& pileup_sequences,
            Pileup& pileup,
            size_t min_index,
            size_t max_index);

    template <class T> void fetch_region(Region& region, T& sequence_reader, Pileup& pileup);
    void fetch_sequence_indexes_from_region(Region& region, vector <tuple <string,int64_t,int64_t> >& read_indexes);
    int64_t find_depth_index(int64_t start_index);
    void parse_insert(Pileup& pileup, int64_t pileup_width_index, int64_t pileup_depth_index, uint64_t cigar_length, AlignedSegment& aligned_segment, vector<float>& read_data);
    void update_insert_column(
        Pileup& pileup,
        int64_t pileup_depth_index,
        uint64_t insert_anchor_index,
        uint64_t insert_index,
        vector<float>& read_data);

    template <class T> void generate_reference_pileup(Pileup& pileup, Pileup& ref_pileup, Region& region, T& ref_reader);

private:
    /// Attributes ///
    deque <pair <int64_t, int64_t> > lowest_free_index_per_depth;
    vector <float> default_data_vector;
    vector <vector <float> > default_insert_column;
    vector <vector <vector <float> > > default_insert_pileup;

    /// Methods ///
    void backfill_insert_columns(Pileup& pileup);
};


template <class T> void PileupGenerator::fetch_region(Region& region,
        T& sequence_reader,
        Pileup& pileup) {

    // Initialize BAM reader and relevant containers
    bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    size_t region_size = region.stop - region.start + 1;

    // Initialize reader containers
    auto read_sequence = sequence_reader.generate_sequence_container();

    read_sequence.generate_default_data_vector(this->default_data_vector);
    this->default_insert_column = vector <vector <float> >(this->maximum_depth, this->default_data_vector);
    this->default_insert_pileup = vector <vector <vector <float> > > (1, this->default_insert_column);

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
    pileup = Pileup(read_sequence.n_channels+2, region_size, this->maximum_depth, this->default_data_vector);

    while (bam_reader.next_alignment(aligned_segment)) {
        pileup.n_alignments++;
//        cerr << "\33[2K\rParsed: "<< aligned_segment.to_string() << "\n";

        sequence_reader.get_sequence(read_sequence, aligned_segment.read_name);
        pileup_depth_index = find_depth_index(aligned_segment.ref_start_index - region.start);

        // Update the occupancy status of this row
        this->lowest_free_index_per_depth[0].second = aligned_segment.infer_reference_stop_position_from_alignment() - region.start;

        while (aligned_segment.next_coordinate(coordinate, cigar) and (pileup_depth_index < maximum_depth)) {
            in_left_bound = (coordinate.ref_index >= int64_t(region.start));
            in_right_bound = (coordinate.ref_index <= int64_t(region.stop));

            if (in_left_bound and in_right_bound){

                // Update the current width index
                pileup_width_index = coordinate.ref_index - region.start;

                read_sequence.get_read_data(read_data, cigar, coordinate, aligned_segment);

                if (cigar.is_ref_move()) {
                    // Update the pileup base
                    pileup.pileup[pileup_width_index][pileup_depth_index] = read_data;

                    // For convenience, update the max_observed_depth variable in the pileup object
                    if (size_t(pileup_depth_index + 1) > pileup.max_observed_depth){
                        pileup.max_observed_depth = pileup_depth_index + 1;
                    }

                    // Also track the coverage at every position in pileup
                    pileup.coverage_per_position[pileup_width_index]++;
                }
                else{
                    // Do insert-related stuff
                    if (cigar.get_cigar_code_as_string() == "I"){
                        this->parse_insert(pileup, pileup_width_index, pileup_depth_index, cigar.length, aligned_segment, read_data);
                    }
                }
            }
        }
    }

    this->backfill_insert_columns(pileup);
}


template <class T> void PileupGenerator::generate_reference_pileup(Pileup& pileup, Pileup& ref_pileup, Region& region, T& ref_reader){
    ///
    /// Separately generate a "pileup" object which contains the reference sequence and placeholders wherever there were
    /// inserts in the (completed) read pileup
    ///
    vector<float> data;
    int64_t ref_index;
    size_t depth = 1;

    auto ref_sequence = ref_reader.generate_sequence_container();
    ref_reader.get_sequence(ref_sequence, region.name);

    ref_sequence.generate_default_data_vector(this->default_data_vector);
    this->default_insert_column = vector <vector <float> >(depth, this->default_data_vector);
    this->default_insert_pileup = vector <vector <vector <float> > > (1, this->default_insert_column);

    size_t region_size = region.stop - region.start + 1;
    ref_pileup = Pileup(ref_sequence.n_channels+2, region_size, depth, this->default_data_vector);

    // For every position in the pileup, place a single pileup element in the reference pileup
    for (size_t width_index = 0; width_index<pileup.pileup.size()-1; width_index++){
        ref_index = region.start + width_index;
        ref_sequence.get_ref_data(data, ref_index);
        ref_pileup.pileup[width_index][depth-1] = data;

        // If there is an insert in the read pileup, add a placeholder insert in the ref pileup
        if (pileup.inserts.count(width_index) > 0) {
            ref_sequence.generate_default_data_vector(data);
            data[Pileup::BASE] = Pileup::INSERT_CODE;

            // Add an insert element to the ref for every insert column in the read pileup
            for (uint64_t i=0; i<pileup.inserts.at(width_index).size(); i++) {
                this->update_insert_column(ref_pileup,
                        depth-1,
                        width_index,
                        i,
                        data);
            }
        }
    }

    ref_pileup.max_observed_depth = 1;
}


#endif //RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
