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
    void print(Pileup& pileup);

    template <class T1, class T2> void fetch_region(Region& region, T1& ref_fasta_reader, T2& sequence_reader, Pileup& pileup);
    int64_t find_depth_index(int64_t start_index);
    void parse_insert(Pileup& pileup, int64_t pileup_width_index, int64_t pileup_depth_index, AlignedSegment& aligned_segment, vector<float>& read_data);

private:
    /// Attributes ///
    deque <pair <int64_t, int64_t> > lowest_free_index_per_depth;
    vector <float> default_data_vector;
    vector <vector <float> > default_insert_column;
    vector <vector <vector <float> > > default_insert_pileup;

    /// Methods ///
    void backfill_insert_columns(Pileup& pileup);
};


template <class T1, class T2> void PileupGenerator::fetch_region(Region& region,
        T1& ref_reader,
        T2& sequence_reader,
        Pileup& pileup) {

    // Initialize BAM reader and relevant containers
    bam_reader.initialize_region(region.name, region.start, region.stop);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    size_t region_size = region.stop - region.start + 1;

    // Initialize reader containers
    auto read_sequence = sequence_reader.generate_sequence_container();
    auto ref_sequence = sequence_reader.generate_sequence_container();

    read_sequence.generate_default_data_vector(this->default_data_vector);
    this->default_insert_column = vector <vector <float> >(this->maximum_depth, this->default_data_vector);
    this->default_insert_pileup = vector <vector <vector <float> > > (1, this->default_insert_column);

//    ref_reader.get_sequence(ref_sequence, region.name);

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

//    cout << "size: " << pileup.pileup.size() << " " << pileup.pileup[0].size() << '\n';

    while (bam_reader.next_alignment(aligned_segment)) {
        pileup.n_alignments++;
        cerr << "\33[2K\rParsed: "<< aligned_segment.to_string() << "\n";

        sequence_reader.get_sequence(read_sequence, aligned_segment.read_name);
        pileup_depth_index = find_depth_index(aligned_segment.ref_start_index - region.start);

//        cout << "---" << aligned_segment.ref_start_index << " " << aligned_segment.ref_start_index - region.start << " " << pileup_depth_index << '\n';
//        print_lowest_free_indexes();

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

    this->backfill_insert_columns(pileup);
}

#endif //RUNLENGTH_ANALYSIS_PILEUPGENERATOR_HPP
