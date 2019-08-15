#ifndef RUNLENGTH_ANALYSIS_READ_H
#define RUNLENGTH_ANALYSIS_READ_H

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <string>
#include <array>
#include <unordered_set>
#include <unordered_map>

using std::string;
using std::array;
using std::unordered_set;
using std::unordered_map;


class Coordinate {
public:
    /// Attributes ///
    int64_t ref_index;          // Index of full reference sequence
    int64_t read_index;         // Index of aligned (possibly clipped) SAM sequence
    int64_t read_true_index;    // Index of true unaligned fasta sequence

    Coordinate();
};


class Cigar{
public:
    /// Attributes ///
    uint64_t length;
    uint8_t code;
    static const array <string,10> cigar_name_key;
    static const unordered_map<string,uint8_t> cigar_code_key;

    /// Methods ///
    Cigar(uint32_t bytes = 0b00000011);
    Cigar(uint8_t cigar_code, uint64_t cigar_length);
    string get_cigar_code_as_string();
    bool is_not_clip();
    bool is_ref_move();
    bool is_read_move();
    string to_string();
};


// Used primarily as a container inside BamReader, which needs metadata to iterate the cigar/alignment data
class AlignedSegment{
public:
    /// Attributes ///
    int64_t ref_start_index;     // Left most position of alignment
    string ref_name;             // Reference contig name (usually chromosome)
    int64_t read_length;         // Length of the read.
    uint8_t* read_sequence;      // DNA sequence
    string read_name;
    uint32_t* cigars;
    uint32_t n_cigar;
    bool reversal;
    bool is_secondary;
    int64_t map_quality;

    // Class constants
    static const array <string, 2> bases;
    static const array <bool, 10> cigar_ref_move;
    static const array <bool, 10> cigar_read_move;
    static const array <bool, 10> cigar_true_read_move;

    // Each nt in the sequence is stored within a uint8 with 8 bits, XXXXYYYY, where XXXX is nt1 and YYYY is nt2
    static const uint8_t bam_sequence_shift = 4;
    static const uint8_t bam_sequence_mask = 15;       // 0b1111
    static const uint8_t bam_cigar_shift = 4;
    static const uint8_t bam_cigar_mask = 15;          // 0b1111

    /// Methods ///
    string to_string();
    void initialize_cigar_iterator();
    string get_read_base(int64_t i);
    uint64_t get_ref_index_increment(Cigar& cigar);
    uint64_t get_read_index_increment(Cigar& cigar);
    int64_t infer_reference_stop_position_from_alignment();
    bool next_cigar_in_bounds();
    bool current_cigar_in_bounds();
    bool next_cigar();
    bool next_valid_cigar(Coordinate& coordinate, Cigar& cigar, unordered_set<uint8_t>& target_cigar_codes);
    bool next_coordinate(Coordinate& coordinate, Cigar& cigar);
    bool next_coordinate(Coordinate& coordinate, Cigar& cigar, unordered_set<uint8_t>& target_cigar_codes);
    void update_containers(Coordinate& coordinate, Cigar& cigar);
    void increment_coordinate(Coordinate& coordinate, Cigar& cigar, uint64_t length=1);


private:
    /// Attributes ///

    // Volatile. For tracking positions during read/ref iteration
    int64_t ref_iterator_start_index;
    int64_t read_iterator_start_index;
    int64_t read_true_iterator_start_index;
    int64_t cigar_iterator_start;
    int8_t increment;
    int64_t ref_increment;
    int64_t read_increment;
    int64_t ref_index;
    int64_t read_index;
    int64_t read_true_index;
    int64_t cigar_index;
    int64_t subcigar_index;
    Cigar current_cigar = Cigar(3);
//    bool valid_iterator;

    /// Methods ///

};

#endif //RUNLENGTH_ANALYSIS_READ_H
