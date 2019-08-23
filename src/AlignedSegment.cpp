#include "AlignedSegment.hpp"
#include "htslib/hts.h"
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

using std::cout;
using std::string;
using std::runtime_error;
using std::unordered_map;
using std::unordered_set;


const array <string,10> Cigar::cigar_name_key = {"M",    // 0 BAM_CMATCH
                                            "I",    // 1 BAM_CINS
                                            "D",    // 2 BAM_CDEL
                                            "N",    // 3 BAM_CREF_SKIP
                                            "S",    // 4 BAM_CSOFT_CLIP
                                            "H",    // 5 BAM_CHARD_CLIP
                                            "P",    // 6 BAM_CPAD
                                            "=",    // 7 BAM_CEQUAL
                                            "X"};   // 8 BAM_CDIF


const unordered_map<string,uint8_t> Cigar::cigar_code_key = {{"M", 0},
                                                             {"I", 1},
                                                             {"D", 2},
                                                             {"N", 3},
                                                             {"S", 4},
                                                             {"H", 5},
                                                             {"P", 6},
                                                             {"=", 7},
                                                             {"X", 8}};


const array <string, 2> AlignedSegment::bases = {"=ACMGRSVTWYHKDBN", "=TGKCYSBAWRDKHVN"};


const array <bool, 10> AlignedSegment::cigar_ref_move = {true,      //BAM_CMATCH      0
                                                          false,     //BAM_CINS        1
                                                          true,      //BAM_CDEL        2
                                                          true,      //BAM_CREF_SKIP   3
                                                          false,     //BAM_CSOFT_CLIP  4
                                                          false,     //BAM_CHARD_CLIP  5
                                                          false,     //BAM_CPAD        6
                                                          true,      //BAM_CEQUAL      7
                                                          true,      //BAM_CDIFF       8
                                                          false};    //BAM_CBACK       9


const array <bool, 10> AlignedSegment::cigar_read_move = {true,     //BAM_CMATCH      0
                                                           true,     //BAM_CINS        1
                                                           false,    //BAM_CDEL        2
                                                           false,    //BAM_CREF_SKIP   3
                                                           true,     //BAM_CSOFT_CLIP  4
                                                           false,    //BAM_CHARD_CLIP  5
                                                           false,    //BAM_CPAD        6
                                                           true,     //BAM_CEQUAL      7
                                                           true,     //BAM_CDIFF       8
                                                           false};   //BAM_CBACK       9


const array <bool, 10> AlignedSegment::cigar_true_read_move = {true,     //BAM_CMATCH      0
                                                              true,     //BAM_CINS        1
                                                              false,    //BAM_CDEL        2
                                                              false,    //BAM_CREF_SKIP   3
                                                              true,     //BAM_CSOFT_CLIP  4
                                                              true,     //BAM_CHARD_CLIP  5 - in the raw sequence no bases are lost in alignment
                                                              false,    //BAM_CPAD        6
                                                              true,     //BAM_CEQUAL      7
                                                              true,     //BAM_CDIFF       8
                                                              false};   //BAM_CBACK       9


Coordinate::Coordinate() = default;


Cigar::Cigar(uint32_t bytes){
    ///
    /// Decompress byte from BAM into the 2 components of a cigar operation: (code, length)
    ///
    this->code = bytes & AlignedSegment::bam_cigar_mask;
    this->length = bytes >> AlignedSegment::bam_cigar_shift;
}


Cigar::Cigar(uint8_t cigar_code, uint64_t cigar_length){
    this->code = cigar_code;
    this->length = cigar_length;
}


string Cigar::get_cigar_code_as_string(){
    return Cigar::cigar_name_key.at(this->code);
}


bool Cigar::is_not_clip(){
    bool valid = false;

    if ((this->code < 4) or (this->code > 6)){
        valid = true;
    }

    return valid;
}


bool Cigar::is_ref_move(){
    return AlignedSegment::cigar_ref_move[this->code];
}


bool Cigar::is_read_move(){
    return AlignedSegment::cigar_read_move[this->code];
}


bool Cigar::is_true_read_move(){
    return AlignedSegment::cigar_true_read_move[this->code];
}


string Cigar::to_string() {
    return "(" + Cigar::cigar_name_key.at(this->code) + "," + std::to_string(this->length) + ")";
}


string AlignedSegment::to_string(){
    string reversal_string;

    if (reversal) {
        reversal_string = "R";
    } else {
        reversal_string = "F";
    }

    string s = this->ref_name + " "
               + std::to_string(this->ref_start_index) + " "
               + reversal_string + " "
               + std::to_string(this->read_length) + " "
               + this->read_name + " "
               + std::to_string(this->n_cigar);

    return s;
}


uint64_t AlignedSegment::get_read_index_increment(Cigar& cigar){
    ///
    /// Determine whether a cigar operation should result in a step along the READ sequence
    ///

    uint64_t index_increment = 0;
    if (AlignedSegment::cigar_read_move[cigar.code]) {
        index_increment = cigar.length;
    }
    return index_increment;
}


uint64_t AlignedSegment::get_ref_index_increment(Cigar& cigar){
    ///
    /// Determine whether a cigar operation should result in a step along the REF sequence
    ///

    uint64_t index_increment = 0;
    if (AlignedSegment::cigar_ref_move[cigar.code]){
        index_increment = cigar.length;
    }
    return index_increment;
}


string AlignedSegment::get_read_base(int64_t i){
    ///
    /// Convert the compressed representation of an aligned sequence into a string
    ///

    uint8_t byte;
    uint8_t base_code;
    string base;

    // Fetch 4 bit base code from the correct 8 bit integer and convert to a char
    if (llabs(this->read_iterator_start_index - i) < this->read_length){
        uint64_t sequence_index = i/2;
        byte = this->read_sequence[sequence_index];

        if (i%2 == 0){
            // Perform bit shift and decode using the standard or complemented base map
            base_code = byte >> AlignedSegment::bam_sequence_shift;
            base = AlignedSegment::bases[this->reversal][base_code];
        }
        else {
            // Perform bit mask and decode using the standard or complemented base map
            base_code = byte & AlignedSegment::bam_sequence_mask;
            base = AlignedSegment::bases[this->reversal][base_code];
        }
    }
    else{
        throw runtime_error("ERROR: read index out of bounds: " + this->read_name + " " + std::to_string(i));
    }
    return base;
}


int64_t AlignedSegment::infer_reference_stop_position_from_alignment(){
    ///
    /// Iterate the cigar of an alignment and find the reference stop position
    ///
    int64_t ref_stop_position = this->ref_start_index;
    for (int64_t i=0; i<this->n_cigar; i++) {
        Cigar cigar = Cigar(this->cigars[i]);

        if (cigar.is_not_clip()) {
            ref_stop_position += get_ref_index_increment(cigar);
        }
    }

    return ref_stop_position - 1;
}


void AlignedSegment::initialize_cigar_iterator(){
    this->ref_iterator_start_index = this->ref_start_index - 1;     // SAMs are 1-based
    this->read_iterator_start_index = 0;
    this->read_true_iterator_start_index = 0;
    this->cigar_iterator_start = 0;
    this->increment = 1;

    // Index the sequence in its F direction (not alignment direction)
    if (this->reversal){
        this->ref_iterator_start_index = this->infer_reference_stop_position_from_alignment() - 1;
        this->read_iterator_start_index = this->read_length - 1;
        this->read_true_iterator_start_index = 0;
        this->cigar_iterator_start = this->n_cigar - 1;
        this->increment = -1;
    }

    // Initialize all the indexes
    this->ref_index = this->ref_iterator_start_index;
    this->read_index = this->read_iterator_start_index;
    this->read_true_index = this->read_true_iterator_start_index;
    this->cigar_index = this->cigar_iterator_start;

    this->subcigar_index = 0;
}


bool AlignedSegment::next_cigar_in_bounds(){
    return (llabs(this->cigar_iterator_start - cigar_index) < this->n_cigar);
}


bool AlignedSegment::current_cigar_in_bounds(){
    return (llabs(this->cigar_iterator_start - (cigar_index))-1 < this->n_cigar);
}


//bool AlignedSegment::no_more_cigars(){
//    return (llabs(this->cigar_iterator_start - cigar_index) >= this->n_cigar);
//}


//bool AlignedSegment::is_valid_cigar(set<uint8_t> target_cigar_codes){
//    return (target_cigar_codes.find(this->current_cigar.code) != target_cigar_codes.end());
//}


bool AlignedSegment::next_cigar(){
    ///
    /// Update the class-level iterator to contain the next cigar operation in this->cigars, if it exists
    ///

    bool valid = false;

    if (this->next_cigar_in_bounds()){
        this->current_cigar = Cigar(this->cigars[cigar_index]);

        // Increment may be negative if read is reverse
        this->cigar_index += this->increment;
        valid = true;
    }

    this->subcigar_index = 0;
    return valid;
}


void AlignedSegment::update_containers(Coordinate& coordinate, Cigar& cigar){
    coordinate.read_index = this->read_index;
    coordinate.read_true_index = this->read_true_index;
    coordinate.ref_index = this->ref_index;
    cigar.length = this->current_cigar.length;
    cigar.code = this->current_cigar.code;
}


void AlignedSegment::increment_coordinate(Coordinate& coordinate, Cigar& cigar, uint64_t length){
    ///
    /// Walk along the alignment in the direction of the read, tracking ref/read indexes
    ///
    // SAM sequence may be reversed
    this->read_index += AlignedSegment::cigar_read_move[current_cigar.code]*this->increment*length;

    // True sequence always walks forward
    this->read_true_index += AlignedSegment::cigar_true_read_move[current_cigar.code]*length;

    // Ref sequence may be reversed to match read direction
    this->ref_index += AlignedSegment::cigar_ref_move[current_cigar.code]*this->increment*length;

    // Increment sub-cigar index
    this->subcigar_index += length;
}


bool AlignedSegment::next_coordinate(Coordinate& coordinate, Cigar& cigar){
    ///
    /// Iterate the entire alignment step wise for each each l in L where L = sum(cigar.length) for all cigars
    ///
    coordinate = {};
    cigar = {};


    if (not this->current_cigar_in_bounds()){
        throw runtime_error("ERROR: AlignedSegment uninitialized or out of bounds.\n"
                            "\tAlignedSegment must be initialized with `initialize_cigar_iterator()` before iterating");
    }

    // Finished a cigar operation on last iteration
    if (this->subcigar_index == (int64_t)this->current_cigar.length) {
        // If there is another cigar, load it and increment iterators
        if (this->next_cigar()){
            this->update_containers(coordinate, cigar);
            this->increment_coordinate(coordinate, cigar);

            return this->current_cigar_in_bounds();
        }
        else{
            return false;
        }
    }

    // In the middle of a cigar operation
    else{
        this->update_containers(coordinate, cigar);
        this->increment_coordinate(coordinate, cigar);

        return this->current_cigar_in_bounds();
    }
}


bool AlignedSegment::next_valid_cigar(Coordinate& coordinate, Cigar& cigar, unordered_set<uint8_t>& target_cigar_codes){
    ///
    /// Jump forward through alignment until a valid cigar operation is found
    ///

    if (this->next_cigar()) {
        // If there is an invalid cigar, continue incrementing ref/read indexes and calling next_cigar
        while ((target_cigar_codes.count(this->current_cigar.code) == 0) and (this->cigar_index <= this->n_cigar)) {
            // Increment by the length of the cigar
            this->increment_coordinate(coordinate, cigar, cigar.length);

            // Load next cigar if it exists
            if (this->next_cigar()) {
                cigar = this->current_cigar;
            } else {
                // No more cigars to load
                return false;
            }
        }
        return this->current_cigar_in_bounds();
    }
    else{
        // No more cigars to load
        return false;
    }
}


bool AlignedSegment::next_coordinate(Coordinate& coordinate, Cigar& cigar, unordered_set<uint8_t>& target_cigar_codes){
    ///
    /// Iterate the entire alignment step wise for each each l in L where L = sum(cigar.length) for all cigars
    ///
    coordinate = {};
    cigar = {};

    if (not this->current_cigar_in_bounds()){
        throw runtime_error("ERROR: AlignedSegment uninitialized or out of bounds.\n"
                            "\tAlignedSegment must be initialized with `initialize_cigar_iterator()` before iterating");
    }

    // Finished a cigar operation on last iteration
    if (this->subcigar_index == (int64_t)this->current_cigar.length) {
        // Fetch the next valid cigar if it exists
        if (this->next_valid_cigar(coordinate, this->current_cigar, target_cigar_codes)) {
            this->update_containers(coordinate, cigar);
            this->increment_coordinate(coordinate, cigar);

            return this->current_cigar_in_bounds();
        }
        else{
            return false;
        }
    }

    // In the middle of a cigar operation
    else{
        this->update_containers(coordinate, cigar);
        this->increment_coordinate(coordinate, cigar);

        return this->current_cigar_in_bounds();
    }
}
