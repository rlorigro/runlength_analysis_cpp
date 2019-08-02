#include "AlignedSegment.hpp"
#include "htslib/hts.h"
#include <string>
#include <stdexcept>

using std::string;
using std::runtime_error;


const array <string,10> Cigar::cigar_key = {"M",    // 0 BAM_CMATCH
                                            "I",    // 1 BAM_CINS
                                            "D",    // 2 BAM_CDEL
                                            "N",    // 3 BAM_CREF_SKIP
                                            "S",    // 4 BAM_CSOFT_CLIP
                                            "H",    // 5 BAM_CHARD_CLIP
                                            "P",    // 6 BAM_CPAD
                                            "=",    // 7 BAM_CEQUAL
                                            "X"};   // 8 BAM_CDIFF

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
    return Cigar::cigar_key[this->code];
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


string Cigar::to_string() {
    return "(" + Cigar::cigar_key[this->code] + "," + std::to_string(this->length) + ")";
}


string AlignedSegment::to_string(){
    string reversal_string = reversal ? "R" : "F";

    string s = this->ref_name + " "
               + std::to_string(this->ref_start_index) + " "
               + reversal_string + " "
               + std::to_string(this->read_length) + " "
               + this->read_name + " "
               + std::to_string(this->n_cigar) + "\n";

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


//void AlignedSegment::get_cigar(Cigar& cigar, int64_t i){
//    ///
//    /// Decompress byte from BAM into the 2 components of a cigar operation: (code, length)
//    ///
//
//    cigar.code = this->cigars[i] & BamReader::bam_cigar_mask;
//    cigar.length = this->cigars[i] >> BamReader::bam_cigar_shift;
//}


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

    this->valid_iterator = true;
}


bool AlignedSegment::next_cigar(){
    ///
    /// Update the class-level iterator to contain the next cigar operation in this->cigars, if it exists
    ///

    bool valid = false;

    if (llabs(this->cigar_iterator_start - cigar_index) < this->n_cigar){
        this->current_cigar = Cigar(this->cigars[cigar_index]);

        // Increment may be negative if read is reverse
        this->cigar_index += this->increment;
        valid = true;
    }

    return valid;
}


void AlignedSegment::increment_coordinate(Coordinate& coordinate, Cigar& cigar){
    coordinate.read_index = this->read_index;
    coordinate.read_true_index = this->read_true_index;
    coordinate.ref_index = this->ref_index;
    cigar.length = this->current_cigar.length;
    cigar.code = this->current_cigar.code;

    this->read_index += AlignedSegment::cigar_read_move[current_cigar.code]*this->increment;
    this->read_true_index += AlignedSegment::cigar_true_read_move[current_cigar.code];
    this->ref_index += AlignedSegment::cigar_ref_move[current_cigar.code]*this->increment;
    this->i_subcigar++;
}


bool AlignedSegment::next_coordinate(Coordinate& coordinate, Cigar& cigar){
    ///
    /// Iterate the entire alignment step wise for each each l in L where L = sum(cigar.length) for all cigars
    ///
    coordinate = {};
    cigar = {};

    if (not this->valid_iterator){
        throw runtime_error("ERROR: AlignedSegment uninitialized or out of bounds.\n"
                            "\tAlignedSegment must be initialized with `initialize_cigar_iterator()` before iterating");
    }

//    cout << cigar.to_string() << " " << this->n_cigar << " " << this->cigar_index <<  " " << this->i_subcigar << "\n";

    // Finished a cigar operation on last iteration
    if (this->i_subcigar == (int64_t)this->current_cigar.length) {
        // If there is another cigar, load it and increment iterators
        if (this->next_cigar()){
            this->i_subcigar = 0;
            this->increment_coordinate(coordinate, cigar);
        }
        // If no next cigar, return false
        else {
            this->valid_iterator = false;
        }
    }

    // In the middle of a cigar operation
    else if (this->i_subcigar < (int64_t)this->current_cigar.length){
        this->increment_coordinate(coordinate, cigar);
    }

    return this->valid_iterator;
}
