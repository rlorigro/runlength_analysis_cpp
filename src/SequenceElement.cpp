
#include "SequenceElement.hpp"
#include "Base.hpp"

//SequenceElement::SequenceElement()=default;

void SequenceElement::get_ref_data(vector<float>& ref_data, int64_t index){
    ref_data = {};
    ref_data.emplace_back(base_to_float(this->sequence[index]));
    ref_data.emplace_back(float(0));
}


void SequenceElement::get_read_data(vector<float>& read_data, Cigar& cigar, Coordinate& coordinate, AlignedSegment& alignment){
    read_data = {};

    if (cigar.is_read_move()) {
        if (size_t(coordinate.read_true_index) > this->sequence.size() || coordinate.read_true_index < 0){
            throw runtime_error("ERROR: out of bounds index (" + to_string(coordinate.read_true_index) + ") for read " + this->name + " of length " + to_string(this->sequence.size()));
        }

        float base = base_to_float(this->sequence[coordinate.read_true_index]);

        // Complement base if necessary
        if (alignment.reversal and is_valid_base_index(base)) {
            base = 3 - base;
        }

        read_data.emplace_back(base);
    }
    else{
        read_data.emplace_back(Pileup::DELETE_CODE);
    }

    //TODO: change this so it uses reversal data
    read_data.emplace_back(float(alignment.reversal));
}


void SequenceElement::generate_default_data_vector(vector<float>& read_data){
    read_data = {};
    read_data.push_back(Pileup::EMPTY);
    read_data.push_back(float(0));
}
