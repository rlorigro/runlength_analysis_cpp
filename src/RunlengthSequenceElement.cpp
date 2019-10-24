
#include "RunlengthSequenceElement.hpp"
#include "Base.hpp"


//SequenceElement::SequenceElement()=default;


void RunlengthSequenceElement::get_read_data(vector<float>& read_data, Cigar& cigar, Coordinate& coordinate){
    read_data = {};

    if (cigar.is_read_move()) {
        read_data.emplace_back(base_to_float(this->sequence[coordinate.read_true_index]));
        read_data.emplace_back(float(this->lengths[coordinate.read_true_index]));
    }
    else{
        read_data.emplace_back(Pileup::DELETE_CODE);
    }
}


void RunlengthSequenceElement::generate_default_data_vector(vector<float>& read_data){
    read_data = {};
    read_data.push_back(Pileup::INSERT_CODE);
    read_data.push_back(0);
}
