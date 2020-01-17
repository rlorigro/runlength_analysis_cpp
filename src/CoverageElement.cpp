
#include "CoverageElement.hpp"

CoverageElement::CoverageElement(char base, uint16_t length, bool reversal, float weight){
    ///
    /// The lowest level object which represents one read's worth of observed alignment data at one position
    ///
    this->base = base;
    this->length = length;
    this->reversal = reversal;
    this->weight = weight;
}


const string CoverageElement::reversal_string = "FR";
const string CoverageElement::reversal_string_plus_minus = "+-";


string CoverageElement::to_string(){
    return this->base + std::to_string(this->length) + this->reversal_string[this->reversal] + std::to_string(this->weight);
}


uint8_t CoverageElement::get_base_index() {
    return base_to_index(this->base);
}


bool CoverageElement::is_conventional_base() {
    return is_valid_base(this->base);
}
