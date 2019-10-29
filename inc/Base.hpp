#ifndef RUNLENGTH_ANALYSIS_BASE_HPP
#define RUNLENGTH_ANALYSIS_BASE_HPP

#include <string>
#include <vector>
#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;
using std::runtime_error;

uint8_t base_to_index(string& base);

float base_to_float(char base);

uint8_t base_to_index(char base);

static const vector index_to_base_map = {"A","C","G","T","*","_"};

string index_to_base(uint8_t index);

string float_to_base(float index);

bool is_valid_base(string base);

bool is_valid_base(char base);

bool is_gap(uint8_t index);

#endif //RUNLENGTH_ANALYSIS_BASE_HPP
