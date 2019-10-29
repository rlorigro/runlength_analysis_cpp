
#ifndef RUNLENGTH_ANALYSIS_MATRIX_HPP
#define RUNLENGTH_ANALYSIS_MATRIX_HPP

#include "boost/multi_array.hpp"
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <string>

using boost::multi_array;
using std::vector;
using std::cout;
using std::tie;
using std::runtime_error;
using std::string;
using std::to_string;


typedef multi_array<float,4> runlength_matrix;
typedef multi_array<float,2> reference_runlength_matrix;

void operator+=(runlength_matrix& matrix_a, runlength_matrix& matrix_b);

void operator+=(runlength_matrix& matrix_a, float increment);

string reference_matrix_to_string(reference_runlength_matrix& matrix, size_t cutoff=0);

string matrix_to_string(runlength_matrix matrix, size_t cutoff=0);

void increment_matrix(runlength_matrix& matrix_a, runlength_matrix& matrix_b);

void increment_matrix(runlength_matrix& matrix_a, float increment);

runlength_matrix sum_matrices(vector<runlength_matrix> matrices);

runlength_matrix sum_reverse_complements(runlength_matrix& matrix);

void update_runlength_matrix_with_weibull_probabilities(runlength_matrix& matrix,
        bool& reversal,
        uint8_t& base_index,
        uint16_t& true_length,
        float& scale,
        float& shape);



#endif //RUNLENGTH_ANALYSIS_MATRIX_HPP
