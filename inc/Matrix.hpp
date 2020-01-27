
#ifndef RUNLENGTH_ANALYSIS_MATRIX_HPP
#define RUNLENGTH_ANALYSIS_MATRIX_HPP

#include "boost/multi_array.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>

#include "DiscreteWeibull.hpp"
#include "Base.hpp"


using boost::multi_array;
using std::vector;
using std::cout;
using std::tie;
using std::runtime_error;
using std::string;
using std::to_string;


typedef multi_array<double,4> rle_length_matrix;
typedef multi_array<double,2> reference_rle_length_matrix;
typedef multi_array<double,3> rle_base_matrix;

class RLEConfusion {
public:
    rle_length_matrix length_matrix;
    rle_base_matrix base_matrix;

    RLEConfusion(uint16_t max_runlength);
};

void operator+=(rle_length_matrix& matrix_a, rle_length_matrix& matrix_b);

void operator+=(rle_length_matrix& matrix_a, float increment);

void operator+=(rle_base_matrix& matrix_a, rle_base_matrix& matrix_b);

void increment_matrix(rle_length_matrix& matrix_a, rle_length_matrix& matrix_b);

void increment_matrix(rle_length_matrix& matrix_a, float increment);

void increment_matrix(rle_base_matrix& matrix_a, rle_base_matrix& matrix_b);

rle_length_matrix sum_matrices(vector<rle_length_matrix> matrices);

rle_base_matrix sum_matrices(vector<rle_base_matrix> matrices);

RLEConfusion sum_matrices(vector<RLEConfusion> matrices);

rle_length_matrix sum_reverse_complements(rle_length_matrix& matrix);

rle_base_matrix sum_reverse_complements(rle_base_matrix& matrix);

void update_runlength_matrix_with_weibull_probabilities(rle_length_matrix& matrix,
                                                        bool& reversal,
                                                        uint8_t& base_index,
                                                        uint16_t& true_length,
                                                        float& scale,
                                                        float& shape);

string reference_matrix_to_string(reference_rle_length_matrix& matrix, size_t cutoff=0);

string matrix_to_string(rle_length_matrix matrix, size_t cutoff=0);

string matrix_to_string(rle_base_matrix matrix);


#endif //RUNLENGTH_ANALYSIS_MATRIX_HPP
