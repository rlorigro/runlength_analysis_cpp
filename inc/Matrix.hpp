
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


typedef multi_array<double,4> runlength_matrix;

void operator+=(runlength_matrix& matrix_a, runlength_matrix& matrix_b);

string matrix_to_string(runlength_matrix matrix, size_t cutoff=0);

void increment_matrix(runlength_matrix& matrix_a, runlength_matrix& matrix_b);

runlength_matrix sum_matrices(vector<runlength_matrix> matrices);

runlength_matrix sum_reverse_complements(runlength_matrix& matrix);

#endif //RUNLENGTH_ANALYSIS_MATRIX_HPP
