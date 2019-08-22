#include "DiscreteWeibull.hpp"
#include "Matrix.hpp"
#include "Base.hpp"
#include <algorithm>

using std::min;


string matrix_to_string(runlength_matrix matrix, size_t cutoff){
    string matrix_string;
    string reversal_string;

    auto shape_a = matrix.shape();

    size_t n_directions = shape_a[0];
    size_t n_bases = shape_a[1];
    size_t n_true_lengths = shape_a[2];
    size_t n_observed_lengths = shape_a[3];

    if (cutoff > 0){
        n_true_lengths = cutoff;
        n_observed_lengths = cutoff;
    }

    for (size_t r=0; r<n_directions; r++) {
        if (r == 0){
            reversal_string = "F";
        }
        else{
            reversal_string = "R";
        }

        for (size_t b=0; b<n_bases; b++) {
            matrix_string += ">" + index_to_base(uint8_t(b)) + "_" + reversal_string + "\n";

            for (size_t y=0; y<n_true_lengths; y++) {
                for (size_t x=0; x<n_observed_lengths; x++) {
                    matrix_string += to_string(matrix[r][b][y][x]);

                    if (x < n_observed_lengths-1){
                        matrix_string += ",";
                    }
                }
                matrix_string += "\n";
            }
            matrix_string += "\n";
        }
    }

    return matrix_string;
}


void operator+=(runlength_matrix& matrix_a, runlength_matrix& matrix_b){
    ///
    /// Increment 'matrix_a' element-wise with values from 'matrix_b'
    ///
    increment_matrix(matrix_a, matrix_b);
}


void increment_matrix(runlength_matrix& matrix_a, runlength_matrix& matrix_b){
    ///
    /// Increment 'matrix_a' element-wise with values from 'matrix_b'
    ///
    auto shape_a = matrix_a.shape();
    auto shape_b = matrix_b.shape();

    for (size_t i=0; i<4; i++){
        if (shape_a[i] != shape_b[i]){
            string shape_a_string = to_string(shape_a[0]) + "," + to_string(shape_a[1]) + "," + to_string(shape_a[2]) + "," + to_string(shape_a[3]);
            string shape_b_string = to_string(shape_b[0]) + "," + to_string(shape_b[1]) + "," + to_string(shape_b[2]) + "," + to_string(shape_b[3]);

            throw runtime_error("ERROR: matrices with unequal sizes cannot be added: " + shape_a_string + " " + shape_b_string);
        }
    }

    size_t n_directions = shape_a[0];
    size_t n_bases = shape_a[1];
    size_t n_true_lengths = shape_a[2];
    size_t n_observed_lengths = shape_a[3];

    for (size_t r=0; r<n_directions; r++) {
        for (size_t b=0; b<n_bases; b++) {
            for (size_t y=0; y<n_true_lengths; y++) {
                for (size_t x=0; x<n_observed_lengths; x++) {
                    matrix_a[r][b][y][x] += matrix_b[r][b][y][x];
                }
            }
        }
    }
}


void operator+=(runlength_matrix& matrix_a, double increment){
    ///
    /// Increment 'matrix_a' element-wise with fixed value
    ///
    increment_matrix(matrix_a, increment);
}


void increment_matrix(runlength_matrix& matrix_a, double increment){
    ///
    /// Increment 'matrix_a' element-wise with fixed value
    ///
    auto shape_a = matrix_a.shape();

    size_t n_directions = shape_a[0];
    size_t n_bases = shape_a[1];
    size_t n_true_lengths = shape_a[2];
    size_t n_observed_lengths = shape_a[3];

    for (size_t r=0; r<n_directions; r++) {
        for (size_t b=0; b<n_bases; b++) {
            for (size_t y=0; y<n_true_lengths; y++) {
                for (size_t x=0; x<n_observed_lengths; x++) {
                    matrix_a[r][b][y][x] += increment;
                }
            }
        }
    }
}


runlength_matrix sum_matrices(vector<runlength_matrix> matrices){
    auto shape = matrices[0].shape();

    size_t n_directions = shape[0];
    size_t n_bases = shape[1];
    size_t n_true_lengths = shape[2];
    size_t n_observed_lengths = shape[3];

    runlength_matrix sum(boost::extents[n_directions][n_bases][n_true_lengths][n_observed_lengths]);

    for (auto& matrix: matrices){
        sum += matrix;
    }

    return sum;
}


runlength_matrix sum_reverse_complements(runlength_matrix& matrix){
    ///
    /// Convert a bidirectional matrix of frequencies into a unidirectional one by summing reverse complements
    ///

    auto shape = matrix.shape();

    size_t n_directions = shape[0];
    size_t n_bases = shape[1];
    size_t n_true_lengths = shape[2];
    size_t n_observed_lengths = shape[3];

    runlength_matrix sum(boost::extents[1][n_bases][n_true_lengths][n_observed_lengths]);

    for (size_t r=0; r<n_directions; r++) {
        for (size_t b=0; b<n_bases; b++) {
            for (size_t y=0; y<n_true_lengths; y++) {
                for (size_t x=0; x<n_observed_lengths; x++) {
                    if (r==0) {
                        sum[0][b][y][x] += matrix[r][b][y][x];
                    }
                    else if (r==1){
                        sum[0][b][y][x] += matrix[r][3-b][y][x];
                    }
                    else{
                        throw runtime_error("ERROR: cannot sum reverse complements, directional dimension size != 2: " + to_string(n_directions));
                    }
                }
            }
        }
    }

    return sum;
}


void update_runlength_matrix_with_weibull_probabilities(runlength_matrix& matrix,
        bool& reversal,
        uint8_t& base_index,
        uint16_t& true_length,
        double& scale,
        double& shape){
    ///
    /// Evaluate the Discrete Weibull Distribution for a given scale and shape, taking an ALLOCATED (!) vector which
    /// defines the range to evaluate over
    ///

    auto matrix_shape = matrix.shape();

    size_t n_observed_lengths = matrix_shape[3];
    size_t n_true_lengths = matrix_shape[2];

    // Skip any counts that are too large
    if (true_length < n_true_lengths - 1){
        for (size_t i = 0; i < n_observed_lengths-1; i++) {
            matrix[reversal][base_index][true_length][i+1] +=
            evaluate_weibull_cdf(i, scale, shape) - evaluate_weibull_cdf(i + 1, scale, shape);
        }
    }
}

