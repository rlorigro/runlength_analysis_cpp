#include "Matrix.hpp"

using std::min;


RLEConfusion::RLEConfusion(uint16_t max_runlength):
    length_matrix(boost::extents[2][4][max_runlength + 1][max_runlength + 1]),   // 0 length included
    base_matrix(boost::extents[2][4][4])
{}


string reference_matrix_to_string(reference_rle_length_matrix& matrix, size_t cutoff){
    string matrix_string;

    auto shape_a = matrix.shape();
    size_t n_bases = shape_a[0];
    size_t max_runlength = shape_a[1];

    if (cutoff == 0) {
        cutoff = max_runlength;
    }

    for (size_t b=0; b<n_bases; b++) {
        cout << ">" + index_to_base(b) + "\n";

        for (size_t i = 0; i < cutoff; i++) {
            matrix_string += to_string(matrix[b][i]) + ",";
        }
        matrix_string += '\n';
    }

    return matrix_string;
}


string matrix_to_string(rle_base_matrix matrix){
    string matrix_string;
    string reversal_string;

    auto shape_a = matrix.shape();

    size_t n_directions = shape_a[0];
    size_t n_true_bases = shape_a[1];
    size_t n_observed_bases = shape_a[2];

    for (size_t r=0; r<n_directions; r++) {
        if (r == 0){
            reversal_string = "F";
        }
        else{
            reversal_string = "R";
        }

        for (size_t b1=0; b1<n_true_bases; b1++) {
            matrix_string += ">" + index_to_base(uint8_t(b1)) + "_" + reversal_string + "\n";

            for (size_t b2=0; b2<n_observed_bases; b2++) {
                matrix_string += to_string(matrix[r][b1][b2]);

                if (b2 < n_observed_bases-1){
                    matrix_string += ",";
                }
                matrix_string += "\n";
            }

            // Write a newline if it is not the last position in the file
            if (not (r == n_directions - 1 and b1 == n_true_bases - 1)){
                matrix_string += "\n";
            }
        }
    }

    return matrix_string;
}


string matrix_to_string(rle_length_matrix matrix, size_t cutoff){
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

            // Write a newline if it is not the last position in the file
            if (not (r == n_directions - 1 and b == n_bases - 1)){
                matrix_string += "\n";
            }
        }
    }

    return matrix_string;
}


void operator+=(rle_length_matrix& matrix_a, rle_length_matrix& matrix_b){
    ///
    /// Increment 'matrix_a' element-wise with values from 'matrix_b'
    ///
    increment_matrix(matrix_a, matrix_b);
}


void operator+=(rle_base_matrix& matrix_a, rle_base_matrix& matrix_b){
    ///
    /// Increment 'matrix_a' element-wise with values from 'matrix_b'
    ///
    increment_matrix(matrix_a, matrix_b);
}


void increment_matrix(rle_length_matrix& matrix_a, rle_length_matrix& matrix_b){
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


void increment_matrix(rle_base_matrix& matrix_a, rle_base_matrix& matrix_b){
    ///
    /// Increment 'matrix_a' element-wise with values from 'matrix_b'
    ///
    auto shape_a = matrix_a.shape();
    auto shape_b = matrix_b.shape();

    for (size_t i=0; i<3; i++){
        if (shape_a[i] != shape_b[i]){
            string shape_a_string = to_string(shape_a[0]) + "," + to_string(shape_a[1]) + "," + to_string(shape_a[2]);
            string shape_b_string = to_string(shape_b[0]) + "," + to_string(shape_b[1]) + "," + to_string(shape_b[2]);

            throw runtime_error("ERROR: matrices with unequal sizes cannot be added: " + shape_a_string + " " + shape_b_string);
        }
    }

    size_t n_directions = shape_a[0];
    size_t n_true_bases = shape_a[1];
    size_t n_observed_bases = shape_a[2];

    for (size_t r=0; r<n_directions; r++) {
        for (size_t b1=0; b1<n_true_bases; b1++) {
            for (size_t b2=0; b2<n_observed_bases; b2++) {
                matrix_a[r][b1][b2] += matrix_b[r][b1][b2];
            }
        }
    }
}


void operator+=(rle_length_matrix& matrix_a, double increment){
    ///
    /// Increment 'matrix_a' element-wise with fixed value
    ///
    increment_matrix(matrix_a, increment);
}


void increment_matrix(rle_length_matrix& matrix_a, float increment){
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


rle_length_matrix sum_matrices(vector<rle_length_matrix> matrices){
    auto shape = matrices[0].shape();

    size_t n_directions = shape[0];
    size_t n_bases = shape[1];
    size_t n_true_lengths = shape[2];
    size_t n_observed_lengths = shape[3];

    rle_length_matrix sum(boost::extents[n_directions][n_bases][n_true_lengths][n_observed_lengths]);

    for (auto& matrix: matrices){
        sum += matrix;
    }

    return sum;
}


rle_base_matrix sum_matrices(vector<rle_base_matrix> matrices){
    auto shape = matrices[0].shape();

    size_t n_directions = shape[0];
    size_t n_true_bases = shape[1];
    size_t n_observed_bases = shape[2];

    rle_base_matrix sum(boost::extents[n_directions][n_true_bases][n_observed_bases]);

    for (auto& matrix: matrices){
        sum += matrix;
    }

    return sum;
}


RLEConfusion sum_matrices(vector<RLEConfusion> matrices){
    auto length_shape = matrices[0].length_matrix.shape();
    size_t n_directions = length_shape[0];

    size_t n_bases = length_shape[1];
    size_t n_true_lengths = length_shape[2];
    size_t n_observed_lengths = length_shape[3];

    auto base_shape = matrices[0].base_matrix.shape();
    size_t n_true_bases = base_shape[1];
    size_t n_observed_bases = base_shape[2];

    cout <<
    "n_directions: " << to_string(n_directions) << '\n' <<
    "n_bases: " << to_string(n_bases) << '\n' <<
    "n_true_lengths: " << to_string(n_true_lengths) << '\n' <<
    "n_observed_lengths: " << to_string(n_observed_lengths) << '\n' <<
    "n_true_bases: " << to_string(n_true_bases) << '\n' <<
    "n_observed_bases: " << to_string(n_observed_bases) << '\n';

    rle_length_matrix length_sum(boost::extents[n_directions][n_bases][n_true_lengths][n_observed_lengths]);
    rle_base_matrix base_sum(boost::extents[n_directions][n_true_bases][n_observed_bases]);

    for (auto& matrix: matrices){
        length_sum += matrix.length_matrix;
        base_sum += matrix.base_matrix;
    }

    RLEConfusion confusion_sum(n_true_lengths-1);
    confusion_sum.length_matrix = length_sum;
    confusion_sum.base_matrix = base_sum;

    return confusion_sum;
}


rle_base_matrix sum_reverse_complements(rle_base_matrix& matrix){
    ///
    /// Convert a bidirectional matrix of frequencies into a unidirectional one by summing reverse complements
    ///

    auto shape = matrix.shape();

    size_t n_directions = shape[0];
    size_t n_true_bases = shape[1];
    size_t n_observed_bases = shape[2];

    rle_base_matrix sum(boost::extents[1][n_true_bases][n_observed_bases]);

    for (size_t r=0; r<n_directions; r++) {
        for (size_t b1=0; b1<n_true_bases; b1++) {
            for (size_t b2=0; b2<n_observed_bases; b2++) {
                if (r==0) {
                    sum[0][b1][b2] += matrix[r][b1][b2];
                }
                else if (r==1){
                    sum[0][3-b1][3-b2] += matrix[r][3-b1][3-b2];
                }
                else{
                    throw runtime_error("ERROR: cannot sum reverse complements, directional dimension size != 2: " + to_string(n_directions));
                }
            }
        }
    }

    return sum;
}


rle_length_matrix sum_reverse_complements(rle_length_matrix& matrix){
    ///
    /// Convert a bidirectional matrix of frequencies into a unidirectional one by summing reverse complements
    ///

    auto shape = matrix.shape();

    size_t n_directions = shape[0];
    size_t n_bases = shape[1];
    size_t n_true_lengths = shape[2];
    size_t n_observed_lengths = shape[3];

    rle_length_matrix sum(boost::extents[1][n_bases][n_true_lengths][n_observed_lengths]);

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


void update_runlength_matrix_with_weibull_probabilities(rle_length_matrix& matrix,
                                                        bool& reversal,
                                                        uint8_t& base_index,
                                                        uint16_t& true_length,
                                                        float& scale,
                                                        float& shape){
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

