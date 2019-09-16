
#ifndef RUNLENGTH_ANALYSIS_REFERENCERUNLENGTH_HPP
#define RUNLENGTH_ANALYSIS_REFERENCERUNLENGTH_HPP

#include "Runlength.hpp"
#include "Matrix.hpp"
#include "Base.hpp"
#include <experimental/filesystem>
#include <cmath>

using std::min;
using std::experimental::filesystem::path;


void measure_runlength_priors_from_reference(path fasta_path, uint16_t max_runlength);


template<class T> void count_runlengths(reference_runlength_matrix& runlength_frequencies, T& sequence){
    if (runlength_frequencies.empty()){
        throw runtime_error("ERROR: uninitialized vector passed to 'evaluate_discrete_weibull', vector must be initialized");
    }

    auto shape_a = runlength_frequencies.shape();

    size_t max_runlength = shape_a[1];

    char current_character = 0;
    double current_length = 1;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (char& character: sequence.sequence){
        if (tolower(character) != tolower(current_character)){
            if (not is_valid_base(character)){
                continue;
            }

            uint8_t base_index = base_to_index(character);
            runlength_frequencies[base_index][min(current_length, double(max_runlength-1))]++;
            current_length = 1;
        }
        else{
            current_length++;
        }

        current_character = character;
    }
}


#endif //RUNLENGTH_ANALYSIS_REFERENCERUNLENGTH_HPP
