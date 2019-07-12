#include <thread>
#include <vector>
#include <string>
#include <iostream>
#include "Runlength.hpp"

using std::vector;
using std::string;
using std::thread;
using std::cout;


void runlength_encode(runlength_sequence_element& runlength_sequence, sequence_element& sequence){
    // First just copy the name
    runlength_sequence.name = sequence.name;

    char current_character = 0;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (auto& character: sequence.sequence){
        if (tolower(character) != tolower(current_character)){
            runlength_sequence.sequence += character;
            runlength_sequence.lengths.push_back(1);
        }
        else{
            runlength_sequence.lengths.back()++;
        }

        current_character = character;
    }
}
