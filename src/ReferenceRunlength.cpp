
#include "ReferenceRunlength.hpp"
#include "FastaReader.hpp"
#include "Runlength.hpp"
#include "Matrix.hpp"
#include <iostream>

using std::cout;

void measure_runlength_priors_from_reference(path fasta_path, uint16_t max_runlength){
    FastaReader reader = FastaReader(fasta_path);
    reader.index();

    reference_rle_length_matrix runlength_frequencies(boost::extents[4][max_runlength + 1]);   // 0 length included

    SequenceElement sequence;
    RunlengthSequenceElement runlength_sequence;
    while (!reader.end_of_file) {
        reader.next_element(sequence);
        cout << sequence.name << '\n';

        count_runlengths(runlength_frequencies, sequence);
    }


    reference_matrix_to_string(runlength_frequencies);


}
