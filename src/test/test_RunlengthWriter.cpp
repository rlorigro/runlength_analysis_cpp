
#include "RunlengthWriter.hpp"
#include "RunlengthReader.hpp"


void write_file(path absolute_output_path){
    RunlengthWriter writer = RunlengthWriter(absolute_output_path);

    vector<string> bases = {"A","C","G","T"};
    vector<uint8_t> lengths = {1,2,3,4};
    RunlengthSequenceElement sequence;

    sequence = {};
    sequence.name = "sequence_a";

    cout << "sequence.name: " << sequence.name << '\n';
    for (size_t i=0; i<lengths.size(); i++){
       sequence.sequence += bases[i%4];
       sequence.lengths.emplace_back(lengths[i]);
       cout << int(lengths[(i)]) << '\n';
    }

    writer.write_sequence(sequence);

    sequence = {};
    sequence.name = "sequence_b";

    cout << "sequence.name: " << sequence.name << '\n';
    for (size_t i=lengths.size(); i>0; i--){
       sequence.sequence += bases[i%4];
       sequence.lengths.emplace_back(lengths[(i-1)]);
       cout << int(lengths[(i-1)]) << '\n';
    }

    writer.write_sequence(sequence);
    writer.write_indexes();

}


int main(){
    path script_path = __FILE__;
    cout << script_path;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    path relative_output_path = "/output/";
    path filename = "test_RunlengthWriter.rlq";
    path absolute_output_path = project_directory / relative_output_path / filename;

    cout << "WRITING: " << absolute_output_path << "\n";

    write_file(absolute_output_path);

    RunlengthReader reader = RunlengthReader(absolute_output_path);

    RunlengthSequenceElement runlength_sequence;
    reader.get_sequence(runlength_sequence, 0);

    cout << runlength_sequence.name << '\n';
    cout << runlength_sequence.sequence << '\n';
    for (auto& length: runlength_sequence.lengths) {
        cout << int(length) << '\n';
    }
    reader.get_sequence(runlength_sequence, 1);

    cout << runlength_sequence.name << '\n';
    cout << runlength_sequence.sequence << '\n';
    for (auto& length: runlength_sequence.lengths) {
        cout << int(length) << '\n';
    }

    return 0;
}
