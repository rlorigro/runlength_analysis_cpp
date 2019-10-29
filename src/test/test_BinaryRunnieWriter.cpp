
#include "BinaryRunnieWriter.hpp"
#include "BinaryRunnieReader.hpp"


void write_file(path absolute_output_path){
    BinaryRunnieWriter writer = BinaryRunnieWriter(absolute_output_path);

    vector<string> bases = {"A","C","G","T"};
    vector<uint8_t> lengths = {1,2,3,4};
    RunnieSequenceElement sequence;

    sequence = {};
    sequence.name = "sequence_a";

    cout << "sequence.name: " << sequence.name << '\n';
    for (size_t i=0; i<lengths.size(); i++){
       sequence.sequence += bases[i%4];
       sequence.scales.emplace_back(lengths[i%4]);
       sequence.shapes.emplace_back(lengths[(i+1)%4]);
       cout << int(lengths[(i)]) << '\n';
    }

    writer.write_sequence(sequence);

    sequence = {};
    sequence.name = "sequence_b";

    cout << "sequence.name: " << sequence.name << '\n';
    for (size_t i=lengths.size(); i>0; i--){
       sequence.sequence += bases[i%4];
       sequence.scales.emplace_back(lengths[(i-1)%4]);
       sequence.shapes.emplace_back(lengths[(i)%4]);
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

    BinaryRunnieReader reader = BinaryRunnieReader(absolute_output_path);

    RunnieSequenceElement runlength_sequence;
    reader.get_sequence(runlength_sequence, 0);

    cout << runlength_sequence.name << '\n';
    cout << runlength_sequence.sequence << '\n';
    for (size_t i=0; i<runlength_sequence.scales.size(); i++) {
        cout << int(runlength_sequence.scales[i]) << " " << int(runlength_sequence.shapes[i]) << '\n';
    }
    reader.get_sequence(runlength_sequence, 1);

    cout << runlength_sequence.name << '\n';
    cout << runlength_sequence.sequence << '\n';
    for (size_t i=0; i<runlength_sequence.scales.size(); i++) {
        cout << int(runlength_sequence.scales[i]) << " " << int(runlength_sequence.shapes[i]) << '\n';
    }

    return 0;
}
