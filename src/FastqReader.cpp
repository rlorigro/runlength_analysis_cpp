
#include "FastqReader.hpp"

using std::to_string;


FastqReader::FastqReader(path file_path){
    this->file_path = file_path;

    ifstream test(file_path);

    if (not test.good()){
        throw runtime_error("ERROR: file read error: " + this->file_path.string());
    }
}


void FastqReader::parse_filter_candidate(
        string& name,
        string& sequence,
        string& qualities,
        IterativeSummaryStats<uint64_t>& stats,
        uint64_t minimum_average_quality,
        ofstream& output_file){

    if (stats.get_mean() >= double(minimum_average_quality)){
        output_file << '@' << name << '\n';
        output_file << sequence << '\n';
        output_file << '+' << '\n';
        output_file << qualities << '\n';
    }

    // Reset the temporary containers
    name.clear();
    sequence.clear();
    qualities.clear();
    stats.clear();
}


void FastqReader::filter_by_quality(uint64_t minimum_average_quality) {
    string output_filename = this->file_path.stem().string();
    output_filename = output_filename + "_filtered_" + to_string(int(minimum_average_quality)) + ".fastq";

    path parent_directory = this->file_path.parent_path();
    path output_path = parent_directory / output_filename;

    ifstream input_file(this->file_path);
    ofstream output_file(output_path);

    if (not output_file.is_open()){
        throw runtime_error("ERROR: file can't be written: " + output_path.string());
    }

    string name;
    string sequence;
    string qualities;
    IterativeSummaryStats<uint64_t> stats;

    uint64_t line_index = 0;
    bool first_char = true;
    char c;

    // Read each char
    while (input_file.get(c)) {
        // When newline found, increment line number
        if (c == '\n'){
            line_index++;
            first_char = true;
        }

        // Header, name
        else if ((line_index % 4) == 0){
            // Skip the header character
            if (first_char){
                if (not sequence.empty()) {
                    this->parse_filter_candidate(
                            name,
                            sequence,
                            qualities,
                            stats,
                            minimum_average_quality,
                            output_file);
                }

                first_char = false;
                continue;
            }

            // Update the name
            name += c;
        }
        else if ((line_index % 4) == 1){
            // Update the sequence
            sequence += c;
        }
        else if ((line_index % 4) == 2){
            // Skip useless garbage (why does this line exist?)
            continue;
        }
        else if ((line_index % 4) == 3){
            // Update the qualities
            qualities += c;
            stats.add((uint64_t(c)-33));
        }
    }

    if (not sequence.empty()){
        this->parse_filter_candidate(
                name,
                sequence,
                qualities,
                stats,
                minimum_average_quality,
                output_file);
    }
}
