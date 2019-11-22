#include <cmath>
#include "ShastaReader.hpp"


ShastaReader::ShastaReader(path directory_path, bool store_length_consensus, bool store_coverage_data){
    this->directory_path = directory_path;
    this->store_length_consensus = store_length_consensus;
    this->store_coverage_data = store_coverage_data;
}


void ShastaReader::index() {
    ///
    /// Load all the filenames of the csv
    ///
    path filename_prefix;
    string read_name;

    for (const path& file_path: directory_iterator(this->directory_path)){
        if (is_regular_file(file_path) and file_path.extension() == ".csv") {
            filename_prefix = file_path.filename();
            filename_prefix.replace_extension("");
            read_name = filename_prefix;
            replace(read_name.begin(), read_name.end(), '.', '_');

            this->file_paths[read_name] = file_path;
        }
    }
}


unordered_map<string,path> ShastaReader::get_index() {
    return this->file_paths;
}


void ShastaReader::set_index(unordered_map<string, path>& file_paths) {
    this->file_paths = file_paths;
}


void ShastaReader::read_file(CoverageSegment& segment, path& file_path) {
    ///
    /// Iterate all the lines in a shasta coverageData output file (CSV)
    ///
    // Clear the container
    segment = {};

    // Open file
    ifstream file = ifstream(file_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not open file " + file_path.string());
    }

    // File iteration variables
    string line;
    uint64_t l = 0;

    // TODO: stop using getline()
    while (getline(file, line)){
        this->parse_coverage_string(segment, line);
        l++;
    }
}


void ShastaReader::fetch_read(CoverageSegment& segment, string& read_name) {
    ///
    /// Fetch a read by its read name, which is derived from its filename
    ///
    if (this->file_paths.empty()){
        cerr << "Index not loaded for marginpolish directory, generating now...\n";
        this->index();
    }

    path file_path = this->file_paths.at(read_name);
    this->read_file(segment, file_path);
}


size_t ShastaReader::parse_consensus(CoverageSegment& shasta_segment, string& line){
    size_t n_commas = 0;
    size_t i = 0;
    string length_string;

    while (n_commas < 3) {
        if (line[i] == ','){
            n_commas++;
        }
        else if (n_commas == 1){
                shasta_segment.sequence += line[i];
        }
        else if (this->store_length_consensus and n_commas == 2){
                length_string += line[i];
        }
        i++;
    }

    // Only keep the consensus length data if specified by class attribute
    if (this->store_length_consensus) {
        uint16_t consensus_length = stoi(length_string);
        shasta_segment.lengths.emplace_back(consensus_length);
//        cout << line << '\n' << consensus_length << '\n';
    }


    return i-1;
}


void ShastaReader::parse_coverage_string(CoverageSegment& segment, string& line) {
    ///
    /// Reads a line of the Shasta runlength CSV and converts it to a vector of CoverageElements
    ///

    // For iterating elements within the tab separated string
    uint64_t start_index = parse_consensus(segment, line);
    uint64_t c = 0;

    // Placeholders for Coverage element
    char base;
    uint16_t length;
    bool reversal;
    float weight;

    // Buffers for elements with indeterminate number of chars
    string reversal_string;
    string weight_string;
    string length_string;
    bool space_found = false;
    float n_coverage = 0;

    vector<CoverageElement> pileup;

    // Iterate observations
    for (uint64_t i=start_index; i<=line.size()-1; i++){
        // Perform element level operations at start of next element and end of line
        if (line[i] == ',' or i == line.size()){
            if (i > start_index){
                reversal = this->parse_reversal_string(reversal_string);
                length = uint16_t(stoi(length_string));
                weight = stof(weight_string);
                pileup.emplace_back(base, length, reversal, weight);
                n_coverage += weight;
            }
            c = 0;
            length_string = "";
            weight_string = "";
            reversal_string = "";
            space_found = false;
        }
        // Parse element members
        else {
            if (c == 1) {
                base = line[i];
            }
            else if (not space_found and c > 1 and line[i] != ' ' and line[i] != '+' and line[i] != '-') {
                length_string += line[i];
            }
            else if (not space_found and c > 2 and line[i] != ' ') {
                reversal_string = line[i];
            }
            else if (line[i] == ' '){
                space_found = true;
            }
            else if (space_found){
                weight_string += line[i];
            }
        }
        c++;
    }

    // Either store all the specific coverage data, or just count the n-fold coverage at this position
    if (this->store_coverage_data) {
        segment.coverage_data.push_back(move(pileup));
    }
    else{
        segment.n_coverage.emplace_back(uint16_t(round(n_coverage)));
//        cout << n_coverage << '\n';
    }
}


bool ShastaReader::parse_reversal_string(string reversal_string) {
    ///
    /// Interpret string encoding of read alignment direction as a bool
    ///
    bool reversal;
    if (reversal_string == "+"){
        reversal = false;
    }
    else if (reversal_string == "-"){
        reversal = true;
    }
    else{
        throw runtime_error("ERROR: Invalid reversal string " + reversal_string);
    }
    return reversal;
}


void ShastaReader::read_consensus_sequence_from_file(CoverageSegment& segment, path& file_path) {
    // Open file
    ifstream file = ifstream(file_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not open file " + file_path.string());
    }

    // File iteration variables
    string line;
    uint64_t l = 0;

    // TODO: stop using getline()
    while (getline(file, line)){
        size_t base_index = line.find_first_of(',') + 1;
        segment.sequence += line[base_index];
        l++;
    }
}


void ShastaReader::fetch_consensus_sequence(CoverageSegment& segment, string& read_name) {
    ///
    /// Fetch a read by its read name (which is derived from its filename) and skip fetching all the coverage data
    ///

    // Clear the container
    segment = {};

    if (this->file_paths.empty()){
        cerr << "Index not loaded for directory, generating now...\n";
        this->index();
    }

    segment.name = read_name;

    path file_path = this->file_paths.at(read_name);
    this->read_consensus_sequence_from_file(segment, file_path);
}
