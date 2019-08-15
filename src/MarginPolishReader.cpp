#include "MarginPolishReader.hpp"
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <experimental/filesystem>

using std::move;
using std::cout;
using std::flush;
using std::to_string;
using std::getline;
using std::replace;
using std::runtime_error;
using std::experimental::filesystem::directory_iterator;
using std::experimental::filesystem::path;


CoverageElement::CoverageElement(string base, uint16_t length, bool reversal, double weight){
    ///
    /// The lowest level object which represents one read's worth of observed alignment data at one position
    ///
    this->base = base;
    this->length = length;
    this->reversal = reversal;
    this->weight = weight;
}


const string CoverageElement::reversal_string = "FR";


string CoverageElement::to_string(){
    return this->base + std::to_string(this->length) + this->reversal_string[this->reversal] + std::to_string(this->weight);
}


uint8_t CoverageElement::get_base_index() {
    return base_to_index(this->base);
}


bool CoverageElement::is_conventional_base() {
    return is_valid_base(this->base);
}


MarginPolishReader::MarginPolishReader(path directory_path){
    this->directory_path = directory_path;
}


void MarginPolishSegment::print(){
    uint64_t i = 0;
    for (auto& coverage_pileup: this->coverage_data){
        cout << to_string(i) << " " << this->sequence[i] << " ";
        for (auto& coverage_element: coverage_pileup){
            cout << coverage_element.to_string() << " ";
        }
        cout << "\n";
        i++;
    }
}


void MarginPolishReader::index(){
    ///
    /// Load all the filenames of the tsv
    ///
    path filename_prefix;
    string read_name;

    for (const path& file_path: directory_iterator(this->directory_path)){
        filename_prefix = file_path.filename();
        filename_prefix.replace_extension("");
        read_name = filename_prefix;
        replace(read_name.begin(), read_name.end(), '.', '_');

        this->file_paths[read_name] = file_path;
    }
}


unordered_map<string,path> MarginPolishReader::get_index(){
    return this->file_paths;
}


void MarginPolishReader::set_index(unordered_map<string,path>& file_paths){
    this->file_paths = file_paths;
}


bool MarginPolishReader::parse_reversal_string(string reversal_string){
    ///
    /// Interpret string encoding of read alignment direction as a bool
    ///
    bool reversal;
    if (reversal_string == "+"){
        reversal = true;
    }
    else if (reversal_string == "-"){
        reversal = false;
    }
    else{
        throw runtime_error("ERROR: Invalid reversal string " + reversal_string);
    }
    return reversal;
}


void MarginPolishReader::parse_coverage_string(MarginPolishSegment& mp_segment, string& line){
    ///
    /// Read a line of the MarginPolish runlength TSV and converts it to a vector of CoverageElements
    ///

    // For iterating elements within the tab separated string
    uint64_t start_index = line.find_first_of('\t') + 2;
    uint64_t c = 0;

    // Append consensus base to MarginPolish segment
    mp_segment.sequence += line[start_index-1];

    // Placeholders for Coverage element
    string base;
    uint16_t length;
    bool reversal;
    double weight;

    // Buffers for elements with indeterminate number of chars
    string reversal_string;
    string weight_string;
    string length_string;
    bool comma_found = false;

    vector<CoverageElement> pileup;

    // Iterate observations
    for (uint64_t i=start_index; i<=line.size(); i++){
        // Perform element level operations at start of next element and end of line
        if (line[i] == '\t' or i == line.size()){
            if (i > start_index){
                reversal = this->parse_reversal_string(reversal_string);
                length = uint16_t(stoi(length_string));
                weight = stod(weight_string);
                CoverageElement coverage_element = CoverageElement(base, length, reversal, weight);
                pileup.push_back(move(coverage_element));
            }
            c = 0;
            length_string = "";
            weight_string = "";
            comma_found = false;
        }
        // Parse element members
        else {
            if (c == 1) {
                base = line[i];
            }
            else if (c == 2) {
                reversal_string = line[i];
            }
            else if (c > 2 and line[i] != ',' and not comma_found) {
                length_string += line[i];
            }
            else if (line[i] == ','){
                comma_found = true;
            }
            else if (comma_found){
                weight_string += line[i];
            }
        }
        c++;
    }
    mp_segment.coverage_data.push_back(move(pileup));
}


void MarginPolishReader::read_file(MarginPolishSegment& mp_segment, path& file_path){
    ///
    /// Iterate all the lines in a marginpolish runlength output file (TSV)
    ///
    // Clear the container
    mp_segment = {};

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
        if (l > 2){
            // Skip the first 2 header lines
            this->parse_coverage_string(mp_segment, line);
        }
        l++;
    }
}


void MarginPolishReader::fetch_read(MarginPolishSegment& mp_segment, string& read_name){
    ///
    /// Fetch a read by its read name, which is derived from its filename
    ///
    if (this->file_paths.empty()){
        this->index();
    }

    path file_path = this->file_paths.at(read_name);
    this->read_file(mp_segment, file_path);
}


void MarginPolishReader::read_consensus_sequence_from_file(MarginPolishSegment& mp_segment, path& file_path){
    ///
    /// Iterate all the lines in a marginpolish runlength output file (TSV)
    ///
    // Open file
    ifstream file = ifstream(file_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not open file " + file_path.string());
    }

    // File iteration variables
    string line;
    uint64_t l = 0;
    uint16_t tab_separator_index;

    while (getline(file, line)){
        if (l > 2){
            tab_separator_index = line.find_first_of('\t');
            mp_segment.sequence += line[tab_separator_index+1];
        }
        l++;
    }
}


void MarginPolishReader::fetch_consensus_sequence(MarginPolishSegment& mp_segment, string& read_name){
    ///
    /// Fetch a read by its read name (which is derived from its filename) and skip fetching all the coverage data
    ///
    // Clear the container
    mp_segment = {};

    if (this->file_paths.empty()){
        this->index();
    }
    mp_segment.name = read_name;

    path file_path = this->file_paths.at(read_name);
    this->read_consensus_sequence_from_file(mp_segment, file_path);
}
