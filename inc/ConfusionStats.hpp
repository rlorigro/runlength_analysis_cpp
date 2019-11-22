
#ifndef RUNLENGTH_ANALYSIS_CONFUSIONSTATS_HPP
#define RUNLENGTH_ANALYSIS_CONFUSIONSTATS_HPP

#include <map>
#include <string>
#include <vector>
#include <experimental/filesystem>

using std::map;
using std::string;
using std::vector;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;


class ConfusionStats {
    ///
    /// Track confusion vs coverage for base and length confusion
    ///

public:
    map <uint16_t,map <uint16_t,uint64_t> > length_match_coverage;
    map <uint16_t,map <uint16_t,uint64_t> > length_mismatch_coverage;

    map <uint8_t,map <uint16_t,uint64_t> > base_match_coverage;
    map <uint8_t,map <uint16_t,uint64_t> > base_mismatch_coverage;

    void update(char true_base,
        char consensus_base,
        uint16_t true_length,
        uint16_t consensus_length,
        uint16_t n_coverage);

    uint16_t find_max_coverage();
    void write_to_file(path output_file_path);
    void write_summary_to_file(path output_file_path);
};


void operator+=(ConfusionStats& matrix_a, ConfusionStats& matrix_b);

void measure_confusion_stats_from_shasta(path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_threads,
        path bed_path=path());


#endif //RUNLENGTH_ANALYSIS_CONFUSIONSTATS_HPP
