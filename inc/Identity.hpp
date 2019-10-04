#ifndef RUNLENGTH_ANALYSIS_IDENTITY_HPP
#define RUNLENGTH_ANALYSIS_IDENTITY_HPP

#include <iostream>
#include <experimental/filesystem>

using std::cout;
using std::experimental::filesystem::path;


void measure_identity_from_fasta(path reads_fasta_path,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_threads);

void measure_identity_from_bam(path bam_path,
        path reference_fasta_path,
        uint16_t max_threads);

#endif //RUNLENGTH_ANALYSIS_IDENTITY_HPP
