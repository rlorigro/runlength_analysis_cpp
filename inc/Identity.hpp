#ifndef RUNLENGTH_ANALYSIS_IDENTITY_HPP
#define RUNLENGTH_ANALYSIS_IDENTITY_HPP

#include <iostream>
#include <string>
#include <experimental/filesystem>
#include "BamReader.hpp"

using std::cout;
using std::string;
using std::experimental::filesystem::path;


CigarStats measure_identity_from_fasta(path reads_fasta_path,
        path reference_fasta_path,
        path output_directory,
        string minimap_preset,
        uint16_t max_threads,
        uint64_t chunk_size=1*1000*1000,
        vector<Region> regions={});

CigarStats measure_identity_from_bam(path bam_path,
        path reference_fasta_path,
        uint16_t max_threads,
        uint64_t chunk_size=1*1000*1000);


#endif //RUNLENGTH_ANALYSIS_IDENTITY_HPP
