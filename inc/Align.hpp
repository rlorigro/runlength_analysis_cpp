
#ifndef RUNLENGTH_ANALYSIS_ALIGN_HPP
#define RUNLENGTH_ANALYSIS_ALIGN_HPP

#include <string>
#include <iostream>
#include <stdexcept>
#include <experimental/filesystem>

using std::string;
using std::cout;
using std::cerr;
using std::runtime_error;
using std::experimental::filesystem::path;


// Call to minimap2 and any desired successive samtools calls. Returns the path of the output BAM/SAM.
path align(path ref_sequence_path,
           path read_sequence_path,
           path output_dir,
           bool sort = true,
           bool index = true,
           bool delete_intermediates = true,
           uint16_t k = 0,                     // Seed or kmer size (overrides any preset)
           string minimap_preset = "map-ont",   // Which preset to use (affects many params)
           bool explicit_mismatch = true,
           uint16_t max_threads = 1);


// Call to minimap2. Returns the path of the SAM file.
path minimap_align(path ref_sequence_path,
                   path read_sequence_path,
                   path output_dir,
                   string minimap_preset,       // Which preset to use (affects many params)
                   bool explicit_mismatch,
                   uint16_t max_threads,
                   uint16_t k=0);                  // Seed or kmer size (overrides any preset), if 0, no override


// Call to samtools sort. Returns the path of the sorted BAM file.
path samtools_sort(path input_path, uint16_t max_threads);


// Call to samtools index. Returns the path of the .bai file.
path samtools_index(path input_path, uint16_t max_threads);


#endif //RUNLENGTH_ANALYSIS_ALIGN_HPP
