
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

void test();

path align(path ref_sequence_path,
           path read_sequence_path,
           path output_dir,
           bool sort = true,
           bool index = true,
           bool delete_intermediates = true,
           uint16_t k = 15,                     // Seed or kmer size (overrides any preset)
           string minimap_preset = "map-ont",   // Which preset to use (affects many params)
           uint16_t max_threads = 1);

#endif //RUNLENGTH_ANALYSIS_ALIGN_HPP
