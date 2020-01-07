
#include "Pileup.hpp"
#include "Kmer.hpp"
#include "Base.hpp"
#include <cmath>

using std::pow;
using std::min;


Pileup::Pileup() {
    this->n_channels = 0;
    this->pileup = vector <vector <vector <float> > >();
}


Pileup::Pileup(size_t n_channels, size_t region_size, size_t maximum_depth, vector<float>& default_read_data) {
    this->n_channels = n_channels;
    this->pileup = vector <vector <vector <float> > >(region_size, vector <vector <float> >(maximum_depth, default_read_data));
    this->coverage_per_position.resize(region_size, 0);
}


PileupReadKmerIterator::PileupReadKmerIterator(size_t window_size, size_t depth_index, uint8_t k){
    if (window_size % 2 != 1){
        throw runtime_error("ERROR: window size must be odd or it has no middle index");
    }

    this->n_operations = {0};
    this->window_size = window_size;
    this->depth_index = depth_index;
    this->width_index = 0;
    this->middle_index = window_size/2;     // index of middle kmer (c++ truncates int division, so -1 not needed)
    this->k = k;
}


void PileupReadKmerIterator::pop_left_kmers(Pileup& pileup){
    if ((this->n_operations.size() <= this->window_size)){
        return;
    }

    // Cycle out the left kmer indexes
    for (uint64_t i=0; i<this->n_operations[0]; i++){
        this->window_kmer_indexes.pop_front();
    }

    // Update the n_operations deque
    this->n_operations.pop_front();
}


void PileupReadKmerIterator::push_right_kmer(uint8_t base){
    if (not is_valid_base_index(base)){
        return;
    }

    uint64_t kmer_index;

    // Update kmer, ignoring anything not canonical
    if (is_valid_base_index(base)){
        this->current_kmer.emplace_back(base);
    }

    if (current_kmer.size() > this->k) {
        this->current_kmer.pop_front();

        // Get current kmer_index
        kmer_index = kmer_to_index(this->current_kmer);

        // And add the right side kmers
        this->window_kmer_indexes.emplace_back(kmer_index);

        // Update kmer count for this ref index
        this->n_operations.back()++;
    }
}


void PileupReadKmerIterator::update_middle_kmers(unordered_map<uint64_t, unordered_set <uint64_t> >& middle_kmers){
    if (this->n_operations.size() == this->window_size and this->n_operations[middle_index] != 0){
        size_t start_index=0;

        // Iterate up to the middle n_operations to find which kmers in window_kmer_indexes are in the middle of the
        // current window
        for (size_t i=0; i<middle_index; i++){
            start_index += this->n_operations[i];
        }

        // Add whichever kmer(s) is/are in the middle, possibly accounting for insert kmers
        for (size_t i=start_index; i<start_index+n_operations[middle_index]; i++){
            middle_kmers.emplace(this->window_kmer_indexes[i], unordered_set<uint64_t>());
        }
    }

}

///
/// want to provide the start and stop index, iterate the NEXT ref column and any insert columns, adding kmers to the
/// deque
///
/// \param pileup
void PileupReadKmerIterator::step(Pileup& pileup, unordered_map<uint64_t, unordered_set <uint64_t> >& middle_kmers){
    uint8_t base;

    // Make a placeholder counter for this step
    this->n_operations.emplace_back(0);

    // Load next ref-aligned base
    base = pileup.pileup[this->width_index][this->depth_index][Pileup::BASE];

    // Update the latest kmer if not in the edge case of forming a kmer
    this->pop_left_kmers(pileup);

    // Update current kmer
    push_right_kmer(base);

    // Load any insert bases that exist
    if (pileup.inserts.count(this->width_index) > 0) {
        for (auto &column: pileup.inserts.at(this->width_index)) {
            base = column[this->depth_index][Pileup::BASE];

            // Update current kmer
            push_right_kmer(base);
        }
    }

    this->update_middle_kmers(middle_kmers);

    this->width_index++;
}


PileupKmerIterator::PileupKmerIterator(Pileup& pileup, size_t window_size, uint8_t k){
    if (window_size % 2 != 1){
        throw runtime_error("ERROR: window size must be odd or it has no middle index");
    }

    // Initialize read iterators
    this->read_iterators = vector <PileupReadKmerIterator>();

    // Initialize all the read iterators for this pileup
    for (size_t i=0; i<pileup.max_observed_depth; i++){
        this->read_iterators.emplace_back(window_size, i, k);
    }

    this->window_size = window_size;
    this->k = k;
}


PileupKmerIterator::PileupKmerIterator(Pileup& pileup, size_t window_size, uint8_t k, KmerStats& kmer_stats){
    if (window_size % 2 != 1){
        throw runtime_error("ERROR: window size must be odd or it has no middle index");
    }

    // Initialize read iterators
    this->read_iterators = vector <PileupReadKmerIterator>();

    // Initialize all the read iterators for this pileup
    for (size_t i=0; i<pileup.max_observed_depth; i++){
        this->read_iterators.emplace_back(window_size, i, k);
    }

    this->window_size = window_size;
    this->k = k;
}


void PileupKmerIterator::step(Pileup& pileup){
    this->middle_kmers = {};

    // Step each of the read iterators
    for (size_t i=0; i<pileup.max_observed_depth; i++){
        this->read_iterators[i].step(pileup, this->middle_kmers);
    }

    // Find the number of rows that contain this kmer (not allowing multiple counts per row)
    for(size_t i=0; i<pileup.max_observed_depth; i++){
        for (auto& kmer_index: this->read_iterators[i].window_kmer_indexes){
            if (this->middle_kmers.count(kmer_index) > 0){
                this->middle_kmers.at(kmer_index).insert(i);  // A set tracks which rows in the pileup contain this kmer
            }
        }
    }

    this->width_index++;
}


void PileupKmerIterator::update_coverage_stats(Pileup& pileup, KmerStats& kmer_stats){
    double quality;
    uint16_t coverage;

    // Calculate the percentage of reads that contain each kmer and add it to the stats
    for (auto& [middle_kmer_index, kmer_coverage_set]: this->middle_kmers){
        coverage = pileup.coverage_per_position[this->width_index - this->read_iterators[0].middle_index];

        // Don't calculate kmer quality if the coverage is too low
        if (coverage < 5){
            continue;
        }

        // Calculate proportion of reads that contain the kmer
        quality = double(kmer_coverage_set.size())/double(coverage);

        // Occasional gaps in coverage smaller than the kmer size create impossible quality >1.0
        quality = min(1.0,quality);

        // Record this instance of quality
        kmer_stats.update_quality(middle_kmer_index, quality);
    }
}


void PileupKmerIterator::update_confusion_stats(PileupKmerIterator& ref_pileup_iterator, KmerConfusionStats& kmer_confusion_stats){
    for (auto& [ref_kmer_index,_]: ref_pileup_iterator.middle_kmers){
        for (auto& [read_kmer_index, supporting_reads]: this->middle_kmers){
//            cout << ref_kmer_index << " " << read_kmer_index << " " << supporting_reads.size() << '\n';
            kmer_confusion_stats.confusion[ref_kmer_index][read_kmer_index] += supporting_reads.size();
        }
    }
}
