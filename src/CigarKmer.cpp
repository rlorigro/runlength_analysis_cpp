#include "CigarKmer.hpp"


CigarKmer::CigarKmer(uint8_t k){
    this->k = k;
    this->cigar_counts = {};
    this->kmer = {};
}


void CigarKmer::update(Cigar& cigar, uint8_t base) {
    // Read move
    if (cigar.is_true_read_move()){

        // If the kmer is full
        if (this->kmer.size() == this->k){

            // Cycle out last base
            this->kmer.pop_front();

            // Cycle out last cigar count
            for(auto& [cigar_index, count]: this->cigar_counts.front()){
                this->total_cigar_counts[cigar_index] -= count;
            }
            this->cigar_counts.pop_front();
        }

        if (is_valid_base(base)) {
            // Add next base
            this->kmer.emplace_back(base);

            // Add next cigar count
            this->cigar_counts.emplace_back();
        }
    }

    if (not this->cigar_counts.empty()) {
        // Update the cigar counts (assume iteration of 1 ref/read base per step)
        this->total_cigar_counts[cigar.code]++;
        this->cigar_counts.back()[cigar.code]++;
    }
}


uint64_t CigarKmer::to_index(){
    return kmer_to_index(this->kmer);
}


KmerIdentities::KmerIdentities()=default;


KmerIdentities::KmerIdentities(uint8_t k){
    // Need as many arrays as there are kmers
    this->cigar_counts_per_kmer.resize(pow(4,k));
    this->k = k;
}


void operator+=(KmerIdentities& a, KmerIdentities& b){
    for (size_t kmer_index=0; kmer_index < b.cigar_counts_per_kmer.size(); kmer_index++){
        for (auto& [cigar_index, count]: b.cigar_counts_per_kmer[kmer_index]) {
            a.cigar_counts_per_kmer[kmer_index][cigar_index] += count;
        }
    }
}


void KmerIdentities::write_to_output(ofstream& output){
    for (size_t kmer_index=0; kmer_index < this->cigar_counts_per_kmer.size(); kmer_index++){
        output << kmer_index_to_string(kmer_index, this->k) << ',';

        for (size_t cigar_index=0; cigar_index < 9; cigar_index++) {
            if (not this->cigar_counts_per_kmer[kmer_index].empty()){
                output << this->cigar_counts_per_kmer[kmer_index][cigar_index];
            }
            else{
                output << to_string(0);
            }
            output << ",";

        }
        output << '\n';
        output << std::flush;
    }
}


void KmerIdentities::update(uint64_t kmer_index, array<uint64_t, 9>& cigar_counts){
    size_t i = 0;
    for (auto& c: cigar_counts){
        this->cigar_counts_per_kmer[kmer_index][i] += c;
        i++;
    }
}


void KmerIdentities::update(CigarKmer& cigar_kmer){
    // Only update if the kmer is full size
    if (cigar_kmer.kmer.size() == this->k){
        uint64_t kmer_index = cigar_kmer.to_index();

        // For every cigar map in the kmer, update the cigar count map for this entire kmer's index
        for (auto& [cigar_index, count]: cigar_kmer.total_cigar_counts) {
            this->cigar_counts_per_kmer[kmer_index][cigar_index] += count;
        }
    }
}
