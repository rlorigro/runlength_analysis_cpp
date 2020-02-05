
#include "Kmer.hpp"
#include "Base.hpp"

using std::runtime_error;
using std::to_string;


KmerStats::KmerStats()=default;


KmerStats::KmerStats(uint8_t k){
    this->qualities.resize(pow(4,k));
    this->k = k;
}


KmerConfusionStats::KmerConfusionStats()=default;

KmerConfusionStats::KmerConfusionStats(uint8_t k){
    this->k = k;
}


void operator+=(KmerConfusionStats& a, KmerConfusionStats& b){
    for (auto& [ref_kmer_index, submap]: b.confusion){
        for (auto& [read_kmer_index, frequency]: submap) {
            a.confusion[ref_kmer_index][read_kmer_index] += frequency;
        }
    }
}


string KmerConfusionStats::to_string(){
    string s;
    for (auto& [ref_kmer_index, submap]: this->confusion){
        double sum = 0;

        for (auto& [_, frequency]: submap) {
            sum += frequency;
        }

        for (auto& [read_kmer_index, frequency]: submap) {
            s += kmer_index_to_string(ref_kmer_index, k) + "," + kmer_index_to_string(read_kmer_index, k) + "," + std::to_string(double(frequency)/sum) + "," + std::to_string(uint64_t(sum)) + "\n";
        }
    }

    return s;
}


string KmerStats::to_string(){
    string s;
    uint64_t middle_kmer_index = 0;

    for (auto& stat: this->qualities) {
        if (not stat.empty()) {
            string kmer_string = kmer_index_to_string(middle_kmer_index, k);
            string mean_string = std::to_string(stat.get_mean());
            string variance_string = std::to_string(stat.get_variance());
            string n_string = std::to_string(stat.n);

            s += std::to_string(stat.n) + "," + kmer_string + "," + mean_string  + "," + variance_string + "," + n_string + '\n';
        }
        middle_kmer_index++;
    }

    return s;
}


void KmerStats::update_quality(uint64_t kmer_index, double quality){
    this->qualities.at(kmer_index).add(quality);
}


void operator+=(KmerStats& a, KmerStats& b){
    size_t n_kmers_a = a.qualities.size();
    size_t n_kmers_b = b.qualities.size();

    if (n_kmers_a != n_kmers_b){
        throw runtime_error("ERROR: cannot add KmerStats objects which do not have the same number of kmers or same k");
    }

    for (size_t i=0; i<n_kmers_a; i++){
        a.qualities[i] += b.qualities[i];
    }
}


uint64_t kmer_to_index(deque<uint8_t>& kmer){
    // First check if kmer is valid
    if (kmer.size() > 20){
        throw runtime_error("ERROR: cannot use kmer size greater than 20: " + to_string(kmer.size()));
    }

    // Initialize
    uint64_t index = 0;
    uint8_t k = 0;

    // For as many bases as are in the kmer, create a binary representation where each pair of bits comes from 1 base.
    // Later, this binary representation is interpreted as an integer;
    for (uint64_t base_index: kmer){
        base_index <<= k*2;
        index |= base_index;
        k++;
    }

    return index;
}


uint64_t kmer_to_index(deque<char>& kmer){
    // First check if kmer is valid
    if (kmer.size() > 20){
        throw runtime_error("ERROR: cannot use kmer size greater than 20: " + to_string(kmer.size()));
    }

    // Initialize
    uint64_t base_index = 0;
    uint64_t kmer_index = 0;
    uint8_t k = 0;

    // For as many bases as are in the kmer, create a binary representation where each pair of bits comes from 1 base.
    // Later, this binary representation is interpreted as an integer;
    for (char base: kmer){
        base_index = base_to_index(base);
        base_index <<= k*2;
        kmer_index |= base_index;
        k++;
    }

    return kmer_index;
}


void index_to_kmer(deque<uint8_t>& kmer, uint64_t index, uint8_t k){
    kmer = {};
    uint8_t base_index = 0;
    uint64_t mask = 3;

    for (uint8_t i=0; i<k; i++){
        base_index = (index >> i*2);
        base_index &= mask;
        kmer.emplace_back(base_index);
    }
}


string kmer_index_to_string(uint64_t kmer_index, uint8_t k){
    deque<uint8_t> kmer;
    string s;

    index_to_kmer(kmer, kmer_index, k);
    for (auto& b: kmer) {
        s += index_to_base(b);
    }

    return s;
}

