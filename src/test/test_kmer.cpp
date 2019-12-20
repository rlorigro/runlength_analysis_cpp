#include "FastaWriter.hpp"
#include "PileupGenerator.hpp"
#include "SequenceElement.hpp"
#include "FastaReader.hpp"
#include "Align.hpp"
#include "Kmer.hpp"
#include "Base.hpp"

using std::cerr;
using std::ifstream;


int main(){
    uint8_t k = 6;
    string s = "ACGTAACCGGTTTTTT";
    deque<uint8_t> kmer;
    deque<uint8_t> kmer2;
    uint64_t index;

    for (size_t i=0; i<s.size(); i++){
        if (kmer.size() == k) {
            kmer.pop_front();
        }

        kmer.emplace_back(base_to_index(s[i]));

        cout << "----\n";
        cout << i << '\n';
        cout << "LENGTH: " << kmer.size() << '\n';
        for (auto& b: kmer){
            cout << int(b) << ' ';
        }
        cout << '\n';

        index = kmer_to_index(kmer);
        cout << index << '\n';

        index_to_kmer(kmer2, index, k);
        for (auto& b: kmer2){
            cout << int(b) << ' ';
        }
        cout << '\n';
    }

    cout << "MAX: " << std::pow(4,k) - 1 << '\n';

    return 0;
}
