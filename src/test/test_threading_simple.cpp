#include <thread>
#include <vector>
#include <iostream>

using std::vector;
using std::thread;
using std::cout;


void append_vector(vector <vector <uint16_t> >& v, uint16_t i){
    v[i].push_back(1);
    v[i].push_back(2);
    v[i].push_back(3);
}

void threaded_vector(uint16_t n_threads){
    vector <vector <uint16_t> > results(n_threads);
    vector <thread> threads;

    for (uint16_t i=0; i<n_threads; i++){
        threads.push_back(thread(append_vector, ref(results), i));
    }

    for (auto& item: threads){
        item.join();
    }

    for (auto& result: results){
        cout << "size " << result.size() << "\n";
        for (auto& item: result){
            cout << item << "\n";
        }
    }
}


int main(){
    uint16_t n_threads = 3;
    threaded_vector(n_threads);
}
