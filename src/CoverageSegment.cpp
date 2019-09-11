
#include "CoverageSegment.hpp"


void CoverageSegment::print(){
    uint64_t i = 0;
    for (auto& coverage_pileup: this->coverage_data){
        cout << to_string(i) << " " << this->sequence[i] << " ";
        for (auto& coverage_element: coverage_pileup){
            cout << coverage_element.to_string() << " ";
        }
        cout << "\n";
        i++;
    }
}
