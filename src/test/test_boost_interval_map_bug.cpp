
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"
#include <utility>
#include <iostream>

using boost::icl::interval_map;
using boost::icl::interval;
using boost::icl::partial_enricher;
using boost::icl::total_enricher;
using std::make_pair;
using std::cout;


int main(){
    auto intervals = interval_map<double,char,total_enricher>();

    auto a0 = interval<double>::right_open(0, 1.9999);
    auto a1 = interval<double>::right_open(1.9999, 4);
    auto a2 = interval<double>::right_open(4, 6);


    cout << "ADDING PAIRS TO MAP:\n";
    cout << a1 << '\n';
    intervals.insert(make_pair(a1,-1));
    cout << intervals << '\n';

    cout << a0 << '\n';
    intervals.insert(make_pair(a0,0));
    cout << intervals << '\n';

    cout << a2 << '\n';
    intervals.insert(make_pair(a2,1));
    cout << intervals << '\n';

    auto r0 = intervals.find(1.0);
    auto r1 = intervals.find(3.0);
    auto r2 = intervals.find(5.0);
    auto r3 = intervals.find(-1.0);

    cout << "\nFINDING PAIRS:\n";

    cout << "finding 1: ";
    if (r0 != intervals.end()){
        cout << r0->first << int(r0->second) << '\n';
    }
    else{
        cout << "NOT FOUND!\n";
    }
    cout << "finding 3: ";
    if (r1 != intervals.end()){
        cout << r1->first << int(r1->second) << '\n';
    }
    else{
        cout << "NOT FOUND!\n";
    }
    cout << "finding 5: ";
    if (r2 != intervals.end()){
        cout << r2->first << int(r2->second) << '\n';
    }
    else{
        cout << "NOT FOUND!\n";
    }
    cout << "finding -1: ";
    if (r3 != intervals.end()){
        cout << r3->first << int(r3->second) << '\n';
    }
    else{
        cout << "NOT FOUND!\n";
    }


    return 0;
}
