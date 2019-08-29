
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"
#include <utility>
#include <iostream>

using boost::icl::interval_map;
using boost::icl::interval;
using std::make_pair;
using std::cout;


void test_base_level_interval_map(interval_map<int,int>& intervals){
    auto r0 = intervals.find(0);
    auto r1 = intervals.find(1);
    auto r2 = intervals.find(2);
    auto r5 = intervals.find(5);
    auto r6 = intervals.find(6);
    auto r7 = intervals.find(7);

    if (r0 != intervals.end()){
        cout << 0 << r0->first << r0->second << '\n';
    }

    if (r1 != intervals.end()){
        cout << 1 << r1->first << r1->second << '\n';
    }

    if (r2 != intervals.end()){
        cout << 2 << r2->first << r2->second << '\n';
    }

    if (r6 != intervals.end()){
        cout << 6 << r6->first << r6->second << '\n';
    }

    if (r5 != intervals.end()){
        cout << 5 << r5->first << r5->second << '\n';
    }

    if (r7 != intervals.end()){
        cout << 7 << r7->first << r7->second << '\n';
    }

}

int main(){
    auto intervals = interval_map<int,int>();

    auto a0 = interval<int>::right_open(0, 2);
    auto a2 = interval<int>::right_open(2, 4);
    auto a4 = interval<int>::right_open(4, 6);

    intervals.insert(make_pair(a0,1));
    intervals.insert(make_pair(a2,2));
    intervals.insert(make_pair(a4,3));

    auto nested_intervals = interval_map<int,interval_map<int,int>>();

    auto b0 = interval<int>::right_open(0, 2);
    auto b2 = interval<int>::right_open(2, 4);
    auto b4 = interval<int>::right_open(4, 6);

    nested_intervals.insert(make_pair(b0, intervals));
    nested_intervals.insert(make_pair(b2, intervals));
    nested_intervals.insert(make_pair(b4, intervals));

    auto r00 = nested_intervals.find(0);
    auto r01 = nested_intervals.find(1);
    auto r05 = nested_intervals.find(5);
    auto r06 = nested_intervals.find(6);
    auto r07 = nested_intervals.find(7);

    if (r00 != nested_intervals.end()){
        cout << 0 << r00->first << r00->second << '\n';
        test_base_level_interval_map(intervals);
    }

    if (r01 != nested_intervals.end()){
        cout << 1 << r01->first << r01->second << '\n';
        test_base_level_interval_map(intervals);
    }

    if (r06 != nested_intervals.end()){
        cout << 6 << r06->first << r06->second << '\n';
        test_base_level_interval_map(intervals);
    }

    if (r05 != nested_intervals.end()){
        cout << 5 << r05->first << r05->second << '\n';
        test_base_level_interval_map(intervals);
    }

    if (r07 != nested_intervals.end()){
        cout << 7 << r07->first << r07->second << '\n';
        test_base_level_interval_map(intervals);
    }

    return 0;
}
