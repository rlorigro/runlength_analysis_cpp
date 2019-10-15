
#include "IterativeSummaryStats.hpp"
#include <iostream>
#include <vector>
#include <cmath>

using std::vector;
using std::cout;

int main(){
    auto stat_counter = IterativeSummaryStats<double>();

    vector<double> observations = {-2,0,2,4,6};

    for (auto& value: observations){
        stat_counter.add(value);
    }

    cout << "mean:\t\t" << stat_counter.get_mean() << '\n';
    cout << "variance:\t" << stat_counter.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter.get_variance()) << '\n';

    stat_counter.remove(-2);
    stat_counter.remove(6);

    cout << "mean:\t\t" << stat_counter.get_mean() << '\n';
    cout << "variance:\t" << stat_counter.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter.get_variance()) << '\n';

    auto stat_counter_int = IterativeSummaryStats<double>();

    for (auto& value: observations){
        stat_counter_int.add(value);
    }

    cout << "mean:\t\t" << stat_counter_int.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_int.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_int.get_variance()) << '\n';

    stat_counter_int.remove(-2);
    stat_counter_int.remove(6);

    cout << "mean:\t\t" << stat_counter_int.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_int.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_int.get_variance()) << '\n';

    return 0;
}
