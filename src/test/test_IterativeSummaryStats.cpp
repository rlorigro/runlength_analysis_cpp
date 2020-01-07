
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

    cout << "\nTESTING ADD FUNCTION:\n";

    auto stat_counter_a = IterativeSummaryStats<double>();
    auto stat_counter_b = IterativeSummaryStats<double>();
    auto stat_counter_c = IterativeSummaryStats<double>();

    for (auto& value: observations){
        stat_counter_a.add(value);
        stat_counter_b.add(value);
        stat_counter_c.add(value);
        stat_counter_c.add(value);
    }

    stat_counter_b += stat_counter_a;

    cout << "A\n";
    cout << "mean:\t\t" << stat_counter_a.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_a.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_a.get_variance()) << '\n';
    cout << "B\n";
    cout << "mean:\t\t" << stat_counter_b.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_b.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_b.get_variance()) << '\n';
    cout << "C\n";
    cout << "mean:\t\t" << stat_counter_c.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_c.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_c.get_variance()) << '\n';

    cout << "\nTESTING ADD FUNCTION v2:\n";

    vector<double> observations_a = {1,1,1,0.5,0.5,0.5};
    vector<double> observations_b = {1,1,1};
    vector<double> observations_c = {0.5,0.5,0.5};

    stat_counter_a = {};
    stat_counter_b = {};
    stat_counter_c = {};
    auto stat_counter_d = IterativeSummaryStats<double>();

    for (auto& value: observations_a){
        stat_counter_a.add(value);
    }
    cout << "a\n";
    cout << "mean:\t\t" << stat_counter_a.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_a.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_a.get_variance()) << '\n';

    for (auto& value: observations_b){
        stat_counter_b.add(value);
    }
    cout << "b\n";
    cout << "mean:\t\t" << stat_counter_b.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_b.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_b.get_variance()) << '\n';

    for (auto& value: observations_c){
        stat_counter_c.add(value);
    }
    cout << "c\n";
    cout << "mean:\t\t" << stat_counter_c.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_c.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_c.get_variance()) << '\n';

    stat_counter_d += stat_counter_b;
    stat_counter_d += stat_counter_c;
    cout << "d\n";
    cout << "mean:\t\t" << stat_counter_d.get_mean() << '\n';
    cout << "variance:\t" << stat_counter_d.get_variance() << '\n';
    cout << "stddev:\t\t" << std::sqrt(stat_counter_d.get_variance()) << '\n';
    cout << "n:\t\t" << stat_counter_b.n << '\n';
    cout << "sum:\t\t" << stat_counter_b.sum << '\n';
    cout << "sum_of_squares:\t\t" << stat_counter_b.sum_of_squares << '\n';
    cout << "n:\t\t" << stat_counter_c.n << '\n';
    cout << "sum:\t\t" << stat_counter_c.sum << '\n';
    cout << "sum_of_squares:\t\t" << stat_counter_c.sum_of_squares << '\n';
    cout << "n:\t\t" << stat_counter_d.n << '\n';
    cout << "sum:\t\t" << stat_counter_d.sum << '\n';
    cout << "sum_of_squares:\t\t" << stat_counter_d.sum_of_squares << '\n';

    vector <IterativeSummaryStats <double> > stats = {stat_counter_b, stat_counter_c};
    cout << "Pooled variance:\t\t" << std::sqrt(pool_variances(stats)) <<'\n';
    cout << "Pooled variance:\t\t" << pool_means(stats) <<'\n';

    return 0;
}
