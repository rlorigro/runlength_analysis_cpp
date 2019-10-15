
#include "MultiDistributionStats.hpp"
#include "DiscreteWeibull.hpp"
#include "RunnieReader.hpp"
#include "QuadTree.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <experimental/filesystem>
#include "boost/program_options.hpp"

using std::cout;
using std::cerr;
using std::flush;
using std::vector;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;

void test_distribution_set(vector <vector<double> > distribution_set){
    MultiDistributionStats distribution_stats(distribution_set.size());

    for (auto& distribution: distribution_set){
        distribution_stats.update(distribution);
    }

    vector<double> means(distribution_set.size());
    vector<double> variances(distribution_set.size());

    cout << "POINTWISE MEAN:\n";
    distribution_stats.get_pointwise_mean(means);
    print_distribution(means);
    cout << '\n';

    cout << "POINTWISE VARIANCE:\n";
    distribution_stats.get_pointwise_variance(variances);
    print_distribution(variances);
    cout << '\n';

    cout << "Weighted variance:\t" << distribution_stats.get_weighted_variance() << '\n';
    cout << "Meta variance:\t" << distribution_stats.get_meta_variance() << '\n';
    cout << "Loss:\t" << distribution_stats.calculate_loss() << '\n';
}


int main(){
    vector <vector <double> > distributions_a = {
            {1,0,0,0,0,0},
            {0,1,0,0,0,0},
            {0,0,1,0,0,0},
            {0,0,0,1,0,0},
            {0,0,0,0,1,0},
            {0,0,0,0,0,1},
    };

    vector <vector <double> > distributions_b = {
            {1,0,0,0,0,0},
            {1,0,0,0,0,0},
            {1,0,0,0,0,0},
            {0,0,0,0,0,1},
            {0,0,0,0,0,1},
            {0,0,0,0,0,1},
    };

    vector <vector <double> > distributions_c = {
            {0,1,0,0,0,0},
            {0,1,0,0,0,0},
            {0,1,0,0,0,0},
            {0,0,0,0,1,0},
            {0,0,0,0,1,0},
            {0,0,0,0,1,0},
    };

    vector <vector <double> > distributions_d = {
            {0,0,1,0,0,0},
            {0,0,1,0,0,0},
            {0,0,1,0,0,0},
            {0,0,0,1,0,0},
            {0,0,0,1,0,0},
            {0,0,0,1,0,0},
    };

    cerr << "\n\nTESTING: distributions_a\n";
    test_distribution_set(distributions_a);

    cerr << "\n\nTESTING: distributions_b\n";
    test_distribution_set(distributions_b);

    cerr << "\n\nTESTING: distributions_c\n";
    test_distribution_set(distributions_c);

    cerr << "\n\nTESTING: distributions_d\n";
    test_distribution_set(distributions_d);

    return 0;
}
