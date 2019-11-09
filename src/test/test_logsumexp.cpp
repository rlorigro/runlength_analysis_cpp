
#include "Miscellaneous.hpp"
#include <cmath>

using std::log10;

int main(){

    double x1 = 1;      // 1
    double x2 = 10;     // 2
    double x3 = 100;    // 3
    double x4 = 1000;   // 4
    double x5 = 10000;  // 5
    double sum = 0;     // 0

    x1 = log10(x1);
    x2 = log10(x2);
    x3 = log10(x3);
    x4 = log10(x4);
    x5 = log10(x5);
    sum = log10(sum);

    cout << sum << " " << pow(10,sum) << '\n';

    sum = log10_sum_exp(sum, x1);
    cout << sum << " " << pow(10,sum) << '\n';

    sum = log10_sum_exp(sum, x2);
    cout << sum << " " << pow(10,sum) << '\n';

    sum = log10_sum_exp(sum, x3);
    cout << sum << " " << pow(10,sum) << '\n';

    sum = log10_sum_exp(sum, x4);
    cout << sum << " " << pow(10,sum) << '\n';

    sum = log10_sum_exp(sum, x5);
    cout << sum << " " << pow(10,sum) << '\n';

    return 0;
}