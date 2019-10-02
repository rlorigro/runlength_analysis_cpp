
#include "QuadTree.hpp"
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>

using std::vector;
using std::pair;
using std::cout;
using std::cerr;
using std::runtime_error;
using std::to_string;


int main(){
    vector <pair <float,float> > test_points = {
            {0.05,0.05}, {0.05,0.14}, {0.14,0.14},
            {0.11,0.02}, {0.11,0.08}, {0.17,0.08},
            {0.14,0.01}, {0.14,0.03}, {0.17,0.03},{0.17,0.01},
            {-0.5,-0.5}, {1.5,1.5},
            {0.0,-0.5}, {1.0,1.5},
            {0.0,0.0}, {1.0,1.0},
            {-0.0,-0.0},
            };

    vector<bool> truth_set = {true, true, true,
                              true, true, true,
                              true, true, true, true,
                              false, false,
                              false, false,
                              true, false,
                              true};

    float x = 0.09;
    float y = 0.09;

    float size = 0.18;

    // Bottom right
    QuadCoordinate center = QuadCoordinate(x, y);
    BoundingBox bounds = BoundingBox(center, size/2);
    QuadTree tree = QuadTree(bounds);

    bool success;
    size_t i = 0;

    for (auto& [x_i, y_i]: test_points){
        cout << '\n' << x_i << " " << y_i << '\n';
        center = QuadCoordinate(x_i, y_i);
        success = tree.insert(center);

        if (not success == truth_set[i]){
            throw runtime_error("FAIL: " + to_string(i));
        }
        i++;
    }
    cerr << "PASS\n";


    return 0;


}
