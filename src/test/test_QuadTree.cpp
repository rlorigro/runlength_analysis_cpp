
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
//            {-0.0,-0.0},
            };

    vector<bool> truth_set = {true, true, true,
                              true, true, true,
                              true, true, true, true,
                              false, false,
                              false, false,
                              true, false,
//                              true
    };

    float x = 0.09;
    float y = 0.09;

    float size = 0.18;

    // Bottom right
    QuadCoordinate center = QuadCoordinate(x, y);
    BoundingBox bounds = BoundingBox(center, size/2);
    QuadTree tree = QuadTree(bounds);

    cerr << "TESTING ADDITION OF POINTS IN/OUT BOUNDS\n";

    bool success;
    size_t i = 0;
    for (auto& [x_i, y_i]: test_points){
        center = QuadCoordinate(x_i, y_i);
        success = tree.insert(center);

        if (not success == truth_set[i]){
            throw runtime_error("FAIL: " + to_string(x_i) + " " + to_string(y_i));
        }
        i++;
    }
    cerr << "PASS\n";

    cerr << "\nTESTING FETCHING INDIVIDUAL POINTS IN/OUT BOUNDS\n";
    i = 0;
    size_t n_true = 0;
    for (auto& [x_i, y_i]: test_points) {
        vector<QuadCoordinate> results;
        center = QuadCoordinate(x_i, y_i);
        BoundingBox test_bounds = BoundingBox(center, 0.0001);

        tree.query_range(results, test_bounds);

//        cout << results.size() << " " << truth_set[i] << '\n';
        if (not (results.size() == truth_set[i])){
            cerr << "FOUND INCORRECT NUMBER OF POINTS\n";
            for (auto& point: results){
                cerr << point.x << " " << point.y << '\n';
            }
            throw runtime_error("FAIL: " + to_string(x_i) + " " + to_string(y_i));
        }

        n_true += truth_set[i];
        i++;
    }
    cerr << "PASS\n";
    cerr << "\nTESTING FETCHING ALL POINTS IN BOUNDS\n";

    vector<QuadCoordinate> results;
    tree.query_range(results, bounds);

    if (not (results.size() == n_true)){
        throw runtime_error("FAIL: Not all points found");
    }

    cerr << "PASS\n";

    tree.write_as_dot("output/");
    tree.write_bounds("output/");

    return 0;
}
