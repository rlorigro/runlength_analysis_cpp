
#include "QuadCompressor.hpp"
#include <vector>

using std::ofstream;
using std::make_shared;
using std::runtime_error;
using std::cout;
using std::cerr;
using std::tie;
using std::to_string;
using std::vector;


QuadCompressor::QuadCompressor(BoundingBox bounds){
    this->boundary = bounds;
}


QuadCompressor::QuadCompressor(){
    this->boundary = BoundingBox(QuadCoordinate(0,0), 0);
}


void QuadCompressor::update_loss_from_range(BoundingBox& bounds, QuadLoss& loss_calculator){
    //TODO: rewrite for non-square rectangles

    // Find boundaries, for constructing corners
    float x_left = this->boundary.center.x - this->boundary.half_size;
    float x_right = this->boundary.center.x + this->boundary.half_size;
    float y_bottom = this->boundary.center.y - this->boundary.half_size;
    float y_top = this->boundary.center.y + this->boundary.half_size;

    // Check if there is no intersection, terminate
    // TODO: change this lazy implementation to be a decision tree
    if (find_quadrant(QuadCoordinate(x_left,y_bottom)) == QuadTree::NOT_FOUND and
        find_quadrant(QuadCoordinate(x_left,y_top)) == QuadTree::NOT_FOUND and
        find_quadrant(QuadCoordinate(x_right,y_bottom)) == QuadTree::NOT_FOUND and
        find_quadrant(QuadCoordinate(x_right,y_top)) == QuadTree::NOT_FOUND){
        return;
    }

    // If there is an intersection and some of the points are in range, add them
    QuadCompressor query_quad(bounds);
    for (auto& point: this->points){
        if (query_quad.find_quadrant(point) != QuadTree::NOT_FOUND){
            loss_calculator.update(point);
        }
    }

    // If there are children search them
    if (this->subtrees[0] == nullptr) {
        return;
    }
    else{
        for (auto& subtree: this->subtrees){
            subtree->update_loss_from_range(bounds, loss_calculator);
        }
    }
}


void QuadCompressor::find_lossiest_leaf(QuadTree& reference_tree,
                                        map<double, map <double, pair <double,QuadTree*> > >& scores,
                                        QuadLoss& loss_calculator){
    size_t i = 0;
    double loss = 0;
    double x = reference_tree.boundary.center.x;
    double y = reference_tree.boundary.center.y;

    // If this is a leaf, calculate loss, and update max_loss
    if (reference_tree.subtrees[0] == nullptr) {
        bool score_exists = scores.count(x) != 0 and scores.at(x).count(y) != 0;

        if (not score_exists) {
            loss_calculator.reset();
            this->update_loss_from_range(reference_tree.boundary, loss_calculator);
            loss = loss_calculator.calculate_loss();
            scores[x][y] = make_pair(loss, &reference_tree);
        }
    }
    else {
        for (auto& subtree: reference_tree.subtrees) {
            this->find_lossiest_leaf(*subtree, scores, loss_calculator);
            i++;
        }
    }
}


shared_ptr<QuadTree> QuadCompressor::generate_child(BoundingBox bounds){
    return make_shared<QuadCompressor>(bounds);
}


void QuadCompressor::compress(uint64_t max_quadrants, QuadLoss& loss_calculator, path output_dir) {
    // Initialize a separate tree to represent custom loss-defined quadrants
    QuadCompressor reference_tree = QuadCompressor(this->boundary);

    map<double, map <double, pair <double,QuadTree*> > > scores;
    uint64_t n_quadrants = 1;
    uint64_t i = 0;
    while (n_quadrants < max_quadrants){
        cerr << n_quadrants << " quadrants calculated\n";
        subdivide_lossiest_leaf(loss_calculator, reference_tree, scores);
//        reference_tree.write_as_dot(output_dir,to_string(i), true);
        reference_tree.write_bounds(output_dir,to_string(i));
        n_quadrants += 3;
        i++;
    }
}


void QuadCompressor::subdivide_lossiest_leaf(QuadLoss& loss_calculator,
                                             QuadCompressor& reference_tree,
                                             map<double, map <double, pair <double,QuadTree*> > >& scores){
    ///
    /// loss_calculator must have an "update()" and "calculate_loss()" function.
    ///

    QuadTree* lossiest_leaf = nullptr;
    QuadTree* leaf;
    double max_loss = 0;
    double loss = 0;
    double x_max = -1;
    double y_max = -1;

    find_lossiest_leaf(reference_tree, scores, loss_calculator);

    for (auto& x_pair: scores){
        for (auto& y_pair: x_pair.second){
            tie(loss, leaf) = y_pair.second;

            if (loss > max_loss){
                max_loss = loss;
                lossiest_leaf = leaf;
                x_max = x_pair.first;
                y_max = y_pair.first;
                cout << x_max << " " << y_max << " " << loss << " " <<  lossiest_leaf << " MAX\n";
            }
        }
    }

    scores.at(x_max).erase(y_max);
    if (scores.at(x_max).empty()){
        scores.erase(x_max);
    }

    if (lossiest_leaf == nullptr) {
        throw runtime_error("ERROR: lossiest leaf not found. More bins than nodes?");
    }
    else {
        lossiest_leaf->subdivide();
    }
}
