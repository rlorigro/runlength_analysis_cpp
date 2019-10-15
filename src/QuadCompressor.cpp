
#include "QuadCompressor.hpp"
#include <vector>

using std::ofstream;
using std::make_shared;
using std::runtime_error;
using std::cout;
using std::to_string;
using std::vector;


QuadCompressor::QuadCompressor(BoundingBox bounds){
    this->boundary = bounds;
}


QuadCompressor::QuadCompressor(){
    this->boundary = BoundingBox(QuadCoordinate(0,0), 0);
}


void QuadCompressor::find_lossiest_leaf(double& max_loss, shared_ptr<QuadCompressor>& lossiest_leaf, QuadLoss& loss_calculator){
    size_t i = 0;
    double loss = 0;

    // If this is a leaf, calculate loss, and update max_loss
    if (not this->subtrees[0]) {
        loss_calculator.reset();
        update_loss_from_range(this->boundary, loss_calculator);
        loss = loss_calculator.calculate_loss();

        if (loss > max_loss){
            max_loss = loss;
            lossiest_leaf = make_shared<QuadCompressor>(*this);
        }
    }
    else {
        for (auto &subtree: this->subtrees) {
            if (not subtree) {
                continue;
            }

            subtree->find_lossiest_leaf(max_loss, lossiest_leaf, loss_calculator);

            i++;
        }
    }
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
    QuadTree query_quad(bounds);
    for (auto& point: this->points){
        if (query_quad.find_quadrant(point) != QuadTree::NOT_FOUND){
            loss_calculator.update(point);
        }
    }

    // If there are children search them
    if (this->subtrees[0] == nullptr){
        return;
    }
    else{
        for (auto& subtree: this->subtrees){
            subtree->update_loss_from_range(bounds, loss_calculator);
        }
    }
}


void QuadCompressor::subdivide_most_lossy_quadrant(QuadLoss& loss_calculator){
    ///
    /// loss_calculator must have an "update()" and "calculate_loss()" function.
    ///

    // Initialize a separate tree to represent custom loss-defined quadrants
    QuadTree reference_tree = QuadTree(this->boundary);
    reference_tree.subdivide();

    double max_loss = 0;
    shared_ptr<QuadCompressor> lossiest_leaf;

    find_lossiest_leaf(max_loss, lossiest_leaf, loss_calculator);

    lossiest_leaf->subdivide();
}
