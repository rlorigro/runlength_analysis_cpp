
#include "QuadTree.hpp"
#include <stdexcept>
#include <iostream>

using std::make_shared;
using std::runtime_error;
using std::cout;
using std::to_string;

QuadCoordinate::QuadCoordinate()=default;


QuadCoordinate::QuadCoordinate(float x, float y){
    this->x = x;
    this->y = y;
}


BoundingBox::BoundingBox()=default;


BoundingBox::BoundingBox(QuadCoordinate center, float half_size){
    this->center = center;
    this->half_size = half_size;
}


QuadTree::QuadTree(BoundingBox boundary){
    this->boundary = boundary;
}


void QuadTree::redistribute_points(){
    uint8_t index;

    while (not this->points.empty()){
        index = this->find_quadrant(this->points.back());
        this->subtrees[index]->insert(this->points.back());
        this->points.pop_back();
    }
}


void QuadTree::subdivide(){
    // Find new half size
    float new_half_size = this->boundary.half_size/2;

    // Find boundaries
    float x_left = this->boundary.center.x - new_half_size;
    float x_right = this->boundary.center.x + new_half_size;
    float y_bottom = this->boundary.center.y - new_half_size;
    float y_top = this->boundary.center.y + new_half_size;

    // Top left
    QuadCoordinate top_left_center = QuadCoordinate(x_left, y_top);
    BoundingBox top_left_bounds = BoundingBox(top_left_center, new_half_size);
    this->subtrees[QuadTree::TOP_LEFT] = make_shared<QuadTree>(QuadTree(top_left_bounds));

    // Top right
    QuadCoordinate top_right_center = QuadCoordinate(x_right, y_top);
    BoundingBox top_right_bounds = BoundingBox(top_right_center, new_half_size);
    this->subtrees[QuadTree::TOP_RIGHT] = make_shared<QuadTree>(QuadTree(top_right_bounds));

    // Bottom left
    QuadCoordinate bottom_left_center = QuadCoordinate(x_left, y_bottom);
    BoundingBox bottom_left_bounds = BoundingBox(bottom_left_center, new_half_size);
    this->subtrees[QuadTree::BOTTOM_LEFT] = make_shared<QuadTree>(QuadTree(bottom_left_bounds));

    // Bottom right
    QuadCoordinate bottom_right_center = QuadCoordinate(x_right, y_bottom);
    BoundingBox bottom_right_bounds = BoundingBox(bottom_right_center, new_half_size);
    this->subtrees[QuadTree::BOTTOM_RIGHT] = make_shared<QuadTree>(QuadTree(bottom_right_bounds));

    this->redistribute_points();
}


uint8_t QuadTree::find_quadrant(QuadCoordinate c){
    uint8_t quadrant = QuadTree::NOT_FOUND;

    // Top
    if (c.y >= this->boundary.center.y){
        // Inside top bounds
        if (c.y < this->boundary.center.y + this->boundary.half_size) {
            // Right
            if (c.x >= this->boundary.center.x){
                // Inside right bounds
                if (c.x < this->boundary.center.x + this->boundary.half_size) {
                    // Top right!
                    quadrant = QuadTree::TOP_RIGHT;
                }
            }

            // Left
            else {
                // Inside left bounds
                if (c.x >= this->boundary.center.x - this->boundary.half_size) {
                    // Top left!
                    quadrant = QuadTree::TOP_LEFT;
                }
            }
        }
    }
    // Bottom
    else {
        // Inside bottom bounds
        if (c.y >= this->boundary.center.y - this->boundary.half_size) {
            // Right
            if (c.x >= this->boundary.center.x){
                // Inside right bounds
                if (c.x < this->boundary.center.x + this->boundary.half_size) {
                    // Bottom right!
                    quadrant = QuadTree::BOTTOM_RIGHT;
                }
            }

            // Left
            else {
                // Inside left bounds
                if (c.x >= this->boundary.center.x - this->boundary.half_size) {
                    // Bottom left!
                    quadrant = QuadTree::BOTTOM_LEFT;
                }
            }
        }
    }

    return quadrant;
}


bool QuadTree::insert(QuadCoordinate c){
    bool success = false;

    // Check whether the point is within one of the quadrants of this tree, if so return the quadrant
    uint8_t index = this->find_quadrant(c);

    if (index != QuadTree::NOT_FOUND) {
        // If this quadtree is below capacity and it has not formed any children, just append the point
        if (this->points.size() < QuadTree::capacity and this->subtrees[0] == nullptr) {
            this->points.push_back(c);
            success = true;
        }

        // If capacity has been exceeded, and children need to be spawned, spawn them and pass the point onward
        else if (this->subtrees[0] == nullptr){
            this->subdivide();
            success = this->subtrees[index]->insert(c);
        }

        // If this tree already has children and exceeded capacity, pass the point onward
        else {
            success = this->subtrees[index]->insert(c);
        }
    }

    return success;
}


void QuadTree::append_dot_string(string& edges, string& labels, uint64_t n){
    if (n>100){
        throw runtime_error("NO.");
    }

    string name;
    string label;
    uint64_t n_child;

    size_t i = 0;
    for (auto& subtree: this->subtrees){
        if (not subtree){
            continue;
        }

        n_child = n + i + 1;
        label = to_string(subtree->points.size());
        name = to_string(n_child);

        edges += "\n\t" + to_string(n) + "->" + name + ";";

        if (not subtree->points.empty()){
            labels += "\t" + name + "[label=\"" + label + "\"];\n";
        }

        subtree->append_dot_string(edges, labels, n_child);

        i++;
    }
}


void QuadTree::write_as_dot(path output_dir){
    string dot_string = "digraph G {\n";
    string edges;
    string labels;

    uint64_t n = 0;

    this->append_dot_string(edges, labels, n);

    dot_string += labels;
    dot_string += edges;
    dot_string += "\n}\n";

    cout << dot_string;
}
