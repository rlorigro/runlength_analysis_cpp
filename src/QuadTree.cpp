
#include "QuadTree.hpp"
#include <stdexcept>
#include <iostream>

using std::make_shared;
using std::runtime_error;
using std::cout;


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

    uint8_t index = this->find_quadrant(c);

    if (index != QuadTree::NOT_FOUND) {
        cout << "quadrant: " << int(index) << '\n';

        if (this->points.size() < QuadTree::capacity and this->subtrees[0] == nullptr) {
            this->points.push_back(c);
            success = true;
        }
        else if (this->subtrees[0] == nullptr){
            cout << "subdividing :) \n";
            this->subdivide();
            success = this->subtrees[index]->insert(c);
        }
        else {
            success = this->subtrees[index]->insert(c);
        }
    }

    return success;
}

