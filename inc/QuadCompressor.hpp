
#ifndef RUNLENGTH_ANALYSIS_QUADCOMPRESSION_H
#define RUNLENGTH_ANALYSIS_QUADCOMPRESSION_H

#include "MultiDistributionStats.hpp"
#include "QuadTree.hpp"
#include "QuadLoss.hpp"
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <experimental/filesystem>

using std::vector;
using std::array;
using std::shared_ptr;
using std::string;
using std::experimental::filesystem::path;


class QuadCompressor: public QuadTree{
public:
    /// Attributes ///
    array<shared_ptr<QuadCompressor>,4> subtrees;

    /// Methods ///
//    using QuadTree::QuadTree;
    QuadCompressor();
    QuadCompressor(BoundingBox bounds);
    void update_loss_from_range(BoundingBox& bounds, QuadLoss& loss_calculator);
    void subdivide_most_lossy_quadrant(QuadLoss& loss_calculator);
    void find_lossiest_leaf(double& max_loss, shared_ptr<QuadCompressor>& lossiest_leaf, QuadLoss& loss_calculator);
};


#endif //RUNLENGTH_ANALYSIS_QUADCOMPRESSION_H
