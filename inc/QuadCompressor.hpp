
#ifndef RUNLENGTH_ANALYSIS_QUADCOMPRESSION_H
#define RUNLENGTH_ANALYSIS_QUADCOMPRESSION_H

#include "MultiDistributionStats.hpp"
#include "QuadTree.hpp"
#include "QuadLoss.hpp"
#include <utility>
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <experimental/filesystem>
#include <map>

using std::pair;
using std::make_pair;
using std::vector;
using std::array;
using std::shared_ptr;
using std::string;
using std::experimental::filesystem::path;
using std::map;


class QuadCompressor: public QuadTree{
public:
    /// Attributes ///
//    array<shared_ptr<QuadCompressor>,4> subtrees;

    /// Methods ///
//    using QuadTree::QuadTree;
    QuadCompressor();
    QuadCompressor(BoundingBox bounds);
    void compress(uint64_t max_quadrants, QuadLoss& loss_calculator, path output_dir);
    void update_loss_from_range(BoundingBox& bounds, QuadLoss& loss_calculator);
    shared_ptr<QuadTree> generate_child(BoundingBox bounds) override;

    void subdivide_lossiest_leaf(QuadLoss& loss_calculator,
                                 QuadCompressor& reference_tree,
                                 map<double, map <double, pair <double,QuadTree*> > >& scores);

    void find_lossiest_leaf(QuadTree& reference_tree,
                            map<double, map <double, pair <double,QuadTree*> > >& scores,
                            QuadLoss& loss_calculator);
};


#endif //RUNLENGTH_ANALYSIS_QUADCOMPRESSION_H
