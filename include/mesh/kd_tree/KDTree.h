#pragma once

#include "mesh/Mesh.h"
#include "mesh/node/NodeContainer.h"
#include "nanoflann.h"

#include <vector>

namespace mesh
{
// An  this is the "dataset to kd-tree" adaptor class:
class NodeContainerAdaptor
{
  public:
    NodeContainerAdaptor(const NodeContainer& nodes)
      : nodes_(nodes)
    {
    }
    // NodeContainerAdaptor(const NodeContainerAdaptor&)     = default;
    NodeContainerAdaptor(NodeContainerAdaptor&&) noexcept = default;

    /// Returns the number of points in the underlying dataset (NodeContainer)
    inline size_t kdtree_get_point_count() const
    {
        return nodes_.size();
    }

    // Returns the indicated component of the position for node with index "point_index"
    inline double kdtree_get_pt(const size_t point_index, const size_t component) const
    {
        // The index can be safely cast to component, since the underlying dataset is a
        // NodeContainer
        switch (component)
        {
            case 0: return nodes_[static_cast<mesh::LocalIndex>(point_index)].x();
            case 1: return nodes_[static_cast<mesh::LocalIndex>(point_index)].y();
            default:
            case 2: return nodes_[static_cast<mesh::LocalIndex>(point_index)].z();
        }
    }

    /// Indicates if the bounding-box is already computed (returns false to force nanoflann to
    /// compute it).
    template<class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }

    /// Reference to underlying NodeContainer on which to build the KD tree
    const NodeContainer& nodes_;
};

/// @brief Provides a KD-tree that is built directly on a NodeContainer (wraps the nanoflann KD-tree
/// library)
/// Usage:
/// \code{cpp}
/// Mesh m; // Assume it is filled
/// NodeContainerKDTree kd_tree(m.nodes());
/// auto closest_node = kd_tree.closest_node({ 0, 0, 0 });
/// \endcode
class NodeContainerKDTree
{
  public:
    NodeContainerKDTree(const NodeContainer& nodes);
    // NodeContainerKDTree(const NodeContainerKDTree& other);
    NodeContainerKDTree(NodeContainerKDTree&&) noexcept = default;

    /// Finds the closest node to the given point and returns an accessor.
    NodeConstAccessor closest_node(const Point& p) const;

    /// Finds the closest 'n' nodes to the given point and returns a subrange iterator.
    IteratorRange<NodeSubrangeConstIterator> closest_nodes(const Point& p, const int n) const;

    /// Finds the nodes within the given distance to the given point and returns a subrange iterator.
    IteratorRange<NodeSubrangeConstIterator> in_radius(const Point& p,
                                                       const double& search_radius) const;

  private:
    using NanoflannKDType = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, NodeContainerAdaptor>,
      NodeContainerAdaptor /* container */,
      3 /* dimensionality */>;

    // Note: adapter_ must appear before index_ due to initializer order for some compilers
    NodeContainerAdaptor adapter_;
    NanoflannKDType index_;
};
}