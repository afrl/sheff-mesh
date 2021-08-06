#pragma once

#include "mesh/Defines.h"
#include "mesh/node/NodeGroup.h"
#include "simple_math/EigenDefines.h"
#include "utility/serialization/NamedValue.h"

#include <array>
#include <assert.h>
#include <set>
#include <vector>

namespace mesh
{
template<typename NC>
class NodeAccessorBase;

template<typename NodeAccessorType>
class NodeIteratorBase;

/// @brief Stores information pertaining to a set of nodes, including indices, points, and mesh
/// relationships. This class along with an ElementContainer object would compose information about
/// a mesh.
class NodeContainer
{
  public:
    using NodeAccessor      = NodeAccessorBase<NodeContainer>;
    using NodeConstAccessor = NodeAccessorBase<const NodeContainer>;
    using NodeIterator      = NodeIteratorBase<NodeAccessor>;
    using NodeConstIterator = NodeIteratorBase<NodeConstAccessor>;

    NodeContainer();

    bool operator==(const NodeContainer& other) const;

    /// @name Accessors
    ///@{
    /// Enables ranged-based loop over nodes.
    NodeIterator begin();
    /// Enables ranged-based loop over nodes (const).
    NodeConstIterator begin() const;
    /// Enables ranged-based loop over nodes.
    NodeIterator end();
    /// Enables ranged-based loop over nodes (const).
    NodeConstIterator end() const;
    /// Access a particular node by local index
    NodeAccessor operator[](LocalIndex i);
    /// Const access a particular node by local index
    NodeConstAccessor operator[](LocalIndex i) const;
    /// Returns the size of the container, which is the number of nodes on this partition.
    LocalIndex size() const;

    /// Access a particular group by its name, creates a new group if none exists by the given name.
    NodeGroup& group(const std::string group_name);
    /// Const access a particular group by its name
    const NodeGroup& group(const std::string group_name) const;
    /// Access the collection of node groups
    std::vector<NodeGroup>& groups();
    /// Const access the collection of node groups
    const std::vector<NodeGroup>& groups() const;
    ///@}

    /// @name Mutators
    ///@{
    /// Adds a node to the end of the container.  Assigns the global index of the new node to the
    /// next highest number (on this partition).
    void add(const Point& p);
    /// Adds a node to the end of the container.  If parent_partitions_ is empty, then the parent
    /// partition parameter will be ignored.
    void add(const Point& p, GlobalIndex i);
    /// Adds a node to the end of the container.
    void add(const Point& p, GlobalIndex i, PartitionIndex parent_partition);

    /// Removes a node by its index from the container without preserving ordering (and therefore
    /// local indices).
    void remove(LocalIndex i);
    /// Removes a node by its index from the container while preserving ordering (and therefore
    /// local indices).
    void remove_stable(LocalIndex i);
    /// Removes a list of nodes by their index from the container while preserving ordering (and
    /// therefore local indices).
    void remove_stable(std::vector<LocalIndex>& v);

    /// Swaps all data with the other container.
    void swap(NodeContainer& other);

    /// Permutes the nodes in the container to rearrange them with the given order.
    /// @param backward_permutation_indices - a map that goes from the new index to the old index of
    ///        the reordering.
    void permute(const std::vector<LocalIndex>& backward_permutation_indices);

    /// @brief Deletes all nodes from the container.
    void clear();
    ///@}

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("n_nodes_on_partition", global_indices_.size()),
          utility::named_value("global_indices", utility::slice(global_indices_)),
          utility::named_value("points", utility::slice(points_)),
          utility::named_value("parent_partitions", utility::slice(parent_partitions_)));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        std::size_t n_nodes;
        s(n_nodes);
        global_indices_.resize(n_nodes);
        points_.resize(n_nodes);
        parent_partitions_.resize(n_nodes);
        s(utility::Slice<GlobalIndex>(global_indices_.data(), global_indices_.size()),
          utility::Slice<Point>(points_.data(), points_.size()),
          utility::Slice<PartitionIndex>(parent_partitions_.data(), parent_partitions_.size()));
    }

    template<typename NC>
    friend class NodeAccessorBase;

  private:
    /// @name Required members
    ///@{
    /// Stores a unique identifier for each node.  The global set of indices is a 0-based,
    /// continuous range of integers.  In a distributed analysis, multiple partitions may have the
    /// same node. Expected size: # of nodes on partition.
    std::vector<GlobalIndex> global_indices_;
    /// Stores coordinates of each node in the global coordinate system used for the analysis.
    /// Expected size: 3 * # of nodes on partition.
    std::vector<Point> points_;
    /// Stores the parent partition of each  node.  Each node is owned by exactly one partition,
    /// while all nodes on a partition that are not owned are considered ghost nodes. Expected size:
    /// # of nodes on partition.
    std::vector<PartitionIndex> parent_partitions_;
    ///@}

    /// @name Group members
    ///@{
    ///
    std::vector<NodeGroup> groups_;
    ///@}
};
}
