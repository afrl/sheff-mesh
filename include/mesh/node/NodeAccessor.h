#pragma once

#include "mesh/Defines.h"
#include "mesh/node/NodeContainer.h"

#include <assert.h>
#include <vector>

namespace mesh
{
/// @brief Provides an interface to access node information.
/// Templated to allow a const and non-const version of the accessor without duplicate code.
template<typename NC>
class NodeAccessorBase
{
  public:
    using NodeContainerType = NC;

    NodeAccessorBase(NC& nodes, LocalIndex starting_index = 0)
      : nodes_(nodes)
      , curr_index_(starting_index)
    {
    }

    /// Returns the index of the current node within the NodeContainer that contains it.
    inline const LocalIndex& local_index() const
    {
        return curr_index_;
    }

    /// Returns the globally unique index for the current node.
    inline auto& global_index() const
    {
        return nodes_.global_indices_[curr_index_];
    }

    /// Returns the first component of the coordinate (in the global coordinate system)
    inline auto& x() const
    {
        return nodes_.points_[curr_index_][0];
    }

    /// Returns the second component of the coordinate (in the global coordinate system)
    inline auto& y() const
    {
        return nodes_.points_[curr_index_][1];
    }

    /// Returns the third component of the coordinate (in the global coordinate system)
    inline auto& z() const
    {
        return nodes_.points_[curr_index_][2];
    }

    /// Returns point of node
    inline Point point() const
    {
        return nodes_.points_[curr_index_];
    }

    /// Returns a map of the coordinates to a 3D Eigen vector. It can be used as a normal Eigen
    /// vector.
    inline Eigen::Map<const Eigen::Vector3d> coord() const
    {
        return Eigen::Map<const Eigen::Vector3d>(&nodes_.points_[curr_index_][0], 3);
    }

    /// Returns the parent partition of the current node
    inline auto& parent_partition() const
    {
        assert(!nodes_.parent_partitions_.empty());
        return nodes_.parent_partitions_[curr_index_];
    }

    /// Templated to allow comparison between const and non-const types
    template<typename OtherAccessorType>
    inline bool operator==(const OtherAccessorType& other) const
    {
        return &nodes_ == &other.nodes_ && curr_index_ == other.curr_index_;
    }

    /// Templated to allow comparison between const and non-const types
    template<typename OtherAccessorType>
    inline bool operator!=(const OtherAccessorType& other) const
    {
        return !(*this == other);
    }

    /// @brief Serializes a single node per the interface.  Note: if serializing a mesh, serialize
    /// the NodeContainer instead of the accessor.
    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("global_index", global_index()),
          utility::named_value("x", x()),
          utility::named_value("y", y()),
          utility::named_value("z", z()));
        if (nodes_.store_parent_partitions_)
        {
            s(utility::named_value("parent_partition", parent_partition()));
        }
    }

    // Note: omitted deserialize because an accessor is not stored (it only provides a view)

    /// Allow friend access to other NodeAccessor types for comparison operators
    template<typename NC>
    friend class NodeAccessorBase;
    /// Allow to iterator accessor to next() method
    template<typename NodeAccessorType>
    friend class NodeIteratorBase;
    /// Allow to iterator accessor to next() and go_to() methods
    template<typename NodeAccessorType>
    friend class NodeSubrangeIteratorBase;

  private:
    NC& nodes_;
    LocalIndex curr_index_;

    /// Moves accessor to the next node in the container (based on local index).
    inline void next()
    {
        curr_index_++;
    }
    /// Moves accessor to the given node (does not check if index is valid).
    inline void go_to(LocalIndex i)
    {
        curr_index_ = i;
    }
};

using NodeAccessor      = NodeAccessorBase<NodeContainer>;
using NodeConstAccessor = NodeAccessorBase<const NodeContainer>;
}