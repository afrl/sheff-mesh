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
    using GlobalIndexType =
      typename std::conditional<std::is_const<NC>::value, const GlobalIndex, GlobalIndex>::type;
    using FloatTypeType =
      typename std::conditional<std::is_const<NC>::value, const FloatType, FloatType>::type;
    using PointType = typename std::conditional<std::is_const<NC>::value, const Point, Point>::type;
    using PartitionIndexType = 
      typename std::conditional<std::is_const<NC>::value, const PartitionIndex, PartitionIndex>::type;
    using MapCoordType = typename std::
      conditional<std::is_const<NC>::value, const Eigen::Vector3d, Eigen::Vector3d>::type;

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
    inline GlobalIndexType& global_index() const
    {
        return nodes_.global_indices_[curr_index_];
    }

    /// Returns the first component of the coordinate (in the global coordinate system)
    inline FloatTypeType& x() const
    {
        return nodes_.points_[curr_index_][0];
    }

    /// Returns the second component of the coordinate (in the global coordinate system)
    inline FloatTypeType& y() const
    {
        return nodes_.points_[curr_index_][1];
    }

    /// Returns the third component of the coordinate (in the global coordinate system)
    inline FloatTypeType& z() const
    {
        return nodes_.points_[curr_index_][2];
    }

    /// Returns point of node
    inline PointType& point() const
    {
        return nodes_.points_[curr_index_];
    }

    /// Returns a component of the point by reference
    inline FloatTypeType& point(int8_t coordinate_axis) const
    {
        return nodes_.points_[curr_index_][coordinate_axis];
    }

    /// Returns a map of the coordinates to a 3D Eigen vector. It can be used as a normal Eigen
    /// vector.
    inline Eigen::Map<MapCoordType> coord() const
    {
        return Eigen::Map<MapCoordType>(&nodes_.points_[curr_index_][0], 3);
    }

    /// Returns the parent partition of the current node
    inline PartitionIndexType& parent_partition() const
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

    std::string str() const
    {
        return std::string(utility::demangle(typeid(*this).name())) + ": { " +
               std::to_string(curr_index_) + " }";
    }

    /// Allow friend access to other NodeAccessor types for comparison operators
    template<typename NC2>
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

    // Methods to help the implementation of iterators

    /// @brief Moves accessor to the next node in the container (based on local index).
    inline void next()
    {
        curr_index_++;
    }

    /// @brief Moves accessor to the given node (does not check if index is valid).
    inline void go_to(LocalIndex i)
    {
        curr_index_ = i;
    }

    /// @brief Returns true if the accessor is at the last node in the container already
    inline bool is_end() const
    {
        return curr_index_ == nodes_.size();
    }
};

using NodeAccessor      = NodeAccessorBase<NodeContainer>;
using NodeConstAccessor = NodeAccessorBase<const NodeContainer>;
}
