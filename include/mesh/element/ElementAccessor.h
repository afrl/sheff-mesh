#pragma once

#include "mesh/common/IteratorRange.h"
#include "mesh/element/ElementContainer.h"
#include "mesh/node/NodeSubrangeIterator.h"

#include <assert.h>

namespace mesh
{
/// @brief Provides an interface to access element information.
/// Templated to allow a const and non-const version of the accessor without duplicate code.
template<typename NC, typename EC>
class ElementAccessorBase
{
  public:
    using NodeContainerType    = NC;
    using ElementContainerType = EC;
    using ConnectivityType =
      typename std::conditional<std::is_const<NC>::value, const LocalIndex, LocalIndex>::type;

    ElementAccessorBase(NC& nodes, EC& elements, LocalIndex starting_index = 0)
      : nodes_(nodes)
      , elements_(elements)
      , curr_index_(starting_index)
    {
    }

    /// Copy constructor takes another accessor type, allowing:
    /// 1) ElementAccessor to be created from another ElementAccessor
    /// 2) ElementConstAccessor to be created from another ElementAccessor
    /// 3) ElementConstAccessor to be created from another ElementConstAccessor
    /// Note: if you try to create a ElementAccessor from a ElementConstAccessor, then a compiler
    /// error will occur about qualifiers and const-ness
    template<typename ElementAccessorType>
    ElementAccessorBase(const ElementAccessorType& other)
      : nodes_(other.nodes_)
      , elements_(other.elements_)
      , curr_index_(other.curr_index_)
    {
    }

    /// Returns the index of the current element within the ElementContainer that contains it.
    inline const LocalIndex& local_index() const
    {
        return curr_index_;
    }

    /// Returns the globally unique index for the current element.
    inline auto& global_index() const
    {
        return elements_.global_indices_[curr_index_];
    }

    /// Returns the type of element (in terms of topology), see ElementType
    inline auto& type() const
    {
        return elements_.types_[curr_index_];
    }

    /// Returns the ordered local indices of nodes that form the connectivity of the current
    /// element.  Note: the (const) slice provides an immutable view of the data.
    inline utility::Slice<ConnectivityType> connectivity() const
    {
        const auto& conn_index_start = elements_.connectivity_offsets_[curr_index_];
        const auto& conn_size = elements_.connectivity_offsets_[curr_index_ + 1] - conn_index_start;
        return utility::Slice<ConnectivityType>(&elements_.connectivities_[conn_index_start],
                                                conn_size);
    }

    /// Allows access to a particular node in the element (node accessor inherits const-ness of the
    /// element accessor).
    inline NodeAccessorBase<NC> node(ConnIndex i) const
    {
        const auto& node_local_index =
          elements_.connectivities_[elements_.connectivity_offsets_[curr_index_] + i];
        return NodeAccessorBase<NC>(nodes_, node_local_index);
    }

    /// Returns a range over the ordered nodes that compose the connectivity of the element.  The
    /// order will match the indices returned by connectivity().
    /// Note: iterators inherit const-ness of the element accessor
    IteratorRange<NodeSubrangeIteratorBase<NodeAccessorBase<NC>>> nodes() const
    {
        const auto& conn_index_start = elements_.connectivity_offsets_[curr_index_];
        const auto& conn_size        = static_cast<LocalIndex>(
          elements_.connectivity_offsets_[curr_index_ + 1] - conn_index_start);
        return {
            NodeSubrangeIteratorBase<NodeAccessorBase<NC>>(
              nodes_,
              &elements_.connectivities_[conn_index_start],
              conn_size), // begin
            NodeSubrangeIteratorBase<NodeAccessorBase<NC>>(
              nodes_,
              &elements_.connectivities_[conn_index_start],
              conn_size,
              conn_size) // end
        };
    }

    /// Returns the partition that owns this element.
    inline auto& parent_partition() const
    {
        assert(!elements_.parent_partitions_.empty());
        return elements_.parent_partitions_[curr_index_];
    }

    /// Returns the object that contains the local coordinate systems throughout the elemetn.
    inline auto& local_coord_sys() const
    {
        assert(!elements_.lcs_.empty());
        return elements_.lcs_[curr_index_];
    }

    /// Returns a reference to the underlying NodeContainer that this accessor refers to.
    inline auto& node_container() const
    {
        return nodes_;
    }

    /// Templated to allow comparison between const and non-const types
    template<typename OtherAccessorType>
    inline bool operator==(const OtherAccessorType& other) const
    {
        return &nodes_ == &other.nodes_ && &elements_ == &other.elements_ &&
               curr_index_ == other.curr_index_;
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
        s(utility::named_value("global_index", global_index()));
    }

    // Note: omitted deserialize because an accessor is not stored (it only provides a view)

    /// Allow friend access to other ElementAccessor types for comparison operators
    template<typename NC, typename EC>
    friend class ElementAccessorBase;
    /// Allow to iterator accessor to next() and goto() methods
    template<typename ElementAccessorType>
    friend class ElementIteratorBase;
    template<typename ElementAccessorType>
    friend class ElementSubrangeIteratorBase;

  private:
    NC& nodes_;
    EC& elements_;
    LocalIndex curr_index_;

    /// Moves accessor to the next element in the container (based on local index).
    inline void next()
    {
        curr_index_++;
    }
    /// Moves accessor to the given element (does not check if index is valid).
    inline void go_to(LocalIndex i)
    {
        curr_index_ = i;
    }
};

using ElementAccessor      = ElementAccessorBase<NodeContainer, ElementContainer>;
using ElementConstAccessor = ElementAccessorBase<const NodeContainer, const ElementContainer>;
}