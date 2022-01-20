#pragma once

#include "mesh/node/NodeAccessor.h"

namespace mesh
{
/// @brief Iterator that sequentially moves through a NodeContainer.
/// @tparam NodeAccessorType - const or nonconst accessor determines if iterator is const or noncost
template<typename NodeAccessorType>
class NodeIteratorBase
{
  public:
    // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = NodeAccessorType;
    using pointer           = const NodeAccessorType*;
    using reference         = const NodeAccessorType&;
    using iterator_category = std::forward_iterator_tag;

    using IteratorType      = NodeIteratorBase<NodeAccessorType>;
    using NodeContainerType = typename NodeAccessorType::NodeContainerType;

    NodeIteratorBase(NodeContainerType& nodes, LocalIndex starting_index = 0)
      : accessor_(nodes, starting_index)
    {
    }

    ~NodeIteratorBase() {}

    /// @brief Pre-increment operator
    IteratorType& operator++()
    {
        accessor_.next();
        return *this;
    }

    /// @brief Post-increment operator
    IteratorType operator++(int)
    {
        IteratorType retval = *this;
        ++(*this);
        return retval;
    }

    bool operator==(const NodeIteratorBase<NodeAccessor>& other) const
    {
        return accessor_ == *other;
    }
    bool operator==(const NodeIteratorBase<NodeConstAccessor>& other) const
    {
        return accessor_ == *other;
    }

    /// Templated to allow comparison between const and non-const types
    template<typename OtherIteratorType>
    bool operator!=(const OtherIteratorType& other) const
    {
        return !(*this == other);
    }

    reference operator*() const
    {
        return accessor_;
    }

    /// @brief Returns true if the iterator is at the last node in the container already
    inline bool is_end() const
    {
        return accessor_.is_end();
    }

    /// @brief Non-const access to underlying NodeAccessor (useful for Python)
    value_type get()
    {
        return accessor_;
    }

  protected:
    NodeAccessorType accessor_;
};

using NodeIterator      = NodeIteratorBase<NodeAccessor>;
using NodeConstIterator = NodeIteratorBase<NodeConstAccessor>;
}