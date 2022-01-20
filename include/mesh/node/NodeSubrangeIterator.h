#pragma once

#include "mesh/node/NodeAccessor.h"

namespace mesh
{
/// @brief Provides an iterator over a subset of a NodeContainer instance.
/// @tparam NodeAccessorType - const or nonconst accessor determines if iterator is const or noncost
template<typename NodeAccessorType>
class NodeSubrangeIteratorBase
{
  public:
      // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = NodeAccessorType;
    using pointer           = const NodeAccessorType*;
    using reference         = const NodeAccessorType&;
    using iterator_category = std::forward_iterator_tag;

    using IteratorType      = NodeSubrangeIteratorBase<NodeAccessorType>;
    using NodeContainerType = typename NodeAccessorType::NodeContainerType;

    /// Constructor for iteration over an arbitrary list of nodes
    /// @param nodes - the NodeContainer for which the nodes to iterate over is a subset
    /// @param nodes_to_visit - list of local indices of nodes to iterate over (can be unsorted)
    /// @param starting_subrange_index - sets the iterator to the given index of nodes_to_visit.  If
    /// nodes_to_visit == nodes_to_visit.size() results in the end iterator for the range.
    NodeSubrangeIteratorBase(NodeContainerType& nodes,
                             const std::vector<LocalIndex>& nodes_to_visit,
                             LocalIndex starting_subrange_index = 0)
      : NodeSubrangeIteratorBase(nodes,
                                 nodes_to_visit.data(),
                                 static_cast<LocalIndex>(nodes_to_visit.size()),
                                 starting_subrange_index)
    {
    }

    /// Constructor for iteration over an element's connectivity
    NodeSubrangeIteratorBase(NodeContainerType& nodes,
                             const LocalIndex* nodes_to_visit,
                             LocalIndex n_nodes_to_visit,
                             LocalIndex starting_subrange_index = 0)
      : accessor_(nodes)
    {
        nodes_to_visit_.resize(n_nodes_to_visit + 1, 0);
        for (LocalIndex i = 0; i < n_nodes_to_visit; i++)
        {
            nodes_to_visit_[i] = nodes_to_visit[i];
        }

        curr_subrange_index_ = starting_subrange_index;
        accessor_.go_to(nodes_to_visit_[curr_subrange_index_]);
    }

    /// @brief Moves iterator to next node in the subset to visit
    IteratorType& operator++()
    {
        curr_subrange_index_++;
        accessor_.go_to(nodes_to_visit_[curr_subrange_index_]);
        return *this;
    }

    /// @brief Post-increment operator
    IteratorType operator++(int)
    {
        IteratorType retval = *this;
        ++(*this);
        return retval;
    }

    bool operator==(const NodeSubrangeIteratorBase<NodeAccessor>& other) const
    {
        return nodes_to_visit_ == other.nodes_to_visit_ &&
               curr_subrange_index_ == other.curr_subrange_index_;
    }
    bool operator==(const NodeSubrangeIteratorBase<NodeConstAccessor>& other) const
    {
        return nodes_to_visit_ == other.nodes_to_visit_ &&
               curr_subrange_index_ == other.curr_subrange_index_;
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
    
    // Iterator like methods

    /// @brief Moves accessor to the next node in the subrange
    inline void next()
    {
        (*this)++;
    }

    /// @brief Moves accessor to the given node (does not check if index is valid).
    inline void go_to(LocalIndex i)
    {
        curr_index_ = i;
    }

  protected:
    NodeAccessorType accessor_;
    /// Includes a 0 at the end so that operator++ does not cause seg fault without need for an
    /// if statement
    std::vector<LocalIndex> nodes_to_visit_;
    LocalIndex curr_subrange_index_;
};

using NodeSubrangeIterator      = NodeSubrangeIteratorBase<NodeAccessor>;
using NodeSubrangeConstIterator = NodeSubrangeIteratorBase<NodeConstAccessor>;

}