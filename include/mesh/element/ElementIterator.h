#pragma once

#include "mesh/element/ElementAccessor.h"

namespace mesh
{
/// @brief Iterator that sequentially moves through an ElementContainer.
/// @tparam ElementAccessorType - const or nonconst accessor determines if iterator is const or
/// noncost
template<typename ElementAccessorType>
class ElementIteratorBase
{
  public:
    // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = ElementAccessorType;
    using pointer           = const ElementAccessorType*;
    using reference         = const ElementAccessorType&;
    using iterator_category = std::forward_iterator_tag;

    using IteratorType         = ElementIteratorBase<ElementAccessorType>;
    using NodeContainerType    = typename ElementAccessorType::NodeContainerType;
    using ElementContainerType = typename ElementAccessorType::ElementContainerType;

    ElementIteratorBase(NodeContainerType& nodes,
                        ElementContainerType& elements,
                        LocalIndex starting_index = 0)
      : accessor_(nodes, elements, starting_index)
    {
    }

    ~ElementIteratorBase() {}

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

    bool operator==(const ElementIteratorBase<ElementAccessor>& other) const
    {
        return accessor_ == *other;
    }
    bool operator==(const ElementIteratorBase<ElementConstAccessor>& other) const
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

  protected:
    ElementAccessorType accessor_;
};

using ElementIterator      = ElementIteratorBase<ElementAccessor>;
using ElementConstIterator = ElementIteratorBase<ElementConstAccessor>;
}