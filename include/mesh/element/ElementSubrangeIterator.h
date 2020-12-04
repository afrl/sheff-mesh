#pragma once 

#include "mesh/element/ElementIterator.h"

namespace mesh
{
/// @brief Provides an iterator over a subset of an ElementContainer instance.
/// @tparam ElementAccessorType - const or nonconst accessor determines if iterator is const or
/// noncost
template<typename ElementAccessorType>
class ElementSubrangeIteratorBase
{
  public:
    // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = ElementAccessorType;
    using pointer           = const ElementAccessorType*;
    using reference         = const ElementAccessorType&;
    using iterator_category = std::forward_iterator_tag;

    using IteratorType         = ElementSubrangeIteratorBase<ElementAccessorType>;
    using NodeContainerType    = typename ElementAccessorType::NodeContainerType;
    using ElementContainerType = typename ElementAccessorType::ElementContainerType;

    /// Constructor for iteration over an arbitrary list of elements
    /// @param nodes - the NodeContainer that the "elements" parameter references
    /// @param elements - the ElementContainer for which the elements to iterate over is a subset
    /// @param elements_to_visit_ - list of local indices of elements to iterate over (can be unsorted)
    /// @param starting_subrange_index - sets the iterator to the given index of elements_to_visit_.
    /// If elements_to_visit_ == elements_to_visit_.size() results in the end iterator for the
    /// range.
    ElementSubrangeIteratorBase(NodeContainerType& nodes,
                                ElementContainerType& elements,
                                const std::vector<LocalIndex>& elements_to_visit_,
                                LocalIndex starting_subrange_index = 0)
      : ElementSubrangeIteratorBase(nodes,
                                    elements,
                                    elements_to_visit_.data(),
                                    static_cast<LocalIndex>(elements_to_visit_.size()),
                                    starting_subrange_index)
    {
    }

    /// Constructor for iteration over an element's connectivity
    ElementSubrangeIteratorBase(NodeContainerType& nodes,
                                ElementContainerType& elements,
                                const LocalIndex* elements_to_visit,
                                LocalIndex n_elements_to_visit,
                                LocalIndex starting_subrange_index = 0)
      : accessor_(nodes, elements)
    {
        elements_to_visit_.resize(n_elements_to_visit + 1, 0);
        for (LocalIndex i = 0; i < n_elements_to_visit; i++)
        {
            elements_to_visit_[i] = elements_to_visit[i];
        }

        curr_subrange_index_ = starting_subrange_index;
        accessor_.go_to(elements_to_visit_[curr_subrange_index_]);
    }

    ~ElementSubrangeIteratorBase() {}

    /// Moves iterator to next element in the subset to visit
    IteratorType& operator++()
    {
        curr_subrange_index_++;
        accessor_.go_to(elements_to_visit_[curr_subrange_index_]);
        return *this;
    }

    /// @brief Post-increment operator
    IteratorType operator++(int)
    {
        IteratorType retval = *this;
        ++(*this);
        return retval;
    }

    bool operator==(const ElementSubrangeIteratorBase<ElementAccessor>& other) const
    {
        return elements_to_visit_ == other.elements_to_visit_ &&
               curr_subrange_index_ == other.curr_subrange_index_;
    }
    bool operator==(const ElementSubrangeIteratorBase<ElementConstAccessor>& other) const
    {
        return elements_to_visit_ == other.elements_to_visit_ &&
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

    /// @brief Returns the index of where the iterator is within the list of elements
    /// to visit (given to the constructor).
    LocalIndex subrange_index() const
    {
        return curr_subrange_index_;
    }

  protected:
    ElementAccessorType accessor_;
    /// Includes a 0 at the end so that operator++ does not cause seg fault without need for an
    /// if statement
    std::vector<LocalIndex> elements_to_visit_;
    LocalIndex curr_subrange_index_;
};

using ElementSubrangeIterator      = ElementSubrangeIteratorBase<ElementAccessor>;
using ElementSubrangeConstIterator = ElementSubrangeIteratorBase<ElementConstAccessor>;
}