#pragma once

#include "mesh/face/FaceAccessor.h"

namespace mesh
{
/// @brief Iterator that sequentially moves through an ElementContainer.
/// @tparam ElementAccessorType - const or nonconst accessor determines if iterator is const or
/// noncost
class FaceConstIterator
{
  public:
    // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = FaceConstAccessor;
    using pointer           = const FaceConstAccessor*;
    using reference         = const FaceConstAccessor&;
    using iterator_category = std::forward_iterator_tag;

    FaceConstIterator(const FaceContainer& faces, ConnIndex starting_index = 0)
      : accessor_(faces, starting_index)
    {
    }

    ~FaceConstIterator() {}

    /// @brief Pre-increment operator
    FaceConstIterator& operator++()
    {
        accessor_.next();
        return *this;
    }

    /// @brief Post-increment operator
    FaceConstIterator operator++(int)
    {
        FaceConstIterator retval = *this;
        ++(*this);
        return retval;
    }

    bool operator==(const FaceConstIterator& other) const
    {
        return accessor_ == *other;
    }

    bool operator!=(const FaceConstIterator& other) const
    {
        return !(*this == other);
    }

    reference operator*() const
    {
        return accessor_;
    }

  protected:
    FaceConstAccessor accessor_;
};
}