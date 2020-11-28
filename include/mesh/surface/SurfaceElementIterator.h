#pragma once

#include "mesh/surface/SurfaceElementAccessor.h"

namespace mesh
{
/// @brief Iterator that sequentially moves through an SurfaceMesh.
class SurfaceElementIterator
{
  public:
    // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = SurfaceElementAccessor;
    using pointer           = const SurfaceElementAccessor*;
    using reference         = const SurfaceElementAccessor&;
    using iterator_category = std::forward_iterator_tag;

    using IteratorType         = SurfaceElementIterator;

    SurfaceElementIterator(SurfaceMesh& surf_mesh,
                        LocalIndex starting_index = 0)
      : accessor_(surf_mesh, starting_index)
    {
    }

    ~SurfaceElementIterator() {}

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

    bool operator==(const SurfaceElementIterator& other) const
    {
        return accessor_ == *other;
    }

    bool operator!=(const SurfaceElementIterator& other) const
    {
        return !(*this == other);
    }

    reference operator*() const
    {
        return accessor_;
    }

  protected:
    SurfaceElementAccessor accessor_;
};
}