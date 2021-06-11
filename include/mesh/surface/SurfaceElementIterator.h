#pragma once

#include "mesh/surface/SurfaceElementAccessor.h"

namespace mesh
{
/// @brief Iterator that sequentially moves through an SurfaceMesh.
/// @tparam SurfaceElementAccessorType - const or nonconst accessor determines if iterator is const
/// or noncost
template<typename SurfaceElementAccessorType>
class SurfaceElementIteratorBase
{
  public:
    // Iterator traits
    using difference_type   = LocalIndex;
    using value_type        = SurfaceElementAccessorType;
    using pointer           = SurfaceElementAccessorType*;
    using reference         = SurfaceElementAccessorType&;
    using const_reference   = const SurfaceElementAccessorType&;
    using iterator_category = std::forward_iterator_tag;

    using IteratorType = SurfaceElementIteratorBase<SurfaceElementAccessorType>;
    using MeshType     = typename SurfaceElementAccessorType::MeshType;

    SurfaceElementIteratorBase(MeshType& surf_mesh, LocalIndex starting_index = 0)
      : accessor_(surf_mesh, starting_index)
    {
    }

    ~SurfaceElementIteratorBase() {}

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

    bool operator==(const SurfaceElementIteratorBase<SurfaceElementAccessorType>& other) const
    {
        return accessor_ == *other;
    }

    bool operator!=(const SurfaceElementIteratorBase<SurfaceElementAccessorType>& other) const
    {
        return !(*this == other);
    }

    reference operator*()
    {
        return accessor_;
    }

    const_reference operator*() const
    {
        return accessor_;
    }

  protected:
    SurfaceElementAccessorType accessor_;
};

using SurfaceElementIterator      = SurfaceElementIteratorBase<SurfaceElementAccessor>;
}