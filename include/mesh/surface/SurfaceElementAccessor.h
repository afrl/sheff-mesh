#pragma once

#include "mesh/common/IteratorRange.h"
#include "mesh/surface/MeshSurface.h"

#include <assert.h>

namespace mesh
{
/// @brief Provides a const interface to access surface element information.
class SurfaceElementAccessor
{
  public:
    // typedef for the benefit of SurfaceElementIterator
    using MeshType = const MeshSurface;

    SurfaceElementAccessor(const MeshSurface& surf_mesh, LocalIndex starting_index = 0);
    SurfaceElementAccessor(const SurfaceElementAccessor& other);

    /// Returns the index of the current element within the ElementContainer that contains it.
    const LocalIndex& local_index() const;

    /// Returns the type of element (in terms of topology), see ElementType
    const ElementType& type() const;

    /// Returns the ordered local indices of nodes that form the connectivity of the current
    /// element.  Note: the (const) slice provides an immutable view of the data.
    utility::Slice<const LocalIndex, ConnIndex> connectivity() const;

    /// Allows access to a particular node in the element (node accessor inherits const-ness of the
    /// element accessor).
    NodeConstAccessor node(ConnIndex i) const;

    /// Returns a range over the ordered nodes that compose the connectivity of the element.  The
    /// order will match the indices returned by connectivity().
    IteratorRange<NodeSubrangeConstIterator> nodes() const;

    bool operator==(const SurfaceElementAccessor& other) const;

    bool operator!=(const SurfaceElementAccessor& other) const;

    /// Allow to iterator accessor to next() and goto() methods
    template<typename SurfaceElementAccessorType>
    friend class SurfaceElementIteratorBase;

  private:
    const MeshSurface& surf_mesh_;
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
}