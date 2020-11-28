#pragma once

#include "mesh/common/IteratorRange.h"
#include "mesh/surface/SurfaceMesh.h"

#include <assert.h>

namespace mesh
{
/// @brief Provides an interface to access surface element information.
class SurfaceElementAccessor
{
  public:

    SurfaceElementAccessor(SurfaceMesh& surf_mesh, LocalIndex starting_index = 0)
      : surf_mesh_(surf_mesh)
      , curr_index_(starting_index)
    {
    }

    SurfaceElementAccessor(const SurfaceElementAccessor& other)
      : surf_mesh_(other.surf_mesh_)
      , curr_index_(other.curr_index_)
    {
    }

    /// Returns the index of the current element within the ElementContainer that contains it.
    inline const LocalIndex& local_index() const
    {
        return curr_index_;
    }

    /// Returns the type of element (in terms of topology), see ElementType
    inline auto& type() const
    {
        return surf_mesh_.types_[curr_index_];
    }

    /// Returns the ordered local indices of nodes that form the connectivity of the current
    /// element.  Note: the (const) slice provides an immutable view of the data.
    inline utility::Slice<const LocalIndex> connectivity() const
    {
        const auto& conn_index_start = surf_mesh_.connectivity_offsets_[curr_index_];
        const auto& conn_size =
          surf_mesh_.connectivity_offsets_[curr_index_ + 1] - conn_index_start;
        return utility::Slice<const LocalIndex>(&surf_mesh_.connectivities_[conn_index_start],
                                                conn_size);
    }

    /// Allows access to a particular node in the element (node accessor inherits const-ness of the
    /// element accessor).
    inline NodeAccessor node(ConnIndex i)
    {
        const auto& node_local_index =
          surf_mesh_.connectivities_[surf_mesh_.connectivity_offsets_[curr_index_] + i];
        return NodeAccessor(surf_mesh_.volume_mesh_.nodes(), node_local_index);
    }

    /// Returns a range over the ordered nodes that compose the connectivity of the element.  The
    /// order will match the indices returned by connectivity().
    IteratorRange<NodeSubrangeIterator> nodes()
    {
        const auto& conn_index_start = surf_mesh_.connectivity_offsets_[curr_index_];
        const auto& conn_size        = static_cast<LocalIndex>(
          surf_mesh_.connectivity_offsets_[curr_index_ + 1] - conn_index_start);
        return {
            NodeSubrangeIterator(surf_mesh_.volume_mesh_.nodes(),
                                 &surf_mesh_.connectivities_[conn_index_start],
                                 conn_size), // begin
            NodeSubrangeIterator(surf_mesh_.volume_mesh_.nodes(),
                                 &surf_mesh_.connectivities_[conn_index_start],
                                 conn_size,
                                 conn_size) // end
        };
    }

    /// Returns a reference to the underlying NodeContainer that this accessor refers to.
    inline auto& node_container() const
    {
        return surf_mesh_.volume_mesh_.nodes();
    }

    inline bool operator==(const SurfaceElementAccessor& other) const
    {
        return &surf_mesh_ == &other.surf_mesh_ && curr_index_ == other.curr_index_;
    }

    inline bool operator!=(const SurfaceElementAccessor& other) const
    {
        return !(*this == other);
    }

    /// Allow to iterator accessor to next() and goto() methods
    friend class SurfaceElementIterator;

  private:
    SurfaceMesh& surf_mesh_;
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