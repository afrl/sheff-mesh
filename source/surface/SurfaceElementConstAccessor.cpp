#include "mesh/surface/SurfaceElementConstAccessor.h"

namespace mesh
{
SurfaceElementConstAccessor::SurfaceElementConstAccessor(const MeshSurface& surf_mesh,
                                                         LocalIndex starting_index)
  : surf_mesh_(surf_mesh)
  , curr_index_(starting_index)
{
}

SurfaceElementConstAccessor::SurfaceElementConstAccessor(const SurfaceElementConstAccessor& other)
  : surf_mesh_(other.surf_mesh_)
  , curr_index_(other.curr_index_)
{
}

const LocalIndex& SurfaceElementConstAccessor::local_index() const
{
    return curr_index_;
}

const ElementType& SurfaceElementConstAccessor::type() const
{
    return surf_mesh_.types_[curr_index_];
}

utility::Slice<const LocalIndex> SurfaceElementConstAccessor::connectivity() const
{
    const auto& conn_index_start = surf_mesh_.connectivity_offsets_[curr_index_];
    const auto& conn_size = surf_mesh_.connectivity_offsets_[curr_index_ + 1] - conn_index_start;
    return utility::Slice<const LocalIndex>(&surf_mesh_.connectivities_[conn_index_start],
                                            conn_size);
}

NodeConstAccessor SurfaceElementConstAccessor::node(ConnIndex i) const
{
    const auto& node_local_index =
      surf_mesh_.connectivities_[surf_mesh_.connectivity_offsets_[curr_index_] + i];
    return NodeConstAccessor(surf_mesh_.volume_mesh_.nodes(), node_local_index);
}

IteratorRange<NodeSubrangeConstIterator> SurfaceElementConstAccessor::nodes() const
{
    const auto& conn_index_start = surf_mesh_.connectivity_offsets_[curr_index_];
    const auto& conn_size =
      static_cast<LocalIndex>(surf_mesh_.connectivity_offsets_[curr_index_ + 1] - conn_index_start);
    return {
        NodeSubrangeConstIterator(surf_mesh_.volume_mesh_.nodes(),
                                  &surf_mesh_.connectivities_[conn_index_start],
                                  conn_size), // begin
        NodeSubrangeConstIterator(surf_mesh_.volume_mesh_.nodes(),
                                  &surf_mesh_.connectivities_[conn_index_start],
                                  conn_size,
                                  conn_size) // end
    };
}

bool SurfaceElementConstAccessor::operator==(const SurfaceElementConstAccessor& other) const
{
    return &surf_mesh_ == &other.surf_mesh_ && curr_index_ == other.curr_index_;
}

bool SurfaceElementConstAccessor::operator!=(const SurfaceElementConstAccessor& other) const
{
    return !(*this == other);
}
}