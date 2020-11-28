#include "mesh/surface/SurfaceMesh.h"

#include "mesh/surface/SurfaceElementAccessor.h"
#include "mesh/surface/SurfaceElementIterator.h"

namespace mesh
{

SurfaceMesh::SurfaceMesh(Mesh& volume_mesh)
  : volume_mesh_(volume_mesh)
  , connectivity_offsets_({ 0 })
{
}

SurfaceMesh::SurfaceMesh(SurfaceMesh&& other) noexcept
  : volume_mesh_(other.volume_mesh_)
{
    types_.swap(other.types_);
    connectivities_.swap(other.connectivities_);
    connectivity_offsets_.swap(other.connectivity_offsets_);
}

SurfaceElementIterator SurfaceMesh::begin()
{
    return SurfaceElementIterator(*this, 0);
}

SurfaceElementIterator SurfaceMesh::end()
{
    return SurfaceElementIterator(*this, size());
}

inline SurfaceElementAccessor SurfaceMesh::operator[](LocalIndex i)
{
    return SurfaceElementAccessor(*this, i);
}

SurfaceElementAccessor SurfaceMesh::back()
{
    return SurfaceElementAccessor(*this, size() - 1);
}

LocalIndex SurfaceMesh::size() const
{
    return static_cast<LocalIndex>(types_.size());
}

Mesh SurfaceMesh::mesh()
{
    Mesh mesh;

    LocalIndex count = 0;
    std::vector<LocalIndex> mesh_to_surf_node_map(volume_mesh_.nodes().size(),
                                                  std::numeric_limits<LocalIndex>::max());
    for (auto& surf_element : *this)
    {
        for (auto& i : surf_element.connectivity())
        {
            if (mesh_to_surf_node_map[i] == std::numeric_limits<LocalIndex>::max())
            {
                mesh.nodes().add(volume_mesh_.nodes()[i].point());
                mesh_to_surf_node_map[i] = count++;
            }
        }
    }

    for (auto& surf_element : *this)
    {
        auto conn = std::vector<LocalIndex>(surf_element.connectivity().begin(),
                                            surf_element.connectivity().end());
        for (auto& i : conn)
        {
            i = mesh_to_surf_node_map[i];
        }
        mesh.elements().add(surf_element.type(), conn);
    }

    return mesh;
}

void SurfaceMesh::add(const FaceConstAccessor& face)
{
    types_.push_back(element_equivalent_for_face(face.type()));
    for (auto& i : face.connectivity())
    {
        connectivities_.push_back(i);
    }
    connectivity_offsets_.push_back(connectivities_.size());
}
}