#include "mesh/surface/MeshSurface.h"

#include "mesh/surface/SurfaceElementAccessor.h"
#include "mesh/surface/SurfaceElementIterator.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

#include <unordered_map>

namespace mesh
{

MeshSurface::MeshSurface(const Mesh& volume_mesh)
  : volume_mesh_(volume_mesh)
  , connectivity_offsets_({ 0 })
{
}

MeshSurface::MeshSurface(MeshSurface&& other) noexcept
  : volume_mesh_(other.volume_mesh_)
{
    types_.swap(other.types_);
    connectivities_.swap(other.connectivities_);
    connectivity_offsets_.swap(other.connectivity_offsets_);
}

MeshSurface::SurfaceElementIterator MeshSurface::begin() const
{
    return SurfaceElementIterator(*this, 0);
}

MeshSurface::SurfaceElementIterator MeshSurface::end() const
{
    return SurfaceElementIterator(*this, size());
}

SurfaceElementAccessor MeshSurface::operator[](LocalIndex i) const
{
    return SurfaceElementAccessor(*this, i);
}

SurfaceElementAccessor MeshSurface::back() const
{
    return SurfaceElementAccessor(*this, size() - 1);
}

LocalIndex MeshSurface::size() const
{
    return static_cast<LocalIndex>(types_.size());
}

Mesh MeshSurface::mesh() const
{
    return mesh_with_map().surface_mesh;
}

SurfaceMeshWithMaps MeshSurface::mesh_with_map() const
{
    Mesh mesh;
    mesh.n_partitions() = volume_mesh_.n_partitions();
    mesh.partition_id() = volume_mesh_.partition_id();

    LocalIndex count = 0;
    std::vector<LocalIndex> vol_to_surf_node_map(volume_mesh_.nodes().size(),
                                                  std::numeric_limits<LocalIndex>::max());
    std::vector<LocalIndex> surf_to_vol_node_map;
    for (auto& surf_element : *this)
    {
        for (auto& i : surf_element.connectivity())
        {
            if (vol_to_surf_node_map[i] == std::numeric_limits<LocalIndex>::max())
            {
                auto node = volume_mesh_.nodes()[i];
                mesh.nodes().add(node.point(), node.global_index());
                vol_to_surf_node_map[i] = count++;
                surf_to_vol_node_map.push_back(i);
            }
        }
    }

    for (auto& surf_element : *this)
    {
        auto conn = std::vector<LocalIndex>(surf_element.connectivity().begin(),
                                            surf_element.connectivity().end());
        for (auto& i : conn)
        {
            i = vol_to_surf_node_map[i];
        }
        mesh.elements().add(surf_element.type(), conn);
    }

    return { mesh, vol_to_surf_node_map, surf_to_vol_node_map };
}

void MeshSurface::add(const FaceConstAccessor& face)
{
    types_.push_back(element_equivalent_for_face(face.type()));
    auto flip_permute = flip_conn_permute(face.type());
    auto conn         = face.connectivity();
    for (auto& i : flip_permute)
    {
        connectivities_.push_back(conn[i]);
    }
    connectivity_offsets_.push_back(connectivities_.size());
}

void MeshSurface::add(const ElementType& type, const utility::Slice<const LocalIndex>& conn)
{
    types_.push_back(type);
    for (auto& i : conn)
    {
        connectivities_.push_back(i);
    }
    connectivity_offsets_.push_back(connectivities_.size());
}

void MeshSurface::add(const ElementType& type, const utility::Slice<const LocalIndex, ConnIndex>& conn)
{
    types_.push_back(type);
    for (auto& i : conn)
    {
        connectivities_.push_back(i);
    }
    connectivity_offsets_.push_back(connectivities_.size());
}

MeshSurface MeshSurface::operator-(const MeshSurface& other) const
{
    std::unordered_map<std::size_t, LocalIndex> hash_to_this_surf_elem_map;
    for (auto& surf_elem : *this)
    {
        auto hash                        = utility::hash_sorted(surf_elem.connectivity());
        hash_to_this_surf_elem_map[hash] = surf_elem.local_index();
    }

    for (auto& surf_elem : other)
    {
        auto hash      = utility::hash_sorted(surf_elem.connectivity());
        auto find_iter = hash_to_this_surf_elem_map.find(hash);
        if (find_iter != hash_to_this_surf_elem_map.end())
        {
            hash_to_this_surf_elem_map.erase(find_iter);
        }
    }

    MeshSurface diff(volume_mesh_);
    for (auto& p : hash_to_this_surf_elem_map)
    {
        auto surf_elem_to_copy = (*this)[p.second];
        diff.add(surf_elem_to_copy.type(), surf_elem_to_copy.connectivity());
    }

    return diff;
}

MeshSurface& MeshSurface::operator+=(const MeshSurface& other)
{
    std::unordered_map<std::size_t, LocalIndex> hash_to_this_surf_elem_map;
    for (auto& surf_elem : *this)
    {
        auto hash                        = utility::hash_sorted(surf_elem.connectivity());
        hash_to_this_surf_elem_map[hash] = surf_elem.local_index();
    }

    for (auto& surf_elem : other)
    {
        auto hash      = utility::hash_sorted(surf_elem.connectivity());
        auto find_iter = hash_to_this_surf_elem_map.find(hash);
        if (find_iter == hash_to_this_surf_elem_map.end())
        {
            this->add(surf_elem.type(), surf_elem.connectivity());
        }
    }

    return *this;
}
}
