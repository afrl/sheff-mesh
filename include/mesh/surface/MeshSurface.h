#pragma once

#include "mesh/Mesh.h"
#include "mesh/face/Face.h"

namespace mesh
{
template<typename SurfaceElementAccessorType>
class SurfaceElementIteratorBase;
class SurfaceElementAccessor;

/// @brief Struct holding the result of MeshSurface::mesh_with_maps
struct SurfaceMeshWithMaps
{
    // A surface mesh containing only nodes that lie on a surface mesh.
    Mesh surface_mesh;
    // A map going from the volume node id to the surface id.  For nodes that appear in the volume
    // mesh but not in the surface mesh, the value will be std::numeric_limits<LocalIndex>::max()
    std::vector<LocalIndex> volume_to_surface_node_id_map;
    // A map going from the surface node id to the volume id.
    std::vector<LocalIndex> surface_to_volume_node_id_map;
};

/// @brief Stores a description of a surface of another volumetric mesh.
///
/// Internally, holds a reference to the mesh for which this surface refers.
/// Note: the connectivity of the surface elements refers to the nodes in the referenced mesh, so if
/// the nodes are rearranged in the referenced mesh, this surface will be invalid.
class MeshSurface
{
  public:
    using SurfaceElementIterator = SurfaceElementIteratorBase<SurfaceElementAccessor>;

    MeshSurface(const Mesh& volume_mesh);
    MeshSurface(MeshSurface&& other) noexcept;

    /// @name Accessors
    ///@{

    /// Enables ranged-based loop over elements.
    SurfaceElementIterator begin() const;
    /// Enables ranged-based loop over elements.
    SurfaceElementIterator end() const;

    /// Access a particular element by local index
    SurfaceElementAccessor operator[](LocalIndex i) const;

    /// Access the last element
    SurfaceElementAccessor back() const;

    /// Returns the size of the container.
    LocalIndex size() const;

    /// Creates a mesh object of the surface that no longer references a volumetric mesh and only
    /// contains nodes that appear in the volume mesh.
    Mesh mesh() const;

    /// Creates a mesh object of the surface that no longer references a volumetric mesh and only
    /// contains nodes that appear in the volume mesh, but also returns maps to go between nodes
    /// indexing in the volume and surface meshes
    SurfaceMeshWithMaps mesh_with_map() const;

    ///@}
    /// @param face

    /// @brief Adds the given face to this surface (reversing orientation to point outward instead
    /// of inward)
    void add(const FaceConstAccessor& face);
    /// @brief Adds a new surface element given the type and connectivity directly.
    /// @param conn - connectivity relative to the referenced volume mesh (assumed to be ordered to
    /// point outward).
    void add(const ElementType& type, const utility::Slice<const LocalIndex>& conn);
    void add(const ElementType& type, const utility::Slice<const LocalIndex, ConnIndex>& conn);

    // Yeilds a temporary object that holds the difference (boolean) of this and the other surface.
    MeshSurface operator-(const MeshSurface& other) const;
    MeshSurface& operator+=(const MeshSurface& other);

    friend class SurfaceElementAccessor;

  private:
    const Mesh& volume_mesh_;
    /// Stores the geometric type of each surface element.
    std::vector<ElementType> types_;
    /// Stores the connectivities of all surface elements on the partition as a 1D array.
    std::vector<LocalIndex> connectivities_;
    /// Offsets for connectivities_ for all of the surface elements
    std::vector<std::size_t> connectivity_offsets_;
};
}