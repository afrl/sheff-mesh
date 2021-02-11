#pragma once

#include "mesh/Mesh.h"
#include "mesh/face/Face.h"

namespace mesh
{
template<typename SurfaceElementAccessorType>
class SurfaceElementIteratorBase;
class SurfaceElementAccessor;
class SurfaceElementConstAccessor;

/// @brief Stores a description of a surface of another volumetric mesh.
///
/// Internally, holds a reference to the mesh for which this surface refers.
/// Note: the connectivity of the surface elements refers to the nodes in the referenced mesh, so if
/// the nodes are rearranged in the referenced mesh, this surface will be invalid.
class MeshSurface
{
  public:
    using SurfaceElementIterator      = SurfaceElementIteratorBase<SurfaceElementAccessor>;
    using SurfaceElementConstIterator = SurfaceElementIteratorBase<SurfaceElementConstAccessor>;

    MeshSurface(Mesh& volume_mesh);
    MeshSurface(MeshSurface&& other) noexcept;

    /// @name Accessors
    ///@{

    /// Enables ranged-based loop over elements.
    SurfaceElementIterator begin();
    /// Enables const ranged-based loop over elements.
    SurfaceElementConstIterator begin() const;
    /// Enables ranged-based loop over elements.
    SurfaceElementIterator end();
    /// Enables const ranged-based loop over elements.
    SurfaceElementConstIterator end() const;

    /// Access a particular element by local index
    SurfaceElementAccessor operator[](LocalIndex i);
    /// Const access a particular element by local index
    SurfaceElementConstAccessor operator[](LocalIndex i) const;

    /// Access the last element
    SurfaceElementAccessor back();
    /// Const access the last element
    SurfaceElementConstAccessor back() const;

    /// Returns the size of the container.
    LocalIndex size() const;

    /// Creates a mesh object of the surface that no longer references a volumetric mesh.
    Mesh mesh();

    ///@}
    /// @param face

    /// @brief Adds the given face to this surface (reversing orientation to point outward instead
    /// of inward)
    void add(const FaceConstAccessor& face);
    /// @brief Adds a new surface element given the type and connectivity directly.
    /// @param conn - connectivity relative to the referenced volume mesh (assumed to be ordered to
    /// point outward).
    void add(const ElementType& type, const utility::Slice<const LocalIndex>& conn);

    // Yeilds a temporary object that holds the difference (boolean) of this and the other surface.
    MeshSurface operator-(const MeshSurface& other) const;
    MeshSurface& operator+=(const MeshSurface& other);

    friend class SurfaceElementAccessor;
    friend class SurfaceElementConstAccessor;

  private:
    Mesh& volume_mesh_;
    /// Stores the geometric type of each surface element.
    std::vector<ElementType> types_;
    /// Stores the connectivities of all surface elements on the partition as a 1D array.
    std::vector<LocalIndex> connectivities_;
    /// Offsets for connectivities_ for all of the surface elements
    std::vector<std::size_t> connectivity_offsets_;
};
}