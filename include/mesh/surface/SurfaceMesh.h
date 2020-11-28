#pragma once

#include "mesh/Mesh.h"
#include "mesh/face/Face.h"

namespace mesh
{
class SurfaceElementIterator;
class SurfaceElementAccessor;

class SurfaceMesh
{
  public:
    SurfaceMesh(Mesh& volume_mesh);
    SurfaceMesh(SurfaceMesh&& other) noexcept;

    /// @name Accessors
    ///@{

    /// Enables ranged-based loop over elements.
    SurfaceElementIterator begin();
    /// Enables ranged-based loop over elements.
    SurfaceElementIterator end();

    /// Const access a particular element by local index
    SurfaceElementAccessor operator[](LocalIndex i);

    /// Const access the last element
    SurfaceElementAccessor back();

    /// Returns the size of the container.
    LocalIndex size() const;

    /// @brief Creates a mesh object of the surface
    Mesh mesh();

    ///@}
    /// @param face

    void add(const FaceConstAccessor& face);

    friend class SurfaceElementAccessor;

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