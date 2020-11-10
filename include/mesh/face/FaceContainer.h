#pragma once

#include "mesh/Defines.h"
#include "mesh/common/FaceTypes.h"
#include "mesh/element/ElementAccessor.h"

#include <memory>
#include <vector>

namespace mesh
{
// Forward declarations
class FaceConstAccessor;
class FaceConstIterator;
// End forward declarations

/// @brief (Work in progress) Stores information about all faces of a element.
class FaceContainer
{
  public:
    /// A reference to the node container held to allow the construction the element
    /// accessors/iterators, but it is not used within this class.
    FaceContainer(ElementConstAccessor element);

    /// @name Accessors
    ///@{
    /// Enables ranged-based loop over elements (const).
    FaceConstIterator begin() const;
    /// Enables ranged-based loop over elements.
    FaceConstIterator end() const;

    /// Const access a particular element by local index
    FaceConstAccessor operator[](ConnIndex i) const;

    /// Const access the last element
    FaceConstAccessor back() const;

    /// Returns the size of the container.
    ConnIndex size() const;
    ///@}

    friend class FaceConstAccessor;

  private:
    ElementConstAccessor element_;

    /// Stores the geometric type of each face.  @see FaceType
    /// Expected size: # of elements on partition.
    std::vector<mesh::FaceType> types_;
    /// Stores the connectivities of all faces in the element as a 1D array.  The connectivity
    /// is described as an ordered list of connectivity (relative to the element) indices of nodes.
    /// The vector connectivity_offsets_ indicates which subrange of the element's connectivity that
    /// corresponds to a given face.  For example, face i has indices in the element's connectivity
    /// of connectivities_[connectivity_offsets_[i]] ...
    /// connectivities_[connectivity_offsets_[i+1]].
    std::vector<mesh::ConnIndex> connectivities_;
    /// Expected size: # of faces in the element + 1.  @see FaceContainer::connectivities_ for
    /// more details
    std::vector<mesh::ConnIndex> connectivity_offsets_;
};
}