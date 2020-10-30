#pragma once

#include "mesh/Defines.h"
#include "mesh/common/EdgeTypes.h"
#include "mesh/common/FaceTypes.h"
#include "ElementAccessor.h"

namespace mesh
{

/// @brief (Work in progress) Stores information about all faces of a element.
class ElementFaces
{
  public:
    std::vector<mesh::FaceType> types_;
    /// Connectivity of faces based on the element connectivity
    std::vector<mesh::ConnIndex> connectivities_;
    /// Size of # of faces + 1
    std::vector<mesh::ConnIndex> connectivity_offsets_;
};

/// @brief (Work in progress) Stores information about all edges of a element.
class ElementEdges
{
    std::vector<mesh::EdgeType> types_;
    /// Connectivity of edges based on the element connectivity
    std::vector<mesh::ConnIndex> connectivities_;
    /// Size of # of edges + 1
    std::vector<mesh::ConnIndex> connectivity_offsets_;
};

/// @brief (Work in progress) Stores information about all faces and edges of a element.
struct ElementSubtopology
{
    // Return range of FaceAccessors
    void get_element_faces() const;
    void get_element_edges() const;
    void get_face_edges() const;
    void get_edge_faces() const;

    ElementFaces faces_;
    ElementEdges edges_;
};

/// @brief (Work in progress) Provides a map from an edge to the faces that contain it.
struct EdgeToFaceMap
{
  public:
    EdgeToFaceMap(const ElementFaces& faces, const ElementEdges& edges);
    // Return range of FaceAccessors
    void get_faces(const mesh::ConnIndex edge_index);

  private:
    const ElementFaces& faces_;
    const ElementEdges& edges_;
    std::vector<std::pair<mesh::ConnIndex, mesh::ConnIndex>> map_;
};

/// @brief (Work in progress) Provides a map from a face to the edges that compose it.
struct FaceToEdgeMap
{
  public:
    FaceToEdgeMap(const ElementFaces& faces, const ElementEdges& edges);
    // Return range of EdgeAccessors
    void get_edges(const mesh::ConnIndex face_index);

  private:
    const ElementFaces& faces_;
    const ElementEdges& edges_;
    std::vector<std::vector<mesh::ConnIndex>> map_;
};

/// @brief (Work in progress) Provides an interface to access face information.
class FaceAccessor
{
  public:
    FaceAccessor(ElementAccessor& element, ElementFaces& faces, LocalIndex starting_index = 0)
      : element_(element)
      , faces_(faces)
      , curr_index_(starting_index)
    {
    }

    /// Returns the index of the current node within the NodeContainer that contains it.
    inline const LocalIndex& local_index() const
    {
        return curr_index_;
    }

    inline NodeAccessor node(ConnIndex i) const
    {
        assert(i >= 0);
        assert(i < faces_.connectivity_offsets_[curr_index_ + 1] -
                     faces_.connectivity_offsets_[curr_index_]);
        const auto& conn_offset           = faces_.connectivity_offsets_[curr_index_];
        const auto& index_in_element_conn = faces_.connectivities_[conn_offset + i];
        return element_.node(index_in_element_conn);
    }

    inline bool operator==(const FaceAccessor& other) const
    {
        return &faces_ == &other.faces_ && curr_index_ == other.curr_index_;
    }
    inline bool operator!=(const FaceAccessor& other) const
    {
        return !(*this == other);
    }

    friend class FaceAccessor;

  private:
    ElementAccessor& element_;
    ElementFaces& faces_;
    /// Index of the current face in the faces_ object.
    ConnIndex curr_index_;

    /// Moves accessor to the next node in the container (based on local index).
    inline void next()
    {
        curr_index_++;
    }
    /// Moves accessor to the given node (does not check if index is valid).
    inline void go_to(LocalIndex i)
    {
        curr_index_ = i;
    }
};


}