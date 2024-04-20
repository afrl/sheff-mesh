#pragma once

#include "mesh/common/IteratorRange.h"
#include "mesh/element/ElementAccessor.h"
#include "mesh/face/FaceContainer.h"
#include "mesh/node/NodeSubrangeIterator.h"

#include <assert.h>

namespace mesh
{
/// @brief Provides an interface to access element information.
class FaceConstAccessor
{
  public:
    FaceConstAccessor(const FaceContainer& faces, ConnIndex starting_index = 0)
      : faces_(faces)
      , curr_face_index_(starting_index)
    {
    }

    /// Returns the index of the current face within the element that contains it.
    inline const ConnIndex& index() const
    {
        return curr_face_index_;
    }

    /// Returns the type of element (in terms of topology), see ElementType
    inline auto& type() const
    {
        return faces_.types_[curr_face_index_];
    }

    /// Returns the ordered local indices of nodes that form the connectivity of the current face.
    inline std::vector<LocalIndex> connectivity() const
    {
        const auto& conn_index_start = faces_.connectivity_offsets_[curr_face_index_];
        const auto& conn_size =
          faces_.connectivity_offsets_[curr_face_index_ + 1] - conn_index_start;

        std::vector<LocalIndex> conn;
        conn.reserve(conn_size);
        auto element_conn = faces_.element_.connectivity();
        for (auto element_conn_index :
             utility::Slice<const mesh::ConnIndex>(&faces_.connectivities_[conn_index_start],
                                                   conn_size))
        {
            conn.push_back(element_conn[element_conn_index]);
        }

        return conn;
    }

    /// @brief Returns the number of nodes that compose the face (size of the connectivity).
    inline ConnIndex n_nodes() const
    {
        const auto& conn_index_start = faces_.connectivity_offsets_[curr_face_index_];
        const auto& conn_size =
          faces_.connectivity_offsets_[curr_face_index_ + 1] - conn_index_start;
        return conn_size;
    }

    /// Returns a slice of the ordered connectivity of the face with indices relative to the parent
    /// element's connectivity.
    inline utility::Slice<const mesh::ConnIndex> element_conn_indices() const
    {
        const auto& conn_index_start = faces_.connectivity_offsets_[curr_face_index_];
        const auto& conn_size =
          faces_.connectivity_offsets_[curr_face_index_ + 1] - conn_index_start;

        return utility::Slice<const mesh::ConnIndex>(&faces_.connectivities_[conn_index_start],
                                                     conn_size);
    }

    /// Allows const-access to a particular node in the face.
    inline NodeConstAccessor node(ConnIndex i) const
    {
        const auto& node_conn_index =
          faces_.connectivities_[faces_.connectivity_offsets_[curr_face_index_] + i];
        return faces_.element_.node(node_conn_index);
    }

    /// Returns a range over the ordered nodes that compose the connectivity of the element.  The
    /// order will match the indices returned by connectivity().
    /// Note: iterators inherit const-ness of the element accessor
    IteratorRange<NodeSubrangeConstIterator> nodes() const
    {
        auto nodes_to_visit = connectivity();
        return {
            NodeSubrangeConstIterator(faces_.element_.node_container(),
                                      nodes_to_visit), // begin
            NodeSubrangeConstIterator(faces_.element_.node_container(),
                                      nodes_to_visit,
                                      static_cast<LocalIndex>(nodes_to_visit.size())) // end
        };
    }

    /// @Returns an accessor to the element that contains the face that was used to initialize the
    /// FaceContainer (which may be only one of the two elements that contain the face).
    ElementConstAccessor parent_element() const
    {
        return faces_.element_;
    }

    inline bool operator==(const FaceConstAccessor& other) const
    {
        return &faces_ == &other.faces_ && curr_face_index_ == other.curr_face_index_;
    }

    /// Allow to iterator accessor to next() and goto() methods
    friend class FaceConstIterator;

  private:
    const FaceContainer& faces_;
    ConnIndex curr_face_index_;

    /// Moves accessor to the next element in the container (based on local index).
    inline void next()
    {
        curr_face_index_++;
    }
    /// Moves accessor to the given element (does not check if index is valid).
    inline void go_to(ConnIndex i)
    {
        curr_face_index_ = i;
    }
};

using ElementAccessor      = ElementAccessorBase<NodeContainer, ElementContainer>;
using ElementConstAccessor = ElementAccessorBase<const NodeContainer, const ElementContainer>;
}