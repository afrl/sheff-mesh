#pragma once

#include "NodeToElementMap.h"
#include "utility/VectorSet.h"

#include <vector>

namespace mesh
{
enum class ElementAdjacencyType
{
    FaceEdgeNode, ///< all types adjacencies
    FaceEdge,     ///< only elements adjacent by a face or edge (not a single node)
    Face,         ///< only elements adjacent by a face (not a single node or an edge)
};

/// @brief Stores the elements that are adjacent to a given element (if two elements share at least
/// one node, they are considered adjacent).
class ElementAdjacencyList
{
  public:
    ElementAdjacencyList();

    /// @brief Update the adjacency list given a mesh.  Internally creates temporary
    /// NodeToElementMap.
    void update(const Mesh& mesh);

    /// @brief Update the adjacency list given a mesh
    void update(const Mesh& mesh, const NodeToElementMap& node_to_element_map);

    /// @brief Returns the list of indices of elements that are adjacent (in any way) to the given
    /// element, while skipping elements were flagged with skip_in_neighbor_search().
    std::vector<LocalIndex> get_adjacent_element_indices_with_skips(
      const LocalIndex element_index) const;

    /// @brief Returns a range that iterates over the elements that are adjacent to the given
    /// element.
    IteratorRange<ElementSubrangeConstIterator> get_adjacent_elements(
      const LocalIndex element_index,
      const ElementAdjacencyType adjacency_type = ElementAdjacencyType::FaceEdgeNode) const;

    /// @brief Returns the indices of the elements that are adjacent to the given element.
    utility::Slice<const mesh::LocalIndex, mesh::ConnIndex> get_adjacent_element_indices(
      const LocalIndex element_index,
      const ElementAdjacencyType adjacency_type = ElementAdjacencyType::FaceEdgeNode) const;

    /// @brief Sets internal flags to skip the given elements in subsequent calls to
    /// ElementAdjacencyList::get_element_indices and ElementAdjacencyList::get_elements
    /// @tparam IndexContainer - any container of LocalIndex obecjts that supports range-based loops
    /// @param element_indices_to_skip - elements to skip in subsequent searchers until
    ///        ElementAdjacencyList::skip_in_neighbor_search is called again or
    ///        ElementAdjacencyList::clear_skip_flags is called.
    template<typename IndexContainer>
    inline void skip_in_neighbor_search(const IndexContainer& element_indices_to_skip)
    {
        for (auto& flag : skip_in_search_)
        {
            flag = false;
        }

        for (auto& i : element_indices_to_skip)
        {
            skip_in_search_[i] = true;
        }
    }

    /// @brief Clears any flags set to skip elements in neighbor searchers from a previous call to
    /// ElementAdjacencyList::skip_in_neighbor_search.
    void clear_skip_flags();

    inline bool initialized() const
    {
        return adjacent_element_list_.size() > 0 ||
               (mesh_ptr_ && mesh_ptr_->elements().size() == 0);
    }

  private:
    /// The unraveled list of adjacent elements for every element in the mesh.  This may be sorted
    /// due to the nature of how it is created but do not count on the order.
    std::vector<LocalIndex> adjacent_element_list_;
    /// The offset in adjacent_element_list_ for a given element (the next offset indicates the end
    /// fo the element's range)
    std::vector<std::size_t> element_offsets_;
    /// Indicates how many elements are share a face with a given element (max # for a given element
    /// is the # of faces for the element and min is 0).
    /// Size: # elements in the mesh
    std::vector<mesh::ConnIndex> n_face_adjacent_elements_;
    /// Indicates how many elements are share an edge with a given element (max # for a given
    /// element is the # of faces for the element and min is 0). Size: # elements in the mesh
    std::vector<mesh::ConnIndex> n_edge_adjacent_elements_;
    /// A flag indicating whether each element should be skipped in neighbor searchers (allows
    /// functionality like skipping ghost elements)
    std::vector<bool> skip_in_search_;

    /// A non-owning pointer that is null until update is called
    /// Needed for get_elements() to return the element accessors
    const Mesh* mesh_ptr_;
};
}