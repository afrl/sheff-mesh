#pragma once

#include "NodeToElementMap.h"
#include "utility/VectorSet.h"

#include <vector>

namespace mesh
{
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

    /// @brief Returns the list of indices of elements that include the given node.
    std::vector<LocalIndex> get_element_indices(const LocalIndex element_index) const;

    /// @brief Returns a range that iterates over the elements that include the given node.
    IteratorRange<ElementSubrangeConstIterator> get_elements(const LocalIndex element_index) const;

    /// @brief Sets internal flags to skip the given elements in subsequent calls to
    /// ElementAdjacencyList::get_element_indices and ElementAdjacencyList::get_elements
    /// @tparam IndexContainer - any container of LocalIndex obecjts that supports range-based loops
    /// @param element_indices_to_skip - elements to skip in subsequent searchers until
    ///        ElementAdjacencyList::skip_in_neighbor_search is called again or
    ///        ElementAdjacencyList::clear_skip_flags is called.
    template<typename IndexContainer>
    inline void skip_in_neighbor_search(const IndexContainer& element_indices_to_skip)
    {
        for (auto& v : adjacency_list_)
        {
            v.skip_neighbor_flag = false;
        }

        for (auto& i : element_indices_to_skip)
        {
            adjacency_list_[i].skip_neighbor_flag = true;
        }
    }

    /// @brief Clears any flags set to skip elements in neighbor searchers from a previous call to
    /// ElementAdjacencyList::skip_in_neighbor_search.
    void clear_skip_flags();

    inline bool initialized() const
    {
        return adjacency_list_.size() > 0 || (mesh_ptr_ && mesh_ptr_->elements().size() == 0);
    }

  private:
    /// @brief A helper struct (nested in ElementAdjacencyList) for storing the adjacency list
    struct AdjacencyVertex
    {
        AdjacencyVertex()
          : skip_neighbor_flag(false)
        {
        }

        /// A flag indicating if the this vertex should be skipped when accumulating neighbors
        /// during a search
        bool skip_neighbor_flag;

        /// The local index of the neighbors that are adjacent to this vertex
        utility::VectorSet<LocalIndex> neighbor_indices;
    };

    std::vector<AdjacencyVertex> adjacency_list_;

    /// A non-owning pointer that is null until update is called
    /// Needed for get_elements() to return the element accessors
    const Mesh* mesh_ptr_;
};
}