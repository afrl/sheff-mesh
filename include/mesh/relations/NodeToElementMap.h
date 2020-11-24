#pragma once

#include "mesh/Mesh.h"
#include "mesh/element/ElementSubrangeIterator.h"

#include <vector>

namespace mesh
{
/// @brief Stores a map from the local index number of a node to the local indices of all of the
/// elements that contain the respective node.
class NodeToElementMap
{
  public:
    NodeToElementMap();

    /// @brief Update the map given a mesh
    void update(const Mesh& mesh);

    /// @brief Returns the list of indices of elements that include the given node.
    const std::vector<LocalIndex>& get_element_indices(const LocalIndex node_index) const;

    /// @brief Returns a range that iterates over the elements that include the given node.
    IteratorRange<ElementSubrangeConstIterator> get_elements(const LocalIndex node_index) const;

    /// @brief Performs a cheap update to account for combining two nodes, avoiding a complete
    /// rebuild of the maps.  It moves the parents from the node_to_replace over to the node_to_keep
    /// leaving the parent list for the node_to_replace empty.
    void adjust_for_node_combination(const LocalIndex node_to_keep,
                                     const LocalIndex node_to_replace);

  private:
    /// Outer vector has the size of the nodes in the mesh, and the inner vector stores the parent
    /// element indices for the respective node.
    std::vector<std::vector<LocalIndex>> parent_element_indices_;

    /// A non-owning pointer that is null until update is called
    const Mesh* mesh_ptr_;
};
}