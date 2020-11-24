#include "mesh/relations/ElementAdjacencyList.h"

namespace mesh
{
ElementAdjacencyList::ElementAdjacencyList()
  : mesh_ptr_(nullptr)
{
}

void ElementAdjacencyList::update(const Mesh& mesh)
{
    NodeToElementMap node_to_element_map;
    node_to_element_map.update(mesh);
    update(mesh, node_to_element_map);
}

void ElementAdjacencyList::update(const Mesh& mesh, const NodeToElementMap& node_to_element_map)
{
    mesh_ptr_ = &mesh;

    adjacency_list_.resize(mesh.elements().size());

    for (auto& node : mesh.nodes())
    {
        auto& parent_indices = node_to_element_map.get_element_indices(node.local_index());

        for (auto i : parent_indices)
        {
            for (auto j : parent_indices)
            {
                if (i != j)
                {
                    adjacency_list_[i].neighbor_indices.insert(j);
                }
            }
        }
    }
}

std::vector<LocalIndex> ElementAdjacencyList::get_element_indices(
  const LocalIndex element_index) const
{
    assert(mesh_ptr_ != nullptr);
    std::vector<LocalIndex> adjacent_elements;
    for (auto& neighbor_index : adjacency_list_[element_index].neighbor_indices)
    {
        if (!adjacency_list_[neighbor_index].skip_neighbor_flag)
        {
            adjacent_elements.push_back(neighbor_index);
        }
    }
    return adjacent_elements;
}

IteratorRange<ElementSubrangeConstIterator> ElementAdjacencyList::get_elements(
  const LocalIndex element_index) const
{
    assert(mesh_ptr_ != nullptr);
    auto adjacent_elements = get_element_indices(element_index);
    return {
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     adjacent_elements), // begin
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     adjacent_elements,
                                     static_cast<LocalIndex>(adjacent_elements.size())) // end
    };
}

inline void ElementAdjacencyList::clear_skip_flags()
{
    for (auto& v : adjacency_list_)
    {
        v.skip_neighbor_flag = false;
    }
}
}
