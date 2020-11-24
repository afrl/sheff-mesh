#include "mesh/relations/NodeToElementMap.h"

namespace mesh
{
NodeToElementMap::NodeToElementMap()
  : mesh_ptr_(nullptr)
{
}

void NodeToElementMap::update(const Mesh& mesh)
{
    mesh_ptr_ = &mesh;

    parent_element_indices_.resize(mesh.nodes().size());

    for (auto& v : parent_element_indices_)
    {
        v.clear();
    }

    for (auto& element : mesh.elements())
    {
        for (auto& node : element.nodes())
        {
            parent_element_indices_[node.local_index()].push_back(element.local_index());
        }
    }
}

const std::vector<LocalIndex>& mesh::NodeToElementMap::get_element_indices(
  const LocalIndex node_index) const
{
    assert(mesh_ptr_ != nullptr);
    return parent_element_indices_[node_index];
}

IteratorRange<ElementSubrangeConstIterator> mesh::NodeToElementMap::get_elements(
  const LocalIndex node_index) const
{
    assert(mesh_ptr_ != nullptr);
    auto& parent_elements = parent_element_indices_[node_index];
    return {
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     parent_elements), // begin
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     parent_elements,
                                     static_cast<LocalIndex>(parent_elements.size())) // end
    };
}
void NodeToElementMap::adjust_for_node_combination(const LocalIndex node_to_keep,
                                                   const LocalIndex node_to_replace)
{
    std::set<LocalIndex> combined_parents;
    for (auto& p : get_element_indices(node_to_keep))
    {
        combined_parents.insert(p);
    }
    for (auto& p : get_element_indices(node_to_replace))
    {
        combined_parents.insert(p);
    }
    parent_element_indices_[node_to_replace].clear();
    parent_element_indices_[node_to_keep].clear();
    parent_element_indices_[node_to_keep].reserve(combined_parents.size());
    for (auto& p : combined_parents)
    {
        parent_element_indices_[node_to_keep].push_back(p);
    }
}
}