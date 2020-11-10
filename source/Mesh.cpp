#include "mesh/Mesh.h"

#include "mesh/common/Permutation.h"

#include <unordered_map>

namespace mesh
{
Mesh::Mesh()
  : nodes_()
  , elements_(nodes_)
{
}

Mesh::Mesh(Mesh&& other) noexcept
  : nodes_()
  , elements_(nodes_)
{
    nodes_.swap(other.nodes_);
    elements_.swap(other.elements_);
}

bool Mesh::operator==(const Mesh& other) const
{
    return nodes_ == other.nodes_ && elements_ == other.elements_;
}

void Mesh::copy(const Mesh& other, const std::set<LocalIndex> element_indices_to_copy)
{
    // Determine which nodes should be kept
    std::set<LocalIndex> node_indices_to_copy;
    for (auto& element_index : element_indices_to_copy)
    {
        for (auto& node : other.elements()[element_index].nodes())
        {
            node_indices_to_copy.insert(node.local_index());
        }
    }

    std::unordered_map<LocalIndex, LocalIndex> other_to_this_node_index_map;
    for (auto i : node_indices_to_copy)
    {
        auto other_node      = other.nodes()[i];
        simple_math::Point p = { other_node.x(), other_node.y(), other_node.z() };
        nodes_.add(p, other_node.global_index(), other_node.parent_partition());
        other_to_this_node_index_map[other_node.local_index()] = nodes_.size() - 1;
    }

    for (auto i : element_indices_to_copy)
    {
        auto new_element = elements_.add_copy(other.elements()[i]);

        // Remap connectivity to the indices used in this mesh
        for (auto& c : new_element.connectivity())
        {
            c = other_to_this_node_index_map[c];
        }
    }
}

void Mesh::copy(const Mesh& other)
{
    LocalIndex node_offset = nodes_.size();
    for (auto& other_node : other.nodes())
    {
        nodes_.add(other_node.point(),
                   node_offset + other_node.global_index(),
                   other_node.parent_partition());
    }

    for (auto& other_element : other.elements())
    {
        auto new_element = elements_.add_copy(other_element);

        // Remap connectivity to the indices used in this mesh
        for (auto& c : new_element.connectivity())
        {
            c = c + node_offset;
        }
    }
}

void Mesh::sort(const PartitionIndex this_partition)
{
    // Local variables to bind in lambdas
    auto& nodes    = nodes_;
    auto& elements = elements_;

    std::vector<LocalIndex> backward_node_permutation(nodes_.size());
    std::iota(backward_node_permutation.begin(), backward_node_permutation.end(), 0);

    std::sort(
      backward_node_permutation.begin(),
      backward_node_permutation.end(),
      [nodes](LocalIndex lhs_index, LocalIndex rhs_index) -> bool {
          return nodes[lhs_index].parent_partition() < nodes[rhs_index].parent_partition() ||
                 (nodes[lhs_index].parent_partition() == nodes[rhs_index].parent_partition() &&
                  nodes[lhs_index].global_index() < nodes[rhs_index].global_index());
      });
    nodes.permute(backward_node_permutation);
    elements.permute_connectivity(detail::get_inverse_map(backward_node_permutation));


    std::vector<LocalIndex> backward_element_permutation(elements_.size());
    std::iota(backward_element_permutation.begin(), backward_element_permutation.end(), 0);
    std::sort(backward_element_permutation.begin(),
              backward_element_permutation.end(),
              [elements, &this_partition](LocalIndex lhs_index, LocalIndex rhs_index) -> bool {
                  auto lhs_element = elements[lhs_index];
                  auto rhs_element = elements[rhs_index];
                  if (lhs_element.parent_partition() == rhs_element.parent_partition())
                  {
                      return lhs_element.global_index() < rhs_element.global_index();
                  }
                  else if (lhs_element.parent_partition() == this_partition)
                  {
                      return true;
                  }
                  else if (rhs_element.parent_partition() == this_partition)
                  {
                      return false;
                  }
                  else
                  {
                      return lhs_element.parent_partition() < rhs_element.parent_partition();
                  }
              });
    elements.permute(backward_element_permutation);
}

void Mesh::clear() {
    nodes_.clear();
    elements_.clear();
}

Mesh& Mesh::operator=(const Mesh& other)
{
    clear();
    copy(other);
    return *this;
}

}
