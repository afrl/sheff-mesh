#include "mesh/Mesh.h"

#include "mesh/common/Permutation.h"
#include "utility/MPIUtils.h"

#include <unordered_map>
#include <unordered_set>

namespace mesh
{
Mesh::Mesh()
  : n_partitions_(1)
  , this_partition_id_(std::numeric_limits<PartitionIndex>::max())
  , nodes_()
  , elements_(nodes_)
{
}

Mesh::Mesh(const Mesh& other)
  : Mesh()
{
    copy(other);
    n_partitions_      = other.n_partitions_;
    this_partition_id_ = other.this_partition_id_;
}

Mesh::Mesh(Mesh&& other) noexcept
  : Mesh()
{
    std::swap(n_partitions_, other.n_partitions_);
    std::swap(this_partition_id_, other.this_partition_id_);
    nodes_.swap(other.nodes_);
    elements_.swap(other.elements_);
}

bool Mesh::operator==(const Mesh& other) const
{
    return n_partitions_ == other.n_partitions_ && this_partition_id_ == other.this_partition_id_ &&
           nodes_ == other.nodes_ && elements_ == other.elements_;
}

bool Mesh::contains_all_partitions() const
{
    return this_partition_id_ == std::numeric_limits<PartitionIndex>::max();
}

const PartitionIndex& Mesh::partition_id() const
{
    return this_partition_id_;
}

ComponentIndex Mesh::n_components() const
{
    return static_cast<ComponentIndex>(component_indices().size());
}

std::vector<ComponentIndex> Mesh::component_indices() const
{
    std::unordered_set<ComponentIndex> component_indices_set;
    for (auto e : elements())
    {
        component_indices_set.insert(e.component_id());
    }
    return std::vector<ComponentIndex>(component_indices_set.begin(), component_indices_set.end());
}

const PartitionIndex& Mesh::n_partitions() const
{
    return n_partitions_;
}

std::unordered_map<LocalIndex, LocalIndex> Mesh::copy(
  const Mesh& other,
  const std::vector<LocalIndex> element_indices_to_copy)
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

    LocalIndex node_offset = nodes_.size();
    std::unordered_map<LocalIndex, LocalIndex> other_to_this_node_index_map;
    for (auto i : node_indices_to_copy)
    {
        auto other_node      = other.nodes()[i];
        simple_math::Point p = { other_node.x(), other_node.y(), other_node.z() };
        nodes_.add(p, other_node.global_index(), other_node.parent_partition());
        other_to_this_node_index_map[other_node.local_index()] = nodes_.size() - 1;
    }

    for (auto& other_group : other.nodes().groups())
    {
        auto& this_group = nodes().group(other_group.name);
        for (auto& i : other_group.node_indices)
        {
            if (other_to_this_node_index_map.count(i) > 0)
            {
                this_group.node_indices.push_back(i + node_offset);
            }
        }
    }

    LocalIndex element_offset = elements_.size();
    for (auto i : element_indices_to_copy)
    {
        auto new_element = elements_.add_copy(other.elements()[i]);

        // Remap connectivity to the indices used in this mesh
        // Note: be careful with polyhedron elements
        if (new_element.type() == ElementType::POLYHEDRON)
        {
            auto conn                       = new_element.connectivity();
            mesh::ConnIndex curr_conn_index = 0;
            const auto& n_facets            = conn[curr_conn_index++];
            for (auto facet_index : utility::range(n_facets))
            {
                const auto& n_points_in_facet = conn[curr_conn_index++];
                for (auto point_index_in_facet : utility::range(n_points_in_facet))
                {
                    conn[curr_conn_index] = other_to_this_node_index_map[conn[curr_conn_index]];
                    curr_conn_index++;
                }
            }
        }
        else
        {
            for (auto& c : new_element.connectivity())
            {
                c = other_to_this_node_index_map[c];
            }
        }
    }

    if (other.elements().groups().size() > 0)
    {
        std::set<mesh::LocalIndex> element_indices_to_copy_set(element_indices_to_copy.begin(),
                                                               element_indices_to_copy.end());
        for (auto& other_group : other.elements().groups())
        {
            auto& this_group = elements().group(other_group.name);
            for (auto& i : other_group.element_indices)
            {
                if (element_indices_to_copy_set.count(i) > 0)
                {
                    this_group.element_indices.push_back(i + element_offset);
                }
            }
        }
    }

    return other_to_this_node_index_map;
}

std::unordered_map<LocalIndex, LocalIndex> Mesh::copy(
  const Mesh& other,
  const std::set<LocalIndex> element_indices_to_copy)
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

    LocalIndex node_offset = nodes_.size();
    std::unordered_map<LocalIndex, LocalIndex> other_to_this_node_index_map;
    for (auto i : node_indices_to_copy)
    {
        auto other_node      = other.nodes()[i];
        simple_math::Point p = { other_node.x(), other_node.y(), other_node.z() };
        nodes_.add(p, other_node.global_index(), other_node.parent_partition());
        other_to_this_node_index_map[other_node.local_index()] = nodes_.size() - 1;
    }

    for (auto& other_group : other.nodes().groups())
    {
        auto& this_group = nodes().group(other_group.name);
        for (auto& i : other_group.node_indices)
        {
            if (other_to_this_node_index_map.count(i) > 0)
            {
                this_group.node_indices.push_back(i + node_offset);
            }
        }
    }

    LocalIndex element_offset = elements_.size();
    for (auto i : element_indices_to_copy)
    {
        auto new_element = elements_.add_copy(other.elements()[i]);

        // Remap connectivity to the indices used in this mesh
        // Note: be careful with polyhedron elements
        if (new_element.type() == ElementType::POLYHEDRON)
        {
            auto conn                       = new_element.connectivity();
            mesh::ConnIndex curr_conn_index = 0;
            const auto& n_facets            = conn[curr_conn_index++];
            for (auto facet_index : utility::range(n_facets))
            {
                const auto& n_points_in_facet = conn[curr_conn_index++];
                for (auto point_index_in_facet : utility::range(n_points_in_facet))
                {
                    conn[curr_conn_index] = other_to_this_node_index_map[conn[curr_conn_index]];
                    curr_conn_index++;
                }
            }
        }
        else
        {
            for (auto& c : new_element.connectivity())
            {
                c = other_to_this_node_index_map[c];
            }
        }
    }

    if (other.elements().groups().size() > 0)
    {
        std::set<mesh::LocalIndex> element_indices_to_copy_set(element_indices_to_copy.begin(),
                                                               element_indices_to_copy.end());
        for (auto& other_group : other.elements().groups())
        {
            auto& this_group = elements().group(other_group.name);
            for (auto& i : other_group.element_indices)
            {
                if (element_indices_to_copy_set.count(i) > 0)
                {
                    this_group.element_indices.push_back(i + element_offset);
                }
            }
        }
    }

    return other_to_this_node_index_map;
}

std::unordered_map<LocalIndex, LocalIndex> Mesh::copy(const Mesh& other,
                                                      const LocalIndex& element_index_to_copy)
{
    return copy(other, std::vector<LocalIndex>({ element_index_to_copy }));
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

    for (auto& other_group : other.nodes().groups())
    {
        auto& this_group = nodes().group(other_group.name);
        for (auto& i : other_group.node_indices)
        {
            this_group.node_indices.push_back(i + node_offset);
        }
    }

    LocalIndex element_offset = elements_.size();
    for (auto& other_element : other.elements())
    {
        auto new_element = elements_.add_copy(other_element);

        // Remap connectivity to the indices used in this mesh
        // Note: be careful with polyhedron elements
        if (new_element.type() == ElementType::POLYHEDRON)
        {
            auto conn                       = new_element.connectivity();
            mesh::ConnIndex curr_conn_index = 0;
            const auto& n_facets            = conn[curr_conn_index++];
            for (auto facet_index : utility::range(n_facets))
            {
                const auto& n_points_in_facet = conn[curr_conn_index++];
                for (auto point_index_in_facet : utility::range(n_points_in_facet))
                {
                    conn[curr_conn_index] = conn[curr_conn_index] + node_offset;
                    curr_conn_index++;
                }
            }
        }
        else
        {
            for (auto& c : new_element.connectivity())
            {
                c = c + node_offset;
            }
        }
    }

    for (auto& other_group : other.elements().groups())
    {
        auto& this_group = elements().group(other_group.name);
        for (auto& i : other_group.element_indices)
        {
            this_group.element_indices.push_back(i + element_offset);
        }
    }
}

void Mesh::move(Mesh&& other)
{
    if (nodes().size() == 0 && elements().size() == 0)
    {
        *this = std::move(other);
    }
    else
    {
        copy(other);
    }
}

void Mesh::sort()
{
    // Local variables to bind in lambdas
    auto& nodes    = nodes_;
    auto& elements = elements_;

    std::vector<LocalIndex> backward_node_permutation(nodes_.size());
    std::iota(backward_node_permutation.begin(), backward_node_permutation.end(), 0);

    std::sort(
      backward_node_permutation.begin(),
      backward_node_permutation.end(),
      [nodes](LocalIndex lhs_index, LocalIndex rhs_index) -> bool
      {
          return nodes[lhs_index].parent_partition() < nodes[rhs_index].parent_partition() ||
                 (nodes[lhs_index].parent_partition() == nodes[rhs_index].parent_partition() &&
                  nodes[lhs_index].global_index() < nodes[rhs_index].global_index());
      });
    nodes.permute(backward_node_permutation);
    elements.permute_connectivity(detail::get_inverse_map(backward_node_permutation));

    PartitionIndex this_partition = this_partition_id_;
    if (contains_all_partitions())
    {
        this_partition = 0;
    }
    std::vector<LocalIndex> backward_element_permutation(elements_.size());
    std::iota(backward_element_permutation.begin(), backward_element_permutation.end(), 0);
    std::sort(backward_element_permutation.begin(),
              backward_element_permutation.end(),
              [elements, &this_partition](LocalIndex lhs_index, LocalIndex rhs_index) -> bool
              {
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

void Mesh::clear()
{
    n_partitions_      = 1;
    this_partition_id_ = 0;
    nodes_.clear();
    elements_.clear();
}

Mesh& Mesh::operator=(const Mesh& other)
{
    clear();
    copy(other);
    n_partitions_      = other.n_partitions_;
    this_partition_id_ = other.this_partition_id_;
    return *this;
}

Mesh& Mesh::operator=(Mesh&& other) noexcept
{
    std::swap(n_partitions_, other.n_partitions_);
    std::swap(this_partition_id_, other.this_partition_id_);
    nodes_.swap(other.nodes_);
    elements_.swap(other.elements_);
    return *this;
}

PartitionIndex& Mesh::n_partitions()
{
    return n_partitions_;
}

PartitionIndex& Mesh::partition_id()
{
    return this_partition_id_;
}


}
