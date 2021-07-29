#include "mesh/relations/ElementAdjacencyList.h"

#include "mesh/face/Face.h"
#include "utility/Hash.h"

#include <unordered_map>
#include <unordered_set>

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


    // Use node to element map to find remaining?

    using RelationTuple = std::tuple<mesh::LocalIndex, mesh::LocalIndex, int8_t>;
    std::vector<RelationTuple> relations;


    // Note: consider also including the local face/edge index so that if you sort by it, you know
    // what face a given element is adjacent through very efficiently
    // A list of local adjacency indices would need to be kept as well (indicating the
    // face/edge/node in the elemnt through which an adjacent element connects to a given element)
    {
        // Create map from face hash to element ids and create map from edge hash to element ids
        std::unordered_map<std::size_t, std::array<mesh::LocalIndex, 2>> face_hash_to_elements_map;
        std::unordered_map<std::size_t, std::vector<mesh::LocalIndex>> edge_hash_to_elements_map;
        for (auto e : mesh.elements())
        {
            FaceContainer faces(e);
            for (auto face : faces)
            {
                auto face_hash = utility::hash_sorted(face.connectivity());
                auto find_iter = face_hash_to_elements_map.find(face_hash);
                if (find_iter == face_hash_to_elements_map.end())
                {
                    face_hash_to_elements_map[face_hash] = {
                        e.local_index(),
                        std::numeric_limits<mesh::LocalIndex>::max()
                    };
                }
                else
                {
                    assert(find_iter->second[1] == std::numeric_limits<mesh::LocalIndex>::max());
                    find_iter->second[1] = e.local_index();
                }
            }

            auto edges = get_edge_info(e.type());
            for (auto edge_index : utility::range(edges.size()))
            {
                auto& edge_local_conn = edges[edge_index].connectivity;
                std::vector<mesh::LocalIndex> edge_conn(edge_local_conn.size());
                for (auto i : utility::range(edge_local_conn.size()))
                {
                    edge_conn[i] = e.connectivity()[edge_local_conn[i]];
                }

                auto edge_hash = utility::hash_sorted(edge_conn);
                auto find_iter = edge_hash_to_elements_map.find(edge_hash);
                if (find_iter == edge_hash_to_elements_map.end())
                {
                    edge_hash_to_elements_map[edge_hash] = { e.local_index() };
                }
                else
                {
                    find_iter->second.push_back(e.local_index());
                }
            }
        }

        for (auto& map_pair : face_hash_to_elements_map)
        {
            if (map_pair.second[0] != std::numeric_limits<mesh::LocalIndex>::max() &&
                map_pair.second[1] != std::numeric_limits<mesh::LocalIndex>::max())
            {
                relations.push_back({ map_pair.second[0], map_pair.second[1], 0 });
                relations.push_back({ map_pair.second[1], map_pair.second[0], 0 });
            }
        }

        for (auto& map_pair : edge_hash_to_elements_map)
        {
            for (auto i : utility::range(map_pair.second.size()))
            {
                for (auto j : utility::range(i + 1, map_pair.second.size()))
                {
                    relations.push_back({ map_pair.second[i], map_pair.second[j], 1 });
                    relations.push_back({ map_pair.second[j], map_pair.second[i], 1 });
                }
            }
        }
    }
    // Add relations for node adjacency
    for (auto& node : mesh.nodes())
    {
        auto& parent_indices = node_to_element_map.get_element_indices(node.local_index());

        for (auto i : parent_indices)
        {
            for (auto j : parent_indices)
            {
                if (i != j)
                {
                    relations.push_back({ i, j, 2 });
                }
            }
        }
    }

    std::sort(
      relations.begin(),
      relations.end(),
      [](const RelationTuple& lhs, const RelationTuple& rhs) -> bool {
          return (std::get<0>(lhs) < std::get<0>(rhs)) ||
                 (std::get<0>(lhs) == std::get<0>(rhs) && std::get<2>(lhs) < std::get<2>(rhs)) ||
                 (std::get<0>(lhs) == std::get<0>(rhs) && std::get<2>(lhs) == std::get<2>(rhs) &&
                  std::get<1>(lhs) < std::get<1>(rhs));
      });

    adjacent_element_list_.clear();
    element_offsets_.clear();
    element_offsets_.resize(mesh.elements().size() + 1, 0);
    n_face_adjacent_elements_.clear();
    n_face_adjacent_elements_.resize(mesh.elements().size(), 0);
    n_edge_adjacent_elements_.clear();
    n_edge_adjacent_elements_.resize(mesh.elements().size(), 0);
    skip_in_search_.clear();
    skip_in_search_.resize(mesh.elements().size(), false);

    mesh::LocalIndex curr_element_index  = 0;
    int8_t curr_n_face_adjacent_elements = 0;
    int8_t curr_n_edge_adjacent_elements = 0;
    std::unordered_set<mesh::LocalIndex> curr_adjacent_elements;
    for (auto& relation : relations)
    {
        // Ensures entries exist for any elements that have no adjacencies
        while (std::get<0>(relation) != curr_element_index)
        {
            n_face_adjacent_elements_[curr_element_index] = curr_n_face_adjacent_elements;
            n_edge_adjacent_elements_[curr_element_index] = curr_n_edge_adjacent_elements;
            element_offsets_[curr_element_index + 1]      = adjacent_element_list_.size();

            curr_element_index            = std::get<0>(relation);
            curr_n_face_adjacent_elements = 0;
            curr_n_edge_adjacent_elements = 0;
            curr_adjacent_elements.clear();
        }

        if (curr_adjacent_elements.count(std::get<1>(relation)) == 0)
        {
            if (std::get<2>(relation) == 0)
            {
                curr_n_face_adjacent_elements++;
            }
            else if (std::get<2>(relation) == 1)
            {
                curr_n_edge_adjacent_elements++;
            }

            adjacent_element_list_.push_back(std::get<1>(relation));
            curr_adjacent_elements.insert(std::get<1>(relation));
        }
    }
    n_face_adjacent_elements_[curr_element_index] = curr_n_face_adjacent_elements;
    n_edge_adjacent_elements_[curr_element_index] = curr_n_edge_adjacent_elements;
    element_offsets_[curr_element_index + 1]      = adjacent_element_list_.size();
}

std::vector<LocalIndex> ElementAdjacencyList::get_adjacent_element_indices_with_skips(
  const LocalIndex element_index) const
{
    assert(mesh_ptr_ != nullptr);
    std::vector<LocalIndex> adjacent_elements;
    for (auto i = element_offsets_[element_index]; i < element_offsets_[element_index + 1]; i++)
    {
        if (!skip_in_search_[i])
        {
            adjacent_elements.push_back(adjacent_element_list_[i]);
        }
    }
    return adjacent_elements;
}

IteratorRange<ElementSubrangeConstIterator> ElementAdjacencyList::get_adjacent_elements(
  const LocalIndex element_index,
  const ElementAdjacencyType adjacency_type) const
{
    auto adjacent_indices = get_adjacent_element_indices(element_index, adjacency_type);

    return {
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     adjacent_indices.data(),
                                     adjacent_indices.size()), // begin
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     adjacent_indices.data(),
                                     adjacent_indices.size(),
                                     adjacent_indices.size()) // end
    };
}

utility::Slice<const mesh::LocalIndex, mesh::ConnIndex>
ElementAdjacencyList::get_adjacent_element_indices(const LocalIndex element_index,
                                                   const ElementAdjacencyType adjacency_type) const
{
    assert(mesh_ptr_ != nullptr);

    switch (adjacency_type)
    {
        case ElementAdjacencyType::FaceEdgeNode:
            return utility::Slice<const mesh::LocalIndex, mesh::ConnIndex>(
              &adjacent_element_list_[element_offsets_[element_index]],
              static_cast<mesh::ConnIndex>(element_offsets_[element_index + 1] -
                                           element_offsets_[element_index]));
        case ElementAdjacencyType::FaceEdge:
            return utility::Slice<const mesh::LocalIndex, mesh::ConnIndex>(
              &adjacent_element_list_[element_offsets_[element_index]],
              static_cast<mesh::ConnIndex>(n_face_adjacent_elements_[element_index] +
                                           n_edge_adjacent_elements_[element_index]));
        case ElementAdjacencyType::Face:
            return utility::Slice<const mesh::LocalIndex, mesh::ConnIndex>(
              &adjacent_element_list_[element_offsets_[element_index]],
              static_cast<mesh::ConnIndex>(n_face_adjacent_elements_[element_index]));
        default:
            throw _UTIL_MSG_EXCEPT_C("A case is not implemented for the given adjacency type.");
    }
}

inline void ElementAdjacencyList::clear_skip_flags()
{
    for (auto& flag : skip_in_search_)
    {
        flag = false;
    }
}
}
