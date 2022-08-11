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

    // Assumptions about the mesh made by this adjacency list algorithm:
    // - No elements are super imposed on each other
    // - Each element contains the same face no more than once
    // - No more than 2 elements can share a face (related to the first)

    // Within this functino, this type defines some adjacency relationship between two
    // elements in the form:
    // < element_1 local index, element #2 local index, int indicating type of adjacency >
    // Note: for the type int, 0 corresponds to face adjacent, 1 corresponds to edge adjacent,
    //       and 2 corresponds to node adjacent.
    using RelationTuple = std::tuple<mesh::LocalIndex, mesh::LocalIndex, int8_t>;
    // This will temporarily hold all element adjacency relations discovered
    std::vector<RelationTuple> relations;

    // Note: consider also including the local face/edge index so that if you sort by it, you know
    // what face a given element is adjacent through very efficiently
    // A list of local adjacency indices would need to be kept as well (indicating the
    // face/edge/node in the element through which an adjacent element connects to a given element)

    // Tracks the hash of all element pairs that have been added to 'relations'
    std::unordered_set<std::size_t> adjacent_element_pairs_added;

    // Loop over the elements and find all elements sharing the same face
    {
        // A map from face hash to element ids
        std::unordered_map<std::size_t, std::array<mesh::LocalIndex, 2>> face_hash_to_elements_map;
        for (auto e : mesh.elements())
        {
            FaceContainer faces(e);
            for (auto face : faces)
            {
                auto face_hash = utility::hash_sorted(face.connectivity());
                auto find_iter = face_hash_to_elements_map.find(face_hash);
                // if: a hash for this face does not exist
                if (find_iter == face_hash_to_elements_map.end())
                {
                    face_hash_to_elements_map[face_hash] = {
                        e.local_index(),
                        std::numeric_limits<mesh::LocalIndex>::max()
                    };
                }
                // else: an element has already been visited with this face
                else
                {
                    _UTIL_ASSERT_C(find_iter->second[1] ==
                                     std::numeric_limits<mesh::LocalIndex>::max(),
                                   "The algorithm used to build the element adjacency list assumes "
                                   "that no more than 2 elements can share the same face.  This "
                                   "assumption is violated for the given mesh.  There may be "
                                   "something wrong, or you might have superimposed elements.");
                    find_iter->second[1] = e.local_index();
                }
            }
        }

        // Add adjacency relations for every pair of elements sharing a face
        for (auto& map_pair : face_hash_to_elements_map)
        {
            // Only consider faces with 2 connected elements
            if (map_pair.second[0] != std::numeric_limits<mesh::LocalIndex>::max() &&
                map_pair.second[1] != std::numeric_limits<mesh::LocalIndex>::max())
            {
                // Add relation for the first element -> second element
                adjacent_element_pairs_added.insert(
                  utility::hash(map_pair.second[0], map_pair.second[1]));
                relations.push_back({ map_pair.second[0], map_pair.second[1], 0 });

                // Add relation for the second element -> first element
                adjacent_element_pairs_added.insert(
                  utility::hash(map_pair.second[1], map_pair.second[0]));
                relations.push_back({ map_pair.second[1], map_pair.second[0], 0 });
            }
        }
    }

    // Loop over the elements and find all elements sharing the same edge
    {
        // A map from edge hash to element ids
        std::unordered_map<std::size_t, std::vector<mesh::LocalIndex>> edge_hash_to_elements_map;
        std::vector<mesh::LocalIndex> edge_conn_local_indices;
        edge_conn_local_indices.reserve(20);
        for (auto e : mesh.elements())
        {
            auto edges = get_edge_info(e.type());
            for (auto edge_index : utility::range(edges.size()))
            {
                // Lists of indices indicating the nodes in the element that compose the edge.
                // Note: indices are relative to the element connectivity list.
                auto& edge_conn_element_indices = edges[edge_index].connectivity;
                // Same list as above but indices are relative to the node's local index in the mesh
                edge_conn_local_indices.resize(edge_conn_element_indices.size());
                for (auto i : utility::range(edge_conn_element_indices.size()))
                {
                    edge_conn_local_indices[i] = e.connectivity()[edge_conn_element_indices[i]];
                }

                // Use the sorted node's local indices for the hash so that the same edge in two
                // elements will have the same hash
                auto edge_hash = utility::hash_sorted(edge_conn_local_indices);
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

        for (auto& map_pair : edge_hash_to_elements_map)
        {
            // Note: loop strategically goes over the upper triangle of the adjacency matrix
            // for the current edge so that we do not need a `if (i != j)` check.
            for (auto i : utility::range(map_pair.second.size()))
            {
                for (auto j : utility::range(i + 1, map_pair.second.size()))
                {
                    // Add relation for the element i -> element j (if there was already an
                    // adjacency relation due to them sharing a face)
                    if (adjacent_element_pairs_added.find(
                          utility::hash(map_pair.second[i], map_pair.second[j])) ==
                        adjacent_element_pairs_added.end())
                    {
                        adjacent_element_pairs_added.insert(
                          utility::hash(map_pair.second[i], map_pair.second[j]));
                        relations.push_back({ map_pair.second[i], map_pair.second[j], 1 });
                    }
                    // Add relation for the element j -> element i (if there was already an
                    // adjacency relation due to them sharing a face)
                    if (adjacent_element_pairs_added.find(
                          utility::hash(map_pair.second[j], map_pair.second[i])) ==
                        adjacent_element_pairs_added.end())
                    {
                        adjacent_element_pairs_added.insert(
                          utility::hash(map_pair.second[j], map_pair.second[i]));
                        relations.push_back({ map_pair.second[j], map_pair.second[i], 1 });
                    }
                }
            }
        }
    }

    // Loop over nodes and use the node->element map to find elements that adjacent via a node
    // that do not share an edge or face (which have already been discovered)
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
                    if (adjacent_element_pairs_added.find(utility::hash(i, j)) ==
                        adjacent_element_pairs_added.end())
                    {
                        // Note: no need to add the relation to adjacent_element_pairs_added,
                        // since this is the last type and elements will not appear multiple
                        // times for the same node in the node->element map.
                        relations.push_back({ i, j, 2 });
                    }
                }
            }
        }
    }
    adjacent_element_pairs_added.clear();

    // Sort relations by 1) the first element, 2) the second element, 3) the relation type (an int
    // being one of: 0, 1, 2)
    std::sort(
      relations.begin(),
      relations.end(),
      [](const RelationTuple& lhs, const RelationTuple& rhs) -> bool
      {
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
    // Note: This algorithm was modified to ensure duplicate relations for element pairs are not
    // added, but I will leave the code relating to `curr_adjacent_elements` that makes sure we
    // do not add duplicate relations in this loop in case there is a weird case like an element
    // contains the same node more than once in its connectivity list.
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
        if (!skip_in_search_[adjacent_element_list_[i]])
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
    for (auto i : utility::range(skip_in_search_.size()))
    {
        skip_in_search_[i] = false;
    }
}
}
