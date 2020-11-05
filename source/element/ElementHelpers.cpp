#include "mesh/element/ElementHelpers.h"

#include "simple_math/Utility.h"
#include "utility/Hash.h"

#include <set>

using namespace simple_math::point_operators;

namespace mesh
{
std::pair<bool, ConnIndex> find_face_index(const mesh::ElementConstAccessor element,
                                           const std::size_t conn_hash)
{
    auto face_infos = get_face_info(element.type());
    for (auto face_index : utility::range(face_infos.size()))
    {
        auto& face_info = face_infos[face_index];
        if (conn_hash == utility::hash_sorted(face_info.connectivity))
        {
            return { true, static_cast<ConnIndex>(face_index) };
        }
    }

    return { false, 0 };
}

std::pair<bool, ConnIndex> find_edge_index(const mesh::ElementConstAccessor element,
                                           const std::size_t conn_hash)
{
    auto edge_infos = get_edge_info(element.type());
    for (auto edge_index : utility::range(edge_infos.size()))
    {
        auto& edge_info = edge_infos[edge_index];
        if (conn_hash == utility::hash_sorted(edge_info.connectivity))
        {
            return { true, static_cast<ConnIndex>(edge_index) };
        }
    }

    return { false, 0 };
}
ElementTopologyComponent find_topology_component(const ElementConstAccessor element,
                                                 const std::vector<ConnIndex>& connectivity)
{
    if (connectivity.size() == 1)
    {
        return { TopologyType::NODE, connectivity.front() };
    }

    auto conn_hash = utility::hash_sorted(connectivity);
    auto result    = find_face_index(element, conn_hash);
    if (result.first == true)
    {
        return { TopologyType::FACE, result.second };
    }

    result = find_edge_index(element, conn_hash);
    if (result.first == true)
    {
        return { TopologyType::EDGE, result.second };
    }

    throw _UTIL_MSG_EXCEPT(
      "There was no sub-topology in the element found that contains the given nodes.");
}

ElementTopologyComponent find_shared_topology_component(
  const ElementConstAccessor element,
  const bool in_first_element,
  const std::vector<std::pair<ConnIndex, ConnIndex>>& shared_connectivities)
{
    std::vector<mesh::ConnIndex> relevant_shared_connectivity_indices;
    if (in_first_element)
    {
        for (auto shared_position : shared_connectivities)
        {
            relevant_shared_connectivity_indices.push_back(shared_position.first);
        }
    }
    else
    {
        for (auto shared_position : shared_connectivities)
        {
            relevant_shared_connectivity_indices.push_back(shared_position.second);
        }
    }

    return find_topology_component(element, relevant_shared_connectivity_indices);
}
std::vector<LocalIndex> find_shared_nodes(const ElementConstAccessor element_1,
                                          const ElementConstAccessor element_2)
{
    std::set<LocalIndex> common_nodes(element_1.connectivity().begin(),
                                      element_1.connectivity().end());
    for (auto other_local_index : element_2.connectivity())
    {
        auto find_iter = common_nodes.find(other_local_index);
        if (find_iter == common_nodes.end())
        {
            common_nodes.erase(find_iter);
        }
    }

    return std::vector<LocalIndex>(common_nodes.begin(), common_nodes.end());
}

std::vector<std::pair<ConnIndex, ConnIndex>> find_shared_node_positions(
  const ElementConstAccessor element_1,
  const ElementConstAccessor element_2)
{
    std::vector<std::pair<ConnIndex, ConnIndex>> shared_node_pos;
    const auto other_conn = element_2.connectivity();
    for (ConnIndex other_pos = 0; other_pos < other_conn.size(); other_pos++)
    {
        for (ConnIndex this_pos = 0; this_pos < element_1.connectivity().size(); this_pos++)
        {
            if (other_conn[other_pos] == element_1.connectivity()[this_pos])
            {
                shared_node_pos.push_back({ this_pos, other_pos });
            }
        }
    }
    return shared_node_pos;
}

Point calc_mean_coordinate(const ElementConstAccessor element)
{
    Point mean                = { 0, 0, 0 };
    const FloatType num_nodes = static_cast<FloatType>(element.connectivity().size());
    for (auto& node : element.nodes())
    {
        mean = mean + node.point();
    }
    mean = (1.0 / num_nodes) * mean;
    return mean;
}
}