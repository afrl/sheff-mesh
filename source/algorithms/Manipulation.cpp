#include "mesh/algorithms/Manipulation.h"

#include "mesh/relations/NodeToElementMap.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

using namespace simple_math::point_operators;

namespace mesh
{
void convert_to_quadratic(Mesh& mesh)
{
    // Note: currently, all of the support element types do not have face centered
    // or body centered nodes in the quadratic version.  If those elements are added
    // then this algorithm needs extra steps for faces (that would be exactly the same
    // as the edge algorithm) and body centered (which is localized to each element) nodes.

    // Create a map of edge hash id's to new nodes.
    std::map<std::size_t, LocalIndex> hash_to_new_node_map;
    for (auto& element : mesh.elements())
    {
        auto edge_info = get_edge_info(element.type());

        for (auto& edge : edge_info)
        {
            std::vector<LocalIndex> conn;
            conn.reserve(edge.connectivity.size());
            for (auto& c : edge.connectivity)
            {
                conn.push_back(element.connectivity()[c]);
            }
            auto edge_hash = utility::hash_sorted(conn);

            if (conn.size() > 2)
            {
                continue;
            }

            auto find_iter = hash_to_new_node_map.find(edge_hash);
            if (find_iter == hash_to_new_node_map.end())
            { // If it doesn't exist, add a midside node
                auto new_coord =
                  (mesh.nodes()[conn[1]].coord() + mesh.nodes()[conn[0]].coord()) / 2.0;
                mesh.nodes().add(simple_math::to_point(new_coord));
                hash_to_new_node_map[edge_hash] = mesh.nodes().size() - 1;
            }
        }
    }

    // Update the element connectivities
    std::vector<ConnIndex> midside_node_indices;
    std::vector<std::pair<ConnIndex, ConnIndex>> linear_node_pair_indices;
    for (auto element : mesh.elements())
    {
        bool is_node_level_orient = element.local_coord_sys().size() > 1;

        // Allocate the new connectivity to null pointers
        auto quadratic_type    = get_quadratic_type(element.type());
        auto n_quadratic_nodes = get_num_nodes_for_type(quadratic_type);

        std::vector<LocalIndex> quadratic_conn(n_quadratic_nodes,
                                               std::numeric_limits<LocalIndex>::max());
        QuaternionContainer quadratic_orientations;
        quadratic_orientations.initialize(element.local_coord_sys().base_type(),
                                          element.connectivity().size());
        if (element.local_coord_sys().size() == 1)
        {
            quadratic_orientations.set(0, element.local_coord_sys()[0]);
        }

        // Get the neccessary information for the conversion
        auto linear_node_indices = get_quadratic_indices_of_linear_nodes(element.type());
        get_quadratic_indices_of_midside_nodes(element.type(),
                                               midside_node_indices,
                                               linear_node_pair_indices);

        // Copy over information for corners
        for (unsigned int i = 0; i < linear_node_indices.size(); i++)
        {
            quadratic_conn[linear_node_indices[i]] = element.connectivity()[i];
            if (is_node_level_orient)
            {
                quadratic_orientations.set(linear_node_indices[i], element.local_coord_sys()[i]);
            }
        }

        // Copy over midside nodes
        for (unsigned int i = 0; i < midside_node_indices.size(); i++)
        {
            std::vector<LocalIndex> edge_conn = {
                quadratic_conn[linear_node_pair_indices[i].first],
                quadratic_conn[linear_node_pair_indices[i].second]
            };
            auto edge_hash                          = utility::hash_sorted(edge_conn);
            quadratic_conn[midside_node_indices[i]] = hash_to_new_node_map[edge_hash];

            // Assign the midside orientations to be an average of the corner nodes it lies between
            if (is_node_level_orient)
            {
                auto avg_quaternion = simple_math::Quaternion(
                  (element.local_coord_sys()[linear_node_pair_indices[i].first].coeffs() +
                   element.local_coord_sys()[linear_node_pair_indices[i].second].coeffs()) /
                  2.0);
                quadratic_orientations[midside_node_indices[i]] = avg_quaternion;
            }
        }

        // Check to make sure there were no null values left
        for (auto& i : quadratic_conn)
        {
            if (i == std::numeric_limits<LocalIndex>::max())
            {
                throw _UTIL_MSG_EXCEPT(
                  "A null value was not set for the quadratic connectivity for element " +
                  std::to_string(element.local_index()));
            }
        }

        element.type() = quadratic_type;
        mesh.elements().replace_connectivity(element.local_index(), utility::slice(quadratic_conn));
        if (quadratic_orientations.size() > 1)
        {
            element.local_coord_sys() = std::move(quadratic_orientations);
        }
    }
}

void convert_to_linear(Mesh& mesh)
{
    std::vector<LocalIndex> new_conn;
    new_conn.reserve(20);
    for (auto element : mesh.elements())
    {
        auto linear_type = get_linearized_type(element.type());
        auto permutation = get_quadratic_indices_of_linear_nodes(linear_type);
        new_conn.resize(permutation.size());
        for (auto& i : permutation)
        {
            new_conn[i] = element.connectivity()[i];
        }

        QuaternionContainer new_orientations;
        new_orientations.initialize(element.local_coord_sys().base_type(),
                                    element.connectivity().size());
        if (element.local_coord_sys().size() == 1)
        {
            new_orientations.set(0, element.local_coord_sys()[0]);
        }
        else
        {
            for (auto i : utility::range_as<ConnIndex>(permutation.size()))
            {
                new_orientations.set(i, element.local_coord_sys()[permutation[i]]);
            }
        }

        element.type() = linear_type;
        mesh.elements().replace_connectivity(element.local_index(), utility::slice(new_conn));
        element.local_coord_sys() = std::move(new_orientations);
    }

    remove_unattached_nodes(mesh);
}

void find_and_fix_errors(Mesh& mesh)
{
    bool unchecked_elements = false;
    for (auto e : mesh.elements())
    {
        if (e.type() == ElementType::T3D4N)
        {
            auto v0 = e.node(1).coord() - e.node(0).coord();
            auto v1 = e.node(2).coord() - e.node(0).coord();
            auto v2 = e.node(3).coord() - e.node(0).coord();
            if (v0.cross(v1).dot(v2) < 0.0)
            {
                std::cout << "Fixing element " << e.global_index() << std::endl;
                std::swap(e.connectivity()[1], e.connectivity()[2]);
            }
        }
        else
        {
            unchecked_elements = true;
        }
    }
    if (unchecked_elements)
    {
        std::cout
          << "Warning: some elements were not checked for errors because a case has not been "
             "implemented for them."
          << std::endl;
    }
}

void combine_nodes(Mesh& mesh, std::vector<std::pair<LocalIndex, LocalIndex>> combinations)
{
    std::sort(combinations.begin(),
              combinations.end(),
              [](const std::pair<LocalIndex, LocalIndex>& a,
                 const std::pair<LocalIndex, LocalIndex>& b) -> bool { return a.first < b.first; });

    LocalIndex nodes_removed_count = 0;
    std::vector<LocalIndex> old_to_new_node_index_map(mesh.nodes().size());
    std::iota(old_to_new_node_index_map.begin(), old_to_new_node_index_map.end(), 0);
    auto curr_map_index = 0;
    for (auto& combination : combinations)
    {
        while (curr_map_index != combination.first)
        {
            old_to_new_node_index_map[curr_map_index++] -= nodes_removed_count;
        }
        mesh.nodes().remove_stable(combination.first - nodes_removed_count++);
        curr_map_index++;
    }
    for (auto i = curr_map_index; i < old_to_new_node_index_map.size(); i++)
    {
        old_to_new_node_index_map[i] -= nodes_removed_count;
    }
    for (auto& combination : combinations)
    {
        old_to_new_node_index_map[combination.first] =
          old_to_new_node_index_map[combination.second];
    }
    mesh.elements().permute_connectivity(old_to_new_node_index_map);
}

void remove_duplicate_nodes(Mesh& mesh, const double tol)
{
    // A lambda comparator that ensure no two nodes with approximately the same position will
    // appear in the set.
    auto cmp = [&mesh, &tol](const LocalIndex& a, const LocalIndex& b)
    {
        const auto& A = mesh.nodes()[a];
        const auto& B = mesh.nodes()[b];
        return (A.x() < B.x() - tol) || ((A.y() < B.y() - tol) && fabs(A.x() - B.x()) < tol) ||
               ((A.z() < B.z() - tol) && fabs(A.x() - B.x()) < tol && fabs(A.y() - B.y()) < tol);
    };

    std::vector<std::pair<LocalIndex, LocalIndex>> combinations; // pair: from -> to
    { // Find which nodes need to be combined into other nodes
        std::set<LocalIndex, decltype(cmp)> kd_set(cmp);
        for (auto node : mesh.nodes())
        {
            auto find_iter = kd_set.find(node.local_index());
            if (find_iter == kd_set.end())
            {
                kd_set.insert(node.local_index());
            }
            else
            {
                combinations.push_back({ node.local_index(), *find_iter });
            }
        }
    }

    combine_nodes(mesh, combinations);
}

void remove_unattached_nodes(Mesh& mesh)
{
    _UTIL_ASSERT(mesh.contains_all_partitions(),
                 "Algorithm only implemented for meshes that contain all partitions.");

    std::vector<LocalIndex> node_connection_counts(mesh.nodes().size(), 0);
    for (auto element : mesh.elements())
    {
        for (auto& node_id : element.connectivity())
        {
            node_connection_counts[node_id]++;
        }
    }

    std::vector<LocalIndex> nodes_to_remove;
    for (auto i : utility::range_as<LocalIndex>(node_connection_counts.size()))
    {
        if (node_connection_counts[i] == 0)
        {
            nodes_to_remove.push_back(i);
        }
    }

    mesh.nodes().remove_stable(nodes_to_remove);
}
}