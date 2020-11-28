#include "mesh/algorithms/Manipulation.h"

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
    for (auto& element : mesh.elements())
    {
        bool is_node_level_orient = element.local_coord_sys().size() > 1;

        // Allocate the new connectivity to null pointers
        auto quadratic_type    = get_quadratic_type(element.type());
        auto n_quadratic_nodes = get_num_nodes_for_type(quadratic_type);

        std::vector<LocalIndex> quadratic_conn(n_quadratic_nodes, std::numeric_limits<LocalIndex>::max());
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

        // Copy over the needed pointers
        for (unsigned int i = 0; i < linear_node_indices.size(); i++)
        {
            quadratic_conn[linear_node_indices[i]] = element.connectivity()[i];
            if (is_node_level_orient)
            {
                quadratic_orientations.set(linear_node_indices[i], element.local_coord_sys()[i]);
            }
        }

        // Copy over midside node pointers
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
}