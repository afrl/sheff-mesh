#include "mesh/algorithms/Extrusion.h"

#include "mesh/common/ElementTypes.h"
#include "mesh/relations/NodeToElementMap.h"

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Helper static functions

static Mesh extrude_2D_meshalong_path_linear(const Mesh& mesh,
                                             const std::vector<double>& z_coordinates,
                                             const std::vector<double>& x_offsets,
                                             const std::vector<double>& y_offsets)
{
    NodeToElementMap node_to_element_map;
    node_to_element_map.update(mesh);

    Mesh new_mesh;

    // These will hold pointers to all the nodes that are attached to an element
    std::vector<LocalIndex> all_attached_nodes;
    for (auto& node : mesh.nodes())
    {
        if (!node_to_element_map.get_element_indices(node.local_index()).empty())
        {
            all_attached_nodes.push_back(node.local_index());
        }
    }
    std::map<LocalIndex, LocalIndex> all_attached_node_to_local_index_map;
    for (auto i : utility::indices(all_attached_nodes))
    {
        all_attached_node_to_local_index_map[all_attached_nodes[i]] = static_cast<LocalIndex>(i);
    }

    // Lay out all the nodes in the new mesh
    for (auto i : utility::indices(z_coordinates))
    {
        // Add layer of all nodes
        for (auto node_index : all_attached_nodes)
        {
            Point p = mesh.nodes()[node_index].point();
            p[0] += x_offsets[i];
            p[1] += y_offsets[i];
            p[2] = z_coordinates[i];
            new_mesh.nodes().add(p);
        }
    }

    // Now that nodes have been placed in the mesh, make the elements
    for (auto i : utility::range(z_coordinates.size() - 1))
    {
        // Loop through each element in the old mesh
        for (auto& element : mesh.elements())
        {
            std::vector<LocalIndex> connectivity;
            // Add nodes on the negative z face of the element (assuming positive direction z
            // extrusion)
            for (auto node : element.nodes())
            {
                auto index = all_attached_node_to_local_index_map[node.local_index()] +
                             i * all_attached_nodes.size();
                connectivity.push_back(static_cast<LocalIndex>(index));
            }
            // Add nodes on the positive z face of the element (assuming positive direction z
            // extrusion)
            for (auto node : element.nodes())
            {
                auto index = all_attached_node_to_local_index_map[node.local_index()] +
                             (i + 1) * all_attached_nodes.size();
                connectivity.push_back(static_cast<LocalIndex>(index));
            }

            new_mesh.elements().add(infer_type(3, static_cast<unsigned int>(connectivity.size())),
                                    connectivity);

            // new_element->set_component_id(e_ptr->get_component_id());
        }
    }

    return new_mesh;
}

static Mesh extrude_2D_meshalong_path_quadratic(const Mesh& mesh,
                                                const std::vector<double>& z_coordinates,
                                                const std::vector<double>& x_offsets,
                                                const std::vector<double>& y_offsets)
{
    throw _UTIL_MSG_EXCEPT("TODO Implement");
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions

Mesh extrude_2D_mesh(const Mesh& mesh, const std::vector<double>& z_coordinates)
{
    return extrude_2D_meshalong_path(mesh,
                                     z_coordinates,
                                     std::vector<double>(z_coordinates.size(), 0.0),
                                     std::vector<double>(z_coordinates.size(), 0.0));
}

Mesh extrude_2D_meshalong_path(const Mesh& mesh,
                               const std::vector<double>& z_coordinates,
                               const std::vector<double>& x_offsets,
                               const std::vector<double>& y_offsets)
{
    if (mesh.nodes().size() == 0)
        throw _UTIL_MSG_EXCEPT("Cannot extrude an empty mesh.");

    if (is_type_linear(mesh.elements()[0].type()))
    {
        return extrude_2D_meshalong_path_linear(mesh, z_coordinates, x_offsets, y_offsets);
    }
    else
    {
        return extrude_2D_meshalong_path_quadratic(mesh, z_coordinates, x_offsets, y_offsets);
    }
}
}