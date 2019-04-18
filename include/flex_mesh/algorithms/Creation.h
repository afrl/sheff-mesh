#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/maps/LocalNodeToElementMap.h"

#include "utility/Range.h"

#include <vector>

namespace flex_mesh
{
    /// Generate a 2D rectangular mesh with a uniform grid
    /// \param x_coordinates - coordinates of the corners of the grid cells along x-axis
    /// \param y_coordinates - coordinates of the corners of the grid cells along y-axis
    /// \param is_quadratic - indicates whether to use quadratic or linear elements
    ///        (defaults to quadratic)
    template<typename MeshType>
    static MeshType make_grid(
        const std::vector<double> & x_coordinates,
        const std::vector<double> & y_coordinates,
        const bool & is_quadratic = true);

    /// Generate a 3D rectangular mesh with a uniform grid
    /// \param x_coordinates - coordinates of the corners of the grid cells along x-axis
    /// \param y_coordinates - coordinates of the corners of the grid cells along y-axis
    /// \param z_coordinates - coordinates of the corners of the grid cells along z-axis
    /// \param is_quadratic - indicates whether to use quadratic or linear elements
    ///        (defaults to quadratic)
    template<typename MeshType>
    static MeshType make_grid(
        const std::vector<double> & x_coordinates,
        const std::vector<double> & y_coordinates,
        const std::vector<double> & z_coordinates,
        const bool & is_quadratic = true);

    /// Extrude a 2D mesh along z direction
    template<typename MeshType>
    static MeshType extrude_2D_mesh(
        MeshType & mesh,
        const std::vector<double> & z_coordinates);

    /// Extrude a 2D mesh along z direction with some offsets in x and y
    template<typename MeshType>
    static MeshType extrude_2D_meshalong_path(
        MeshType & mesh,
        const std::vector<double> & z_coordinates,
        const std::vector<double> & x_offsets,
        const std::vector<double> & y_offsets);

    namespace detail
    {
        template<typename MeshType>
        static MeshType extrude_2D_meshalong_path_linear(
            MeshType & mesh,
            const std::vector<double> & z_coordinates,
            const std::vector<double> & x_offsets,
            const std::vector<double> & y_offsets);

        template<typename MeshType>
        static MeshType extrude_2D_meshalong_path_quadratic(
            MeshType & mesh,
            const std::vector<double> & z_coordinates,
            const std::vector<double> & x_offsets,
            const std::vector<double> & y_offsets);
    }

    /////////////////
    // Definitions //
    /////////////////

    template<typename MeshType>
    MeshType make_grid(const std::vector<double>& x_coordinates, const std::vector<double>& y_coordinates, const bool & is_quadratic)
    {
        using NodeVector = std::vector<typename MeshType::NodeTypename *>;

        MeshType mesh;
        if (is_quadratic)
        {
            // Create nodes
            NodeVector indexed_nodes;
            for (auto i : utility::indices(x_coordinates))
            {
                // For all i's but the first one, create the row of midside nodes
                if (i != 0)
                {
                    const double x_mid = (x_coordinates[i] + x_coordinates[i - 1]) / 2.0;
                    for (auto j : utility::indices(y_coordinates))
                    {
                        auto n_mid = mesh.add_node_by_coord({ x_mid, y_coordinates[j], 0 });
                        indexed_nodes.push_back(n_mid);
                    }
                }

                // Along rows of corner nodes, we need midside and corner nodes
                const double x_corner = x_coordinates[i];
                for (auto j : utility::indices(y_coordinates))
                {
                    // Create midside node
                    if (j != 0)
                    {
                        const double y_mid = (y_coordinates[j] + y_coordinates[j - 1]) / 2.0;
                        auto n_mid = mesh.add_node_by_coord({ x_corner, y_mid, 0 });
                        indexed_nodes.push_back(n_mid);
                    }

                    // Create corner node
                    const double y_corner = y_coordinates[j];
                    auto n_corner = mesh.add_node_by_coord({ x_corner, y_corner, 0 });
                    indexed_nodes.push_back(n_corner);
                }
            }

            // Create elements
            size_t n_elements_x = x_coordinates.size() - 1;
            size_t n_elements_y = y_coordinates.size() - 1;
            for (size_t i : utility::range<size_t>(0, n_elements_x))
            {
                for (size_t j : utility::range<size_t>(0, n_elements_y))
                {
                    NodeVector elem_connectivity(8);
                    elem_connectivity[0] = indexed_nodes[i*(n_elements_y * 3 + 2) + j * 2];
                    elem_connectivity[1] = indexed_nodes[i*(n_elements_y * 3 + 2) + n_elements_y * 2 + 1 + j];
                    elem_connectivity[2] = indexed_nodes[(i + 1)*(n_elements_y * 3 + 2) + j * 2];
                    elem_connectivity[3] = indexed_nodes[(i + 1)*(n_elements_y * 3 + 2) + j * 2 + 1];
                    elem_connectivity[4] = indexed_nodes[(i + 1)*(n_elements_y * 3 + 2) + j * 2 + 2];
                    elem_connectivity[5] = indexed_nodes[i*(n_elements_y * 3 + 2) + n_elements_y * 2 + 1 + j + 1];
                    elem_connectivity[6] = indexed_nodes[i*(n_elements_y * 3 + 2) + j * 2 + 2];
                    elem_connectivity[7] = indexed_nodes[i*(n_elements_y * 3 + 2) + j * 2 + 1];
                    auto e_ptr = mesh.elements().add();
                    e_ptr->initialize(ElementType::Q2D8N, elem_connectivity);
                }
            }
        }
        // Linear grid
        else
        {
            // Create nodes
            NodeVector indexed_nodes;
            for (auto i : utility::indices(x_coordinates))
            {
                for (auto j : utility::indices(y_coordinates))
                {
                    auto n_corner = mesh.add_node_by_coord({ x_coordinates[i], y_coordinates[j], 0 });
                    indexed_nodes.push_back(n_corner);
                }
            }

            // Create elements
            size_t n_elements_x = x_coordinates.size() - 1;
            size_t n_elements_y = y_coordinates.size() - 1;
            for (size_t i : utility::range<size_t>(0, n_elements_x))
            {
                for (size_t j : utility::range<size_t>(0, n_elements_y))
                {
                    NodeVector elem_connectivity(4);
                    elem_connectivity[0] = indexed_nodes[i + j * x_coordinates.size()];
                    elem_connectivity[1] = indexed_nodes[(i + 1) + j * x_coordinates.size()];
                    elem_connectivity[2] = indexed_nodes[(i + 1) + (j + 1) * x_coordinates.size()];
                    elem_connectivity[3] = indexed_nodes[i + (j + 1) * x_coordinates.size()];
                    auto e_ptr = mesh.elements().add();
                    e_ptr->initialize(ElementType::Q2D4N, elem_connectivity);
                }
            }
        }

        mesh.update_global_numbering();
        return mesh;
    }

    template<typename MeshType>
    MeshType make_grid(const std::vector<double>& x_coordinates, const std::vector<double>& y_coordinates, const std::vector<double>& z_coordinates, const bool & is_quadratic)
    {
        auto mesh_2D = make_grid<MeshType>(x_coordinates, y_coordinates, is_quadratic);
        return extrude_2D_mesh(mesh_2D, z_coordinates);
    }

    template<typename MeshType>
    MeshType extrude_2D_mesh(MeshType & mesh, const std::vector<double> & z_coordinates)
    {
        return extrude_2D_meshalong_path(mesh,
            z_coordinates,
            std::vector<double>(z_coordinates.size(), 0.0),
            std::vector<double>(z_coordinates.size(), 0.0));
    }

    template<typename MeshType>
    MeshType extrude_2D_meshalong_path(MeshType & mesh, const std::vector<double>& z_coordinates, const std::vector<double>& x_offsets, const std::vector<double>& y_offsets)
    {
        if (mesh.nodes().size() == 0)
            throw _UTIL_MSG_EXCEPT("Cannot extrude an empty mesh.");

        if (is_type_linear(mesh.elements().front()->get_type()))
        {
            return detail::extrude_2D_meshalong_path_linear(mesh,
                z_coordinates,
                x_offsets,
                y_offsets);
        }
        else
        {
            return detail::extrude_2D_meshalong_path_quadratic(mesh,
                z_coordinates,
                x_offsets,
                y_offsets);
        }
    }

    namespace detail
    {
        template<typename MeshType>
        MeshType extrude_2D_meshalong_path_linear(MeshType & mesh, const std::vector<double>& z_coordinates, const std::vector<double>& x_offsets, const std::vector<double>& y_offsets)
        {
            mesh.update_global_numbering();
            LocalNodeToElementMap<MeshType> node_to_element_map;
            node_to_element_map.update(mesh);

            using NodePtr = typename MeshType::NodeTypename *;
            MeshType new_mesh;

            // These will hold pointers to all the nodes that are attached to an element
            std::vector<NodePtr> all_attached_nodes;
            for (auto& n_ptr : mesh.nodes())
            {
                if (!node_to_element_map.get_element_indices(n_ptr->get_local_number()).empty())
                {
                    all_attached_nodes.push_back(n_ptr.get());
                }
            }
            std::map<NodePtr, GlobalIndex> all_attached_node_to_local_index_map;
            for (auto i : utility::indices(all_attached_nodes))
            {
                all_attached_node_to_local_index_map[all_attached_nodes[i]] = static_cast<GlobalIndex>(i);
            }

            // Lay out all the nodes in the new mesh
            for (auto i : utility::indices(z_coordinates))
            {
                //Add layer of all nodes
                for (auto n_ptr : all_attached_nodes)
                {
                    Point p = n_ptr->get_coord();
                    p[0] += x_offsets[i];
                    p[1] += y_offsets[i];
                    p[2] = z_coordinates[i];
                    new_mesh.add_node_by_coord(p);
                }
            }

            new_mesh.update_global_numbering();

            // Now that nodes have been placed in the mesh, make the elements
            for (auto i : utility::range_from_zero(z_coordinates.size() - 1))
            {
                // Loop through each element in the old mesh
                for (auto& e_ptr : mesh.elements())
                {
                    std::vector<NodePtr> connectivity;
                    // Add nodes on the negative z face of the element (assuming positive direction z extrusion)
                    for (auto n_ptr : *e_ptr)
                    {
                        auto index = all_attached_node_to_local_index_map[n_ptr] + i * all_attached_nodes.size();
                        connectivity.push_back(new_mesh.nodes()[index]);
                    }
                    // Add nodes on the positive z face of the element (assuming positive direction z extrusion)
                    for (auto n_ptr : *e_ptr)
                    {
                        auto index = all_attached_node_to_local_index_map[n_ptr] + (i + 1) * all_attached_nodes.size();
                        connectivity.push_back(new_mesh.nodes()[index]);
                    }

                    auto new_element = new_mesh.elements().add();
                    new_element->initialize(infer_type(3, static_cast<unsigned int>(connectivity.size())), connectivity);
                    new_element->set_component_id(e_ptr->get_component_id());
                }
            }

            new_mesh.update_global_numbering();
            return new_mesh;
        }

        template<typename MeshType>
        MeshType extrude_2D_meshalong_path_quadratic(MeshType & mesh, const std::vector<double>& z_coordinates, const std::vector<double>& x_offsets, const std::vector<double>& y_offsets)
        {
            throw _UTIL_MSG_EXCEPT("TODO Implement");
        }
    }
}