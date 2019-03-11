#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"

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
}