#include "mesh/algorithms/Creation.h"

#include "mesh/algorithms/Extrusion.h"
#include "simple_math/Utility.h"

namespace mesh
{
Mesh make_linear_grid(const std::vector<double>& x_coordinates,
                      const std::vector<double>& y_coordinates)
{
    Mesh mesh;

    // Create nodes
    for (auto i : utility::indices(y_coordinates))
    {
        for (auto j : utility::indices(x_coordinates))
        {
            mesh.nodes().add({ x_coordinates[j], y_coordinates[i], 0 });
        }
    }

    // Create elements
    auto n_elements_x = static_cast<LocalIndex>(x_coordinates.size() - 1);
    auto n_elements_y = static_cast<LocalIndex>(y_coordinates.size() - 1);
    auto n_nodes_x    = static_cast<LocalIndex>(x_coordinates.size());
    for (auto i : utility::range<LocalIndex>(0, n_elements_y))
    {
        for (auto j : utility::range<LocalIndex>(0, n_elements_x))
        {
            std::vector<LocalIndex> elem_connectivity(4);
            elem_connectivity[0] = j + i * n_nodes_x;
            elem_connectivity[1] = (j + 1) + i * n_nodes_x;
            elem_connectivity[2] = (j + 1) + (i + 1) * n_nodes_x;
            elem_connectivity[3] = j + (i + 1) * n_nodes_x;
            mesh.elements().add(ElementType::Q2D4N, elem_connectivity);
        }
    }

    return mesh;
}

Mesh make_linear_grid(const std::vector<double>& x_coordinates,
                      const std::vector<double>& y_coordinates,
                      const std::vector<double>& z_coordinates)
{
    auto mesh_2D = make_linear_grid(x_coordinates, y_coordinates);
    return extrude_2D_mesh(mesh_2D, z_coordinates);
}

Mesh make_uniform_linear_grid(const Point& min,
                              const Point& max,
                              const std::array<int, 3>& n_elements)
{
    auto x_coords = simple_math::linspace(min[0], max[0], n_elements[0]);
    auto y_coords = simple_math::linspace(min[1], max[1], n_elements[1]);
    auto z_coords = simple_math::linspace(min[2], max[2], n_elements[2]);

    return make_linear_grid(x_coords, y_coords, z_coords);
}

Mesh make_quadratic_grid(const std::vector<double>& x_coordinates,
                         const std::vector<double>& y_coordinates)
{
    Mesh mesh;
    
    throw _UTIL_MSG_EXCEPT("TODO swap order of x and y to preserve +z normals.");

    // Create nodes
    for (auto i : utility::indices(x_coordinates))
    {
        // For all i's but the first one, create the row of midside nodes
        if (i != 0)
        {
            const double x_mid = (x_coordinates[i] + x_coordinates[i - 1]) / 2.0;
            for (auto j : utility::indices(y_coordinates))
            {
                mesh.nodes().add({ x_mid, y_coordinates[j], 0 });
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
                mesh.nodes().add({ x_corner, y_mid, 0 });
            }

            // Create corner node
            const double y_corner = y_coordinates[j];
            mesh.nodes().add({ x_corner, y_corner, 0 });
        }
    }

    // Create elements
    auto n_elements_x = static_cast<LocalIndex>(x_coordinates.size() - 1);
    auto n_elements_y = static_cast<LocalIndex>(y_coordinates.size() - 1);
    for (auto i : utility::range<LocalIndex>(0, n_elements_x))
    {
        for (auto j : utility::range<LocalIndex>(0, n_elements_y))
        {
            std::vector<LocalIndex> elem_connectivity(8);
            elem_connectivity[0] = i * (n_elements_y * 3 + 2) + j * 2;
            elem_connectivity[1] = i * (n_elements_y * 3 + 2) + n_elements_y * 2 + 1 + j;
            elem_connectivity[2] = (i + 1) * (n_elements_y * 3 + 2) + j * 2;
            elem_connectivity[3] = (i + 1) * (n_elements_y * 3 + 2) + j * 2 + 1;
            elem_connectivity[4] = (i + 1) * (n_elements_y * 3 + 2) + j * 2 + 2;
            elem_connectivity[5] = i * (n_elements_y * 3 + 2) + n_elements_y * 2 + 1 + j + 1;
            elem_connectivity[6] = i * (n_elements_y * 3 + 2) + j * 2 + 2;
            elem_connectivity[7] = i * (n_elements_y * 3 + 2) + j * 2 + 1;
            mesh.elements().add(ElementType::Q2D8N, elem_connectivity);
        }
    }
    
    return mesh;
}
}