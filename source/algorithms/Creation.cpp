#include "mesh/algorithms/Creation.h"

#include "mesh/algorithms/Extrusion.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/algorithms/Search.h"
#include "mesh/algorithms/Transformation.h"
#include "simple_math/Utility.h"

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Helpers

/// @brief Creates a small block of mixed elements with all quarilateral faces on the outside.
///
/// It is composed of 10 pyramid, 4 tetrahedron, and 2 wedge elements.
///             ^ z
///             |
///            4O-----------O7
///            /|          /|
///           / |         / | (block of pyramids and 2 tets)
///         5O-----------O6 |
///         /|  |    X8 /|  |
///        / | 0X------/-|--O3 --> y
///     11O-----------O12| /
///      /|  |/      /|  |/   (middle layer is wedges)
///     / | 1X------/-|--O2
///  15O-----------O16| /
///    |  |/    X17|  |/
///    | 9X--------|--O10
///    | /         | /
///    |/          |/  (block of pyramids and 2 tets)
///  13O-----------O14
///      /
///     x
/// Diagonal connecting 0 to 6
/// Diagonal connecting 9 to 12
static Mesh make_mixed_element_block()
{
    Mesh block;
    block.nodes().add({ 0, 0, 0 });       // 0
    block.nodes().add({ 1, 0, 0 });       // 1
    block.nodes().add({ 1, 1, 0 });       // 2
    block.nodes().add({ 0, 1, 0 });       // 3
    block.nodes().add({ 0, 0, 1 });       // 4
    block.nodes().add({ 1, 0, 1 });       // 5
    block.nodes().add({ 1, 1, 1 });       // 6
    block.nodes().add({ 0, 1, 1 });       // 7
    block.nodes().add({ 0.5, 0.5, 0.5 }); // 8
    block.nodes().add({ 2, 0, 0 });       // 9
    block.nodes().add({ 2, 1, 0 });       // 10
    block.nodes().add({ 2, 0, 1 });       // 11
    block.nodes().add({ 2, 1, 1 });       // 12
    block.nodes().add({ 3, 0, 0 });       // 13
    block.nodes().add({ 3, 1, 0 });       // 14
    block.nodes().add({ 3, 0, 1 });       // 15
    block.nodes().add({ 3, 1, 1 });       // 16
    block.nodes().add({ 2.5, 0.5, 0.5 }); // 17

    block.elements().add(ElementType::P3D5N, { 0, 1, 2, 3, 8 });
    block.elements().add(ElementType::P3D5N, { 2, 6, 7, 3, 8 });
    block.elements().add(ElementType::P3D5N, { 6, 5, 4, 7, 8 });
    block.elements().add(ElementType::P3D5N, { 0, 4, 5, 1, 8 });
    block.elements().add(ElementType::P3D5N, { 0, 3, 7, 4, 8 });
    block.elements().add(ElementType::T3D4N, { 1, 5, 6, 8 });
    block.elements().add(ElementType::T3D4N, { 1, 6, 2, 8 });
    block.elements().add(ElementType::W3D6N, { 1, 6, 5, 9, 12, 11 });
    block.elements().add(ElementType::W3D6N, { 1, 2, 6, 9, 10, 12 });
    block.elements().add(ElementType::T3D4N, { 9, 10, 12, 17 });
    block.elements().add(ElementType::T3D4N, { 12, 11, 9, 17 });
    block.elements().add(ElementType::P3D5N, { 9, 11, 15, 13, 17 });
    block.elements().add(ElementType::P3D5N, { 13, 15, 16, 14, 17 });
    block.elements().add(ElementType::P3D5N, { 10, 14, 16, 12, 17 });
    block.elements().add(ElementType::P3D5N, { 9, 13, 14, 10, 17 });
    block.elements().add(ElementType::P3D5N, { 11, 12, 16, 15, 17 });

    write_vtk(block, "test.vtk", false);

    return block;
};

/// @brief Creates a 3x3x3 block of mixed elements
static Mesh make_super_mixed_element_block()
{
    auto hex_3_stack   = make_uniform_linear_grid({ 0, 0, 0 }, { 3, 1, 1 }, { 3, 1, 1 });
    auto mixed_3_stack = make_mixed_element_block();

    Mesh mixed_grid;
    // Layer 0
    mixed_grid.copy(mixed_3_stack);
    mixed_grid.copy(translate(hex_3_stack, { 0, 1, 0 }));
    mixed_grid.copy(translate(mixed_3_stack, { 0, 2, 0 }));
    // Layer 1
    mixed_grid.copy(translate(hex_3_stack, { 0, -1, 1 }));
    mixed_grid.copy(translate(mixed_3_stack, { 0, -1, 1 }));
    mixed_grid.copy(translate(hex_3_stack, { 0, 2, 0 }));
    // Layer 2
    mixed_grid.copy(translate(mixed_3_stack, { 0, -1, 1 }));
    mixed_grid.copy(translate(hex_3_stack, { 0, -1, 1 }));
    mixed_grid.copy(translate(mixed_3_stack, { 0, 2, 0 }));

    remove_duplicate_nodes(mixed_grid);
    write_vtk(mixed_grid, "test2.vtk", false);

    return mixed_grid;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions for signatures in header

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

Mesh make_mixed_element_grid(const Point& min,
                             const Point& max,
                             const std::array<uint32_t, 3>& n_blocks_along_axes)
{
    Mesh mixed_grid;
    auto mixed_element_block = make_super_mixed_element_block();
    auto i_vec               = simple_math::Vector3D(3, 0, 0);
    auto j_vec               = simple_math::Vector3D(0, 3, 0);
    auto k_vec               = simple_math::Vector3D(0, 0, 3);
    for (auto i : utility::range(n_blocks_along_axes[0]))
    {
        for (auto j : utility::range(n_blocks_along_axes[1]))
        {
            for (auto k : utility::range(n_blocks_along_axes[2]))
            {
                auto t = i * i_vec + j * j_vec + k * k_vec;
                Mesh tmp_block(mixed_element_block);
                mixed_grid.copy(translate(tmp_block, simple_math::to_point(t)));
            }
        }
    }

    remove_duplicate_nodes(mixed_grid);

    for (auto node : mixed_grid.nodes())
    {
        node.global_index() = node.local_index();
    }
    for (auto element : mixed_grid.elements())
    {
        element.global_index() = element.local_index();
    }

    auto extent = find_extent(mixed_grid);
    std::array<double, 3> scale_factors;
    for (auto i : utility::range(3))
    {
        scale_factors[i] = (max[i] - min[i]) / (extent.second[i] - extent.first[i]);
    }
    scale(mixed_grid, scale_factors);
    translate(mixed_grid, min);

    return mixed_grid;
}
}
