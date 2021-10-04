#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// Generate a 2D grid composed of linear quadrilateral elements.
/// @param x_coordinates - coordinates of the corners of the elements along x-axis
/// @param y_coordinates - coordinates of the corners of the elements along y-axis
Mesh make_linear_grid(const std::vector<double>& x_coordinates,
                      const std::vector<double>& y_coordinates);

/// Generate a 3D grid composed of linear quadrilateral elements.
/// @param x_coordinates - coordinates of the corners of the elements along x-axis
/// @param y_coordinates - coordinates of the corners of the elements along y-axis
/// @param z_coordinates - coordinates of the corners of the elements along z-axis
Mesh make_linear_grid(const std::vector<double>& x_coordinates,
                      const std::vector<double>& y_coordinates,
                      const std::vector<double>& z_coordinates);

/// Generates a 3D uniform grid in the bounding box with linear hexhedron elements
Mesh make_uniform_linear_grid(const Point& min,
                              const Point& max,
                              const std::array<int, 3>& n_elements);

/// Generates a 3D uniform grid in the bounding box with linear tetrahedron elements
Mesh make_uniform_linear_tet_grid(const Point& min,
                                  const Point& max,
                                  const std::array<int, 3>& n_edges_along_axis);

/// Generate a 2D rectangular mesh with a uniform grid
/// @param x_coordinates - coordinates of the corners of the grid cells along x-axis
/// @param y_coordinates - coordinates of the corners of the grid cells along y-axis
Mesh make_quadratic_grid(const std::vector<double>& x_coordinates,
                         const std::vector<double>& y_coordinates);

/// @brief Generates a 3D uniform grid composed of a mix of elements
/// @param n_blocks_along_axes - an integer indicating how many 3x3x3 mixed element units to tile
///        along each axis.
Mesh make_mixed_element_grid(const Point& min,
                             const Point& max,
                             const std::array<uint32_t, 3>& n_blocks_along_axes);
}