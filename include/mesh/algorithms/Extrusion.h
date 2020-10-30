#pragma once

#include "mesh/Mesh.h"

namespace mesh
{

/// Extrude a 2D mesh along z direction
Mesh extrude_2D_mesh(const Mesh& mesh, const std::vector<double>& z_coordinates);

/// Extrude a 2D mesh along z direction with some offsets in x and y
Mesh extrude_2D_meshalong_path(const Mesh& mesh,
                               const std::vector<double>& z_coordinates,
                               const std::vector<double>& x_offsets,
                               const std::vector<double>& y_offsets);
}