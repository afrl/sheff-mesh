#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief Convert a mesh of linear elements to quadratic equivalents
void convert_to_quadratic(Mesh& mesh);

/// @brief Removes duplicate nodes from the mesh based on position.
//void remove_duplicate_nodes(Mesh& mesh);
}