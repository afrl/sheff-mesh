#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief Convert a mesh of linear elements to quadratic equivalents
void convert_to_quadratic(Mesh& mesh);

/// @brief Given a list of nodes that should be combined into another node, performs the combination
/// and removes the nodes that were combined.
/// @param mesh - mesh object to modify
/// @param combinations - a list of pairs in the form { index of node to be removed, index of node
///        to replace it }.  Note that the combinations must be unique and non-circular (these
///        requirements are not checked inside in the algorithm).
void combine_nodes(Mesh& mesh, std::vector<std::pair<LocalIndex, LocalIndex>> combinations);

/// @brief Removes duplicate nodes from the mesh based on the nodes' positions.
void remove_duplicate_nodes(Mesh& mesh, const double tol = 1e-12);
}