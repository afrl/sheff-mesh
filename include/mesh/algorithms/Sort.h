#pragma once

#include "mesh/Mesh.h"
#include "mesh/algorithms/Search.h"

namespace mesh
{
/// @brief Sorts a list of node indices by coordinate (x, then y, then z)
/// @param mesh - the mesh to which node_indices refers
/// @param node_indices - list of local node indices (modified)
void sort_node_indices_by_coord(const Mesh& mesh, std::vector<mesh::LocalIndex>& node_indices);

void sort_element_face_pairs_by_coord(const Mesh& mesh,
                                      std::vector<ElementFacePair>& element_face_pairs);
}