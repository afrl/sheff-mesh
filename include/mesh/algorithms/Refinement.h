#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief Corefines the mesh (splitting every edge).
/// 
/// Refinement occurs as follows:
/// * Each tetrahedron is divided into 5 new tetrahedra.
/// TODO explain others
void refine_uniform(Mesh& mesh);
}