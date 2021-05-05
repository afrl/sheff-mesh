#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief Creates a mesh with a single element with the global coordinate system matching the
/// parametric/natural coordinate system for Gauss quadrature.
Mesh get_unit_element(const ElementType& type);
}