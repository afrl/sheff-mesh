#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief Performs a rigid body translation of a mesh.
/// @return Reference to the given mesh to allow operations to be included on the same line.
Mesh& translate(Mesh& mesh, const simple_math::Point& translation_vector);

/// @brief Performs a linear scaling along all three of a mesh about a point.
/// @return Reference to the given mesh to allow operations to be included on the same line.
Mesh& scale(Mesh& mesh, const double& scale_factor, const simple_math::Point& origin = { 0, 0, 0 });

/// @brief Performs a scaling along each coordinate axis a mesh about a point.
/// @return Reference to the given mesh to allow operations to be included on the same line.
Mesh& scale(Mesh& mesh,
            const std::array<double, 3>& scale_factors,
            const simple_math::Point& origin = { 0, 0, 0 });
}