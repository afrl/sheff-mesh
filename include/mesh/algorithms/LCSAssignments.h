#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief Assigns the local coordinate system for a set/group of elements.
///
/// @param element_set_name - name of element set/group to assign the LCS (error will not be throw
///        if group does not exist)
/// @param x_vector - local x-axis expressed in terms of the global coordinate system.
/// @param y_vector - local y-axis expressed in terms of the global coordinate system.
void assign_lcs(Mesh& m,
                const std::string& element_set_name,
                const simple_math::Point& x_vector,
                const simple_math::Point& y_vector);
}