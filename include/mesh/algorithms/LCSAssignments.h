#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
void assign_lcs(Mesh& m,
                const std::string& element_set_name,
                const simple_math::Point& x_vector,
                const simple_math::Point& y_vector);
}