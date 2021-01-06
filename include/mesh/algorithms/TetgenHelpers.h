#pragma once

#include "mesh/Mesh.h"
#include "second_party_mesh/tetgen/tetgen.h"

namespace mesh
{
/// Creates a mesh that describes the input to tetgen.
Mesh get_mesh_from_tetgen_input(tetgen::tetgenio& in);

/// Creates a mesh that describes the output from tetgen.
Mesh get_mesh_from_tetgen_output(tetgen::tetgenio& out);
}