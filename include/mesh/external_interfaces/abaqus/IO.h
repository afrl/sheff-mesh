#pragma once

#include "mesh/Mesh.h"

#include <string>

namespace mesh
{
/// Reads an Abaqus .inp file (limited to simple .inp files right now) and populates a mesh.
/// \param file_path - relative (to cwd) or absolute path to .inp file
Mesh read_abaqus(const std::string& file_path);

/* To be implemented
/// Write the mesh to Abaqus .inp file
/// @param mesh - mesh to write (note: if distributed, this only writes the part of
/// the mesh in this partition, so be sure to provide a unique path for each partition)
/// @param file_path - relative (to cwd) or absolute path (include the .inp extension)
void write_abaqus(Mesh& mesh,
                  const std::string& file_path);
*/

}