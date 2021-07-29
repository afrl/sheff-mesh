#pragma once

#include "mesh/Mesh.h"

#include <string>

namespace mesh
{
/// Reads an BSAM classic .in file and populates a mesh from the input file.
/// @param file_path - relative (to cwd) or absolute path to .in file
Mesh read_bsam_classic(const std::string& file_path);

/// Reads an BSAM classic .in file and populates a mesh for each cluster in the input file.
/// \param file_path - relative (to cwd) or absolute path to .in file
std::vector<Mesh> read_bsam_classic_by_cluster(const std::string& file_path);

/// Reads a .ele mesh file, which are used in BSAM classic, and populates a mesh.
/// @param file_path - relative (to cwd) or absolute path to .in file
Mesh read_ele(const std::string& file_path);
}