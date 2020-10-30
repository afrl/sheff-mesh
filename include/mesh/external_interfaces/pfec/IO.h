#pragma once

#include "mesh/Mesh.h"
#include "mesh/common/ElementTypes.h"

#include <string>
#include <vector>

namespace mesh
{
/// @brief Reads a PFEC binary format and populates a mesh.
/// @param meta_file_path - relative (to cwd) or absolute path to the .pfecx file
/// @param partition - the partition to read
/// The PFEC file format was developed by Keith Ballard for his dissertation work.  There are still
/// existing tools to create fiber/matrix, multiscale laminated composite, and 3D textile models
/// that use the format.
Mesh read_pfec(const std::string& meta_file_path, const PartitionIndex partition = 0);
}