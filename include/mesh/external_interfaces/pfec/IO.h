#pragma once

#include "mesh/Mesh.h"
#include "mesh/common/ElementTypes.h"
#include "mesh/common/MeshAssociatedData.h"

#include <string>
#include <vector>

namespace mesh
{
/// TODO write doc string
std::pair<Mesh, MeshAssociatedData> read_pfec_all_parts(const std::string& meta_file_path,
  const int time_step_id,
  const bool convert_node_data_to_element_data = false);

/// @brief Reads a PFEC binary format, populates a mesh, and reads all available data.
/// @param meta_file_path - relative (to cwd) or absolute path to the .pfecx file
/// @param partition - the partition to read
///
/// The PFEC file format was developed by Keith Ballard for his dissertation work.  There are still
/// existing tools to create fiber/matrix, multiscale laminated composite, and 3D textile models
/// that use the format.
std::pair<Mesh, MeshAssociatedData> read_pfec(const std::string& meta_file_path,
                                              const PartitionIndex partition = 0,
                                              const int time_step_id         = 0);

/// @brief Reads a PFEC binary format and populates a mesh.
/// @param meta_file_path - relative (to cwd) or absolute path to the .pfecx file
/// @param partition - the partition to read
/// @param time_step_id - the timestep to read (starts at 0)
///
/// The PFEC file format was developed by Keith Ballard for his dissertation work.  There are still
/// existing tools to create fiber/matrix, multiscale laminated composite, and 3D textile models
/// that use the format.
Mesh read_pfec_mesh(const std::string& meta_file_path,
                    const PartitionIndex partition = 0,
                    const int time_step_id         = 0);

/// @brief Reads all data from a file in the PFEC binary format.
/// @param meta_file_path - relative (to cwd) or absolute path to the .pfecx file
/// @param partition - the partition to read
/// @param time_step_id - the timestep to read (starts at 0)
///
/// The PFEC file format was developed by Keith Ballard for his dissertation work.  There are still
/// existing tools to create fiber/matrix, multiscale laminated composite, and 3D textile models
/// that use the format.
MeshAssociatedData read_pfec_data(const std::string& meta_file_path,
                                  const PartitionIndex partition = 0,
                                  const int time_step_id         = 0);

MeshAssociatedData read_pfec_data(
  const std::string& meta_file_path,
  const utility::Slice<const std::string, std::size_t>& node_data_names_to_read,
  const utility::Slice<const std::string, std::size_t>& element_data_names_to_read,
  const PartitionIndex partition = 0,
  const int time_step_id         = 0);
}