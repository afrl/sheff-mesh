#pragma once

#include "fea_io/mesh/MeshInfoStructs.h"
#include "mesh/Defines.h"
#include "mesh/Mesh.h"
#include "mesh/common/ElementTypes.h"

#include <string>
#include <vector>

namespace mesh
{

/// @brief A helper that creates a fea_io::MeshInfo object from information about the mesh.
///
/// This is used in write_fea_io and write_distributed_fea_io and when exporting data.
///
/// @param mesh - a mesh instance containing either a single partition or the global mesh.
/// @param component_names - a map from the component ID to the name of the component.  If none is
///        provided then components will automatically be called "component_<id>".
fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_mesh_info(
  const Mesh& mesh,
  const ElementTypeCategory& element_category_to_write = ElementTypeCategory::ANY,
  const std::map<ComponentIndex, std::string> component_name_map =
    std::map<ComponentIndex, std::string>());

/// This is used in write_fea_io and write_distributed_fea_io and when exporting data.
/// Computes the partition information needed for writing in parallel with fea_io
fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_partition_info(
  const Mesh& mesh,
  const fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType>& mesh_info);

}