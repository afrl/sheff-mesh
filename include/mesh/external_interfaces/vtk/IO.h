#pragma once

#include "mesh/Mesh.h"
#include "mesh/common/ElementTypes.h"

#include <string>
#include <vector>

namespace mesh
{
/// Returns the integer used to indicate VTK cell types for a given ElementType
int get_vtk_element_type(const ElementType& geom_type);
/// Returns the permuation vector to go from the internal connectivity to the
/// connectivity used in VTK.  All types are not mappable to a VTK element,
/// in those cases, it will throw.
std::vector<ConnIndex> get_vtk_permutation(const ElementType& geom_type, const ConnIndex& n_nodes);

/// Write the mesh to VTK file as an unstructured grid
/// \param mesh - mesh to write to the VTK file (note: if distributed, this only writes the part of
/// the mesh of this partition, so be sure to provide a unique path for each partition)
/// \param file_path - relative (to cwd) or absolute path (include the .vtk extension)
/// \param write_lcs - indicates whether to write each axis of the lcs as vectors in the data
/// \param point_data - additional point data to include
void write_vtk(
  Mesh& mesh,
  const std::string& file_path,
  const bool write_lcs                                      = true,
  const std::vector<double>& point_data                     = std::vector<double>(),
  const std::vector<double>& element_data                   = std::vector<double>(),
  const std::vector<unsigned long long>& point_attributes   = std::vector<unsigned long long>(),
  const std::vector<unsigned long long>& element_attributes = std::vector<unsigned long long>());
}