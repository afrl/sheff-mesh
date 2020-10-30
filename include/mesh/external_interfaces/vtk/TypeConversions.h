#pragma once

#include "mesh/common/ElementTypes.h"

#include <vector>

namespace mesh
{
/// Returns the integer used to indicate VTK cell types for a given ElementType
int get_vtk_element_type(const ElementType& geom_type);
/// Returns the permuation vector to go from the internal connectivity to the
/// connectivity used in VTK.  All types are not mappable to a VTK element,
/// in those cases, it will throw.
std::vector<ConnIndex> get_vtk_permutation(const ElementType& geom_type, const ConnIndex& n_nodes);
}