#pragma once

#include "mesh/common/Selection.h"

namespace mesh
{
/// @brief Adds the nodes to a group based on a selection object.
/// @param group_name - the name of the group to add the selected nodes to (does not require the
///        group to exist before this call)
void add_to_group(Mesh& mesh, const std::string group_name, NodeSelection& selection);

/// @brief Adds the elements to a group based on a selection object.
/// @param group_name - the name of the group to add the selected elements to (does not require the
///        group to exist before this call)
void add_to_group(Mesh& mesh, const std::string group_name, ElementSelection& selection);
}