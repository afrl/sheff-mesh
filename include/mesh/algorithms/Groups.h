#pragma once

#include "mesh/common/Selection.h"

#include <vector>

namespace mesh
{
/// @brief Creates groups for nodes along the extent (coincident with the bounding box).
///
/// Creates the following groups: "<prefix_>neg_x" "<prefix_>pos_x" "<prefix_>neg_y"
///                               "<prefix_>pos_y" "<prefix_>neg_z" "<prefix_>pos_z"
/// @param name_prefix - a prefix for the group names (a "_" will be appended to the prefix to form
///                      the name if the prefix is not empty)
/// @param rel_tol - tolerance to use when finding nodes on the extent (relative to diagonal
/// 	   of the bounding box of the mesh)
void create_groups_for_nodes_on_extent(Mesh& mesh,
                                       const std::string& name_prefix = "",
                                       const double& rel_tol          = 1e-12);

/// @brief Adds the nodes to a group based on a selection object.
/// @param group_name - the name of the group to add the selected nodes to (does not require the
///        group to exist before this call)
void add_to_group(Mesh& mesh, const std::string group_name, NodeSelection selection);

/// @brief Adds the elements to a group based on a selection object.
/// @param group_name - the name of the group to add the selected elements to (does not require the
///        group to exist before this call)
void add_to_group(Mesh& mesh, const std::string group_name, ElementSelection selection);

/// @brief Adds the list of nodes to a group (ensuring no duplicates)
/// @param group_name - the name of the group to add to (does not require the group to exist before
///        this call)
void add_nodes_to_group(Mesh& mesh,
                        const std::string group_name,
                        std::vector<mesh::LocalIndex>& local_ids_to_add);

/// @brief Adds the list of elements to a group (ensuring no duplicates)
/// @param group_name - the name of the group to add to (does not require the group to exist before
///        this call)
void add_elements_to_group(Mesh& mesh,
                           const std::string group_name,
                           std::vector<mesh::LocalIndex>& local_ids_to_add);
}
