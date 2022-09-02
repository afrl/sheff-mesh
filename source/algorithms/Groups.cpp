#include "mesh/algorithms/Groups.h"

#include "mesh/algorithms/Search.h"
#include "simple_math/Utility.h"

namespace mesh
{
void create_groups_for_nodes_on_extent(Mesh& mesh,
                                       const std::string& name_prefix,
                                       const double& rel_tol)
{
    auto extent      = find_extent(mesh);
    auto& min_coord  = extent.first;
    auto& max_coord  = extent.second;
    const double tol = rel_tol * simple_math::distance(min_coord, max_coord);

    std::string prefix = "";
    if (name_prefix.size() > 0)
    {
        prefix = name_prefix + "_";
    }

    add_to_group(mesh,
                 prefix + "neg_x",
                 NodeSelection([val = min_coord[0], t = tol](const NodeConstAccessor& n) -> bool {
                     return fabs(n.x() - val) < t;
                 }));
    add_to_group(mesh,
                 prefix + "pos_x",
                 NodeSelection([val = max_coord[0], t = tol](const NodeConstAccessor& n) -> bool {
                     return fabs(n.x() - val) < t;
                 }));
    add_to_group(mesh,
                 prefix + "neg_y",
                 NodeSelection([val = min_coord[1], t = tol](const NodeConstAccessor& n) -> bool {
                     return fabs(n.y() - val) < t;
                 }));
    add_to_group(mesh,
                 prefix + "pos_y",
                 NodeSelection([val = max_coord[1], t = tol](const NodeConstAccessor& n) -> bool {
                     return fabs(n.y() - val) < t;
                 }));
    add_to_group(mesh,
                 prefix + "neg_z",
                 NodeSelection([val = min_coord[2], t = tol](const NodeConstAccessor& n) -> bool {
                     return fabs(n.z() - val) < t;
                 }));
    add_to_group(mesh,
                 prefix + "pos_z",
                 NodeSelection([val = max_coord[2], t = tol](const NodeConstAccessor& n) -> bool {
                     return fabs(n.z() - val) < t;
                 }));
}
void add_to_group(Mesh& mesh, const std::string group_name, NodeSelection selection)
{
    auto& group = mesh.nodes().group(group_name);
    selection.run_selection(mesh);
    auto selected_indices = selection.get_selected_indices();
    // Add any existing nodes into the set
    for (auto& i : group.node_indices)
    {
        selected_indices.insert(i);
    }
    group.node_indices = std::vector<LocalIndex>(selected_indices.begin(), selected_indices.end());
}

void add_to_group(Mesh& mesh, const std::string group_name, ElementSelection selection)
{
    auto& group = mesh.elements().group(group_name);
    selection.run_selection(mesh);
    auto selected_indices = selection.get_selected_indices();
    // Add any existing nodes into the set
    for (auto& i : group.element_indices)
    {
        selected_indices.insert(i);
    }
    group.element_indices =
      std::vector<LocalIndex>(selected_indices.begin(), selected_indices.end());
}

void add_nodes_to_group(Mesh& mesh,
                        const std::string group_name,
                        std::vector<mesh::LocalIndex>& local_ids_to_add)
{
    auto& group = mesh.nodes().group(group_name);
    std::set<mesh::LocalIndex> ids_in_group(group.node_indices.begin(), group.node_indices.end());
    for (auto& i : local_ids_to_add)
    {
        ids_in_group.insert(i);
    }
    group.node_indices = std::vector<mesh::LocalIndex>(ids_in_group.begin(), ids_in_group.end());
}

void add_elements_to_group(Mesh& mesh,
                           const std::string group_name,
                           std::vector<mesh::LocalIndex>& local_ids_to_add)
{
    auto& group = mesh.elements().group(group_name);
    std::set<mesh::LocalIndex> ids_in_group(group.element_indices.begin(),
                                            group.element_indices.end());
    for (auto& i : local_ids_to_add)
    {
        ids_in_group.insert(i);
    }
    group.element_indices = std::vector<mesh::LocalIndex>(ids_in_group.begin(), ids_in_group.end());
}
}
