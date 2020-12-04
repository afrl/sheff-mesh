#include "Groups.h"

namespace mesh
{
void add_to_group(Mesh& mesh, const std::string group_name, NodeSelection& selection)
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

void add_to_group(Mesh& mesh, const std::string group_name, ElementSelection& selection)
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
}
