#include "mesh/external_interfaces/fea_io/MeshInfo.h"

#include "mesh/algorithms/Search.h"
#include "mesh/external_interfaces/fea_io/TypeConversions.h"
#include "utility/MPIUtils.h"
#include "utility/Range.h"

#include <map>
#include <unordered_set>

namespace mesh
{
#ifdef USE_FEAIO

fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_mesh_info(
  const Mesh& mesh,
  const ElementTypeCategory& element_category_to_write,
  const std::map<ComponentIndex, std::string> component_names)
{
    // Create the meta info the mesh file
    fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> mesh_info;
    mesh_info.dimensionality      = 3;
    mesh_info.connectivity_scheme = fea_io::ConnectivityScheme::FlexMesh;
    if (mesh.contains_all_partitions())
    {
        mesh_info.total_num_nodes    = mesh.nodes().size();
        mesh_info.total_num_elements = mesh.elements().size();
        auto extent                  = find_extent(mesh);
        mesh_info.extent             = fea_io::Extent<FloatType>(extent.first, extent.second);
    }
    else
    {
        GlobalIndex n_owned_nodes = 0;
        for (auto& node : mesh.nodes())
        {
            if (node.parent_partition() == mesh.partition_id())
            {
                n_owned_nodes++;
            }
        }
        mesh_info.total_num_nodes = utility::MPIUtils::sum_across_ranks(n_owned_nodes);

        GlobalIndex n_owned_elements = 0;
        for (auto& element : mesh.elements())
        {
            if (element.parent_partition() == mesh.partition_id())
            {
                n_owned_elements++;
            }
        }
        mesh_info.total_num_elements = utility::MPIUtils::sum_across_ranks(n_owned_elements);

        auto extent                  = find_global_extent(mesh);
        mesh_info.extent             = fea_io::Extent<FloatType>(extent.first, extent.second);
    }

    auto n_components = find_global_max_component_id(mesh);
    for (auto component_id : utility::range(n_components))
    {
        // Ensure the orientation type for all the elements is the same
        LCSAssignmentType first_assignment_type = LCSAssignmentType::None;
        for (auto& element : mesh.elements())
        {
            if (element.component_id() == component_id)
            {
                first_assignment_type = element.local_coord_sys().base_type();
                break;
            }
        }
        for (auto& element : mesh.elements())
        {
            if (element.component_id() == component_id)
            {
                if (element.local_coord_sys().base_type() != first_assignment_type)
                {
                    throw _UTIL_MSG_EXCEPT(
                      "The elements within the same component must have the same "
                      "LCSAssignment to be written to the fea_io format.");
                }
            }
        }

        mesh_info.components.push_back({ component_id,
                                         "Component_" + std::to_string(component_id),
                                         convert_lcs_assignment(first_assignment_type),
                                         fea_io::LCSParameterization::Quaternion });
    }
    mesh_info.total_num_partitions = mesh.n_partitions();
    return mesh_info;
}

fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_partition_info(
  const Mesh& mesh,
  const fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType>& mesh_info)
{
    // Create the partition info for the mesh file
    fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> partition_info;

    // Figure out what each subcomponent contains in this partition
    std::vector<fea_io::SubcomponentInfo<LocalIndex>> subcomponent_infos;
    for (auto& component_info : mesh_info.components)
    {
        std::unordered_set<LocalIndex> nodes_in_subcomponent;
        LocalIndex n_elements_in_subcomponent = 0;
        for (auto& element : mesh.elements())
        {
            if (element.component_id() == component_info.id)
            {
                for (auto& node_index : element.connectivity())
                {
                    nodes_in_subcomponent.insert(node_index);
                }
                n_elements_in_subcomponent++;
            }
        }

        subcomponent_infos.push_back({ static_cast<LocalIndex>(component_info.id),
                                       static_cast<LocalIndex>(nodes_in_subcomponent.size()),
                                       n_elements_in_subcomponent });
    }

    LocalIndex total_num_connectivity_entries = 0;
    for (auto& element : mesh.elements())
    {
        total_num_connectivity_entries += static_cast<LocalIndex>(element.connectivity().size());
    }

    partition_info.partition_id                   = mesh.partition_id();
    partition_info.num_nodes                      = static_cast<LocalIndex>(mesh.nodes().size());
    partition_info.num_elements                   = static_cast<LocalIndex>(mesh.elements().size());
    partition_info.total_num_connectivity_entries = total_num_connectivity_entries;
    auto extent                                   = find_extent(mesh);
    partition_info.extent            = fea_io::Extent<FloatType>(extent.first, extent.second);
    partition_info.num_subcomponents = static_cast<LocalIndex>(subcomponent_infos.size());
    partition_info.subcomponents     = subcomponent_infos;

    return partition_info;
}

#endif
}