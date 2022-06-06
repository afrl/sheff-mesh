#include "mesh/external_interfaces/fea_io/IO.h"

#ifdef USE_FEAIO

#include "fea_io/common/Filesystem.h"
#include "fea_io/data/DataWriter.h"
#include "fea_io/group/GroupReader.h"
#include "fea_io/group/GroupWriter.h"
#include "fea_io/mesh/MeshParWriter.h"
#include "fea_io/mesh/MeshReader.h"
#include "fea_io/mesh/MeshWriter.h"
#include "fea_io/meta/MetaReader.h"
#include "fea_io/meta/MetaWriter.h"
#include "mesh/algorithms/Partitioning.h"
#include "mesh/external_interfaces/fea_io/MeshInfo.h"
#include "mesh/external_interfaces/fea_io/TypeConversions.h"
#include "mesh/relations/ElementAdjacencyList.h"
#include "utility/MPIUtils.h"
#include "utility/Range.h"

#include <unordered_set>

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions

static std::string write_meta_file(const std::string& base_directory,
                                   const std::string& mesh_name,
                                   const std::string meta_file_name,
                                   bool include_groups,
                                   bool include_data)
{
    if (meta_file_name != "")
    {
        fea_io::MetaWriter meta_writer;
        std::string meta_filepath = (fs::path(base_directory) / meta_file_name).string();
        meta_writer.open(meta_filepath);
        meta_writer.get_analysis_meta().name = meta_file_name;
        meta_writer.get_analysis_meta().levels_of_detail_names.push_back("original");

        fea_io::LevelOfDetail lod;
        lod.level                    = 0;
        lod.relationship_with_parent = fea_io::LODParentRelation::None;
        fea_io::MeshObject mesh_object;
        mesh_object.type      = fea_io::MeshObjectType::Volume;
        mesh_object.mesh_path = mesh_name;
        if (include_groups)
        {
            mesh_object.group_paths.push_back(mesh_name);
        }
        if (include_data)
        {
            mesh_object.data_paths.push_back(mesh_name);
        }
        lod.mesh_objects.push_back(mesh_object);

        meta_writer.get_analysis_meta().input_info.levels_of_detail.push_back(lod);
    }

    std::string file_stem = (fs::path(base_directory) / mesh_name).string();
    return file_stem;
}

static Point min_coord(const NodeContainer& nodes)
{
    if (nodes.size() > 0)
    {
        Point min = nodes[0].point();
        for (auto& node : nodes)
        {
            auto p = node.point();
            for (auto i = 0; i < 3; i++)
            {
                if (p[i] < min[i])
                {
                    min[i] = p[i];
                }
            }
        }
        return min;
    }
    else
    {
        return { 0, 0, 0 };
    }
}

static std::array<FloatType, 3> max_coord(const NodeContainer& nodes)
{
    if (nodes.size() > 0)
    {
        Point max = nodes[0].point();
        for (auto& node : nodes)
        {
            auto p = node.point();
            for (auto i = 0; i < 3; i++)
            {
                if (p[i] > max[i])
                {
                    max[i] = p[i];
                }
            }
        }
        return max;
    }
    else
    {
        return { 0, 0, 0 };
    }
}

/// Extract group information for writing to the fea_io format
std::vector<fea_io::GroupInfo> get_group_infos(const Mesh& mesh)
{
    std::vector<fea_io::GroupInfo> infos;
    for (auto& group : mesh.nodes().groups())
    {
        fea_io::GroupInfo info;
        info.initialize<LocalIndex>(group.name, fea_io::GroupType::Node, {});
        infos.push_back(info);
    }
    for (auto& group : mesh.elements().groups())
    {
        fea_io::GroupInfo info;
        info.initialize<LocalIndex>(group.name, fea_io::GroupType::Element, {});
        infos.push_back(info);
    }
    return infos;
}

/// Gets the fea_io::Partition object that contains all of the data needed to write the
/// mesh for this parititon.
/// @attention nodes_in_partition is non-const because their local numbering will be updated.
static fea_io::Partition<GlobalIndex, LocalIndex, FloatType> get_fea_io_partition(
  const Mesh& mesh,
  const PartitionIndex partition_id,
  const std::vector<LocalIndex>& nodes_in_partition,
  const std::vector<LocalIndex>& elements_in_partition,
  const fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType>& mesh_info)
{
    fea_io::Partition<GlobalIndex, LocalIndex, FloatType> partition;
    partition.partition_id             = partition_id;
    const auto n_nodes_in_partition    = nodes_in_partition.size();
    const auto n_elements_in_partition = elements_in_partition.size();

    // Note: nodes that are not in the partition get the max limit
    std::vector<LocalIndex> mesh_to_partition_node_map(mesh.nodes().size(),
                                                       std::numeric_limits<LocalIndex>::max());
    for (auto i : utility::indices(nodes_in_partition))
    {
        mesh_to_partition_node_map[nodes_in_partition[i]] = static_cast<LocalIndex>(i);
    }

    // Populate node information
    LocalIndex local_num = 0;
    partition.node_global_ids.reserve(n_nodes_in_partition);
    partition.node_parent_partitions.reserve(n_nodes_in_partition);
    partition.node_parent_partitions.reserve(n_nodes_in_partition);
    partition.node_coords.reserve(3 * n_nodes_in_partition);
    for (auto& n_local_id : nodes_in_partition)
    {
        auto node = mesh.nodes()[n_local_id];

        partition.node_global_ids.push_back(node.global_index());
        partition.node_parent_partitions.push_back(node.parent_partition());
        for (auto i = 0; i < 3; i++)
        {
            partition.node_coords.push_back(node.point()[i]);
        }
        local_num++;
    }

    // Populate element information
    std::vector<LocalIndex> local_connectivity;
    local_connectivity.reserve(40);
    partition.element_global_ids.reserve(n_elements_in_partition);
    partition.element_parent_partitions.reserve(n_elements_in_partition);
    partition.element_types.reserve(n_elements_in_partition);
    partition.element_connectivities.reserve(n_elements_in_partition);
    for (auto& e_local_id : elements_in_partition)
    {
        auto element = mesh.elements()[e_local_id];

        partition.element_global_ids.push_back(element.global_index());
        partition.element_parent_partitions.push_back(element.parent_partition());
        partition.element_types.push_back(convert_element_type(element.type()));
        // This method returns the local node numbers for most elements, but returns
        // a face stream for polyhedra.
        std::vector<LocalIndex> local_connectivity;
        local_connectivity.reserve(element.connectivity().size());
        // Transform the connectivity to the partition numbering
        for (auto& v : element.connectivity())
        {
            assert(mesh_to_partition_node_map[v] != std::numeric_limits<LocalIndex>::max());
            local_connectivity.push_back(mesh_to_partition_node_map[v]);
        }
        assert(std::all_of(local_connectivity.begin(),
                           local_connectivity.end(),
                           [node_size = local_num](const LocalIndex& i) -> bool
                           { return i < node_size; }));

        partition.element_connectivities.push_back(local_connectivity);
    }

    // Figure out what each subcomponent contains in this partition
    std::vector<fea_io::Subcomponent<LocalIndex, FloatType>> subcomponents;
    for (auto& component_info : mesh_info.components)
    {
        fea_io::Subcomponent<LocalIndex, FloatType> subcomponent;
        subcomponent.component_id = component_info.id;

        // The subcomponent needs to store the list of nodes and element indices that appear in the
        // subcomponent based on the numbering of nodes_in_partition and elements_in_partition,
        // so a map is needed from the mesh-based index to the partition-based index for ndoes.
        {
            std::vector<LocalIndex> mesh_to_partition_node_index_map(
              mesh.nodes().size(),
              std::numeric_limits<LocalIndex>::max());
            for (auto i : utility::indices(nodes_in_partition))
            {
                mesh_to_partition_node_index_map[nodes_in_partition[i]] =
                  static_cast<LocalIndex>(i);
            }

            std::vector<LocalIndex> mesh_to_partition_element_index_map(
              mesh.elements().size(),
              std::numeric_limits<LocalIndex>::max());
            for (auto i : utility::indices(elements_in_partition))
            {
                mesh_to_partition_element_index_map[elements_in_partition[i]] =
                  static_cast<LocalIndex>(i);
            }

            // Note: indices are in the partitiong numbering (not neccessarily the total mesh)
            std::set<LocalIndex> nodes_in_subcomponent_set;
            for (auto i : utility::indices(elements_in_partition))
            {
                auto& e_local_id = elements_in_partition[i];
                auto element     = mesh.elements()[e_local_id];

                if (element.component_id() == component_info.id)
                {
                    subcomponent.elements_in_subcomponent.push_back(
                      mesh_to_partition_element_index_map[element.local_index()]);
                    for (auto& node_index : element.connectivity())
                    {
                        nodes_in_subcomponent_set.insert(
                          mesh_to_partition_node_index_map[node_index]);
                    }
                }
            }
            subcomponent.nodes_in_subcomponent.reserve(nodes_in_subcomponent_set.size());
            for (auto& i : nodes_in_subcomponent_set)
            {
                subcomponent.nodes_in_subcomponent.push_back(i);
            }
        }

        //  Note: this reserve is an estimate, but may be insufficient.
        subcomponent.lcs_parameters_in_subcomponent.reserve(
          4 * subcomponent.elements_in_subcomponent.size());
        for (auto& j : subcomponent.elements_in_subcomponent)
        {
            // Recall that the indices in subcomponent.elements_in_subcomponent correspond to
            // elements_in_partition.  So, be sure to map the index back to the mesh-based index
            // when accessing an element.
            auto& lcs = mesh.elements()[elements_in_partition[j]].local_coord_sys();
            for (auto i = 0; i < lcs.size(); i++)
            {
                const Quaternion& quatern = lcs[i];
                subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.w());
                subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.x());
                subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.y());
                subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.z());
            }
        }

        // Only include subcomponents for components that actually appear in the partition.
        if (!subcomponent.elements_in_subcomponent.empty())
        {
            subcomponents.push_back(subcomponent);
        }
    }
    partition.subcomponents = subcomponents;

    return partition;
}

void initialize_data_writer(
  const MeshAssociatedData& data,
  const std::string& file_stem,
  const fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType>& mesh_info,
  const fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType>& partition_info,
  fea_io::DataWriter<GlobalIndex, LocalIndex, FloatType>& data_writer)
{
    GlobalIndex n_quads_in_partition = 0;
    if (!data.get<DataFieldType::QUAD, double>().empty())
    {
        n_quads_in_partition = data.get<DataFieldType::QUAD, double>().front().data().rows();
    }

    std::vector<fea_io::ArrayInfo> data_infos;
    for (auto& d : data.get<DataFieldType::NODE, double>())
    {
        fea_io::ArrayInfo info;
        info.initialize<double>(d.name, fea_io::ArrayType::Node, d.data().cols(), d.component_names);
        data_infos.push_back(info);
    }
    for (auto& d : data.get<DataFieldType::ELEMENT, double>())
    {
        fea_io::ArrayInfo info;
        info.initialize<double>(d.name,
                                fea_io::ArrayType::Element,
                                d.data().cols(),
                                d.component_names);
        data_infos.push_back(info);
    }
    for (auto& d : data.get<DataFieldType::NODE_BY_COMPONENT, double>())
    {
        fea_io::ArrayInfo info;
        info.initialize<double>(d.name,
                                fea_io::ArrayType::NodeByComponent,
                                d.data().cols(),
                                d.component_names);
        data_infos.push_back(info);
    }
    for (auto& d : data.get<DataFieldType::NODE_BY_ELEMENT, double>())
    {
        fea_io::ArrayInfo info;
        info.initialize<double>(d.name,
                                fea_io::ArrayType::NodeByElement,
                                d.data().cols(),
                                d.component_names);
        data_infos.push_back(info);
    }
    for (auto& d : data.get<DataFieldType::QUAD, double>())
    {
        fea_io::ArrayInfo info;
        info.initialize<double>(d.name,
                                fea_io::ArrayType::QuadPointByElement,
                                d.data().cols(),
                                d.component_names);
        data_infos.push_back(info);
    }

    if (data_infos.size() > 0)
    {
        data_writer.initialize(file_stem,
                               mesh_info,
                               partition_info,
                               n_quads_in_partition,
                               data_infos);
    }
}

static void write_serial_single_partition(const Mesh& mesh,
                                          const MeshAssociatedData& data,
                                          const std::string file_stem)
{
    // This method treats the mesh as a single partition.
    const PartitionIndex n_partitions = 1;
    const PartitionIndex partition_id = 0;

    fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer;
    auto mesh_info = get_fea_io_mesh_info(mesh);
    mesh_writer.initialize(file_stem, mesh_info);

    // Get the nodes and elements that should appear in the partition corresponding to
    // partition_id, which is all of them.
    std::vector<LocalIndex> node_indices(mesh.nodes().size());
    std::iota(node_indices.begin(), node_indices.end(), 0);
    std::vector<LocalIndex> element_indices(mesh.elements().size());
    std::iota(element_indices.begin(), element_indices.end(), 0);

    // For a serial write, we do not need to queue up all of the partitions for writing, since
    // the file does not have to be preallocated for serial writing.
    auto partition =
      get_fea_io_partition(mesh, partition_id, node_indices, element_indices, mesh_info);

    // Overwrite the parent partitions
    for (auto& n_part : partition.node_parent_partitions)
    {
        n_part = 0;
    }
    for (auto& e_part : partition.element_parent_partitions)
    {
        e_part = 0;
    }

    mesh_writer.write_partition(partition);

    auto group_infos = get_group_infos(mesh);
    if (group_infos.size() > 0)
    {
        fea_io::GroupWriter group_writer;
        group_writer.initialize(file_stem, group_infos);
        for (auto& group_info : group_infos)
        {
            if (group_info.type == fea_io::GroupType::Node)
            {
                auto& nodes_in_group = mesh.nodes().group(group_info.name).node_indices;
                group_writer.write_list(group_info.name,
                                        0,
                                        static_cast<LocalIndex>(nodes_in_group.size()),
                                        nodes_in_group.data());
            }
            else if (group_info.type == fea_io::GroupType::Element)
            {
                auto& elements_in_group = mesh.elements().group(group_info.name).element_indices;
                group_writer.write_list(group_info.name,
                                        0,
                                        static_cast<LocalIndex>(elements_in_group.size()),
                                        elements_in_group.data());
            }
        }
        group_writer.close();
    }


    fea_io::DataWriter<GlobalIndex, LocalIndex, FloatType> data_writer;
    initialize_data_writer(data,
                           file_stem,
                           mesh_info,
                           get_fea_io_partition_info(mesh, mesh_info),
                           data_writer);

    int array_id = 0;
    for (auto& d : data.get<DataFieldType::NODE, double>())
    {
        data_writer.write(0, array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::ELEMENT, double>())
    {
        data_writer.write(0, array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::NODE_BY_COMPONENT, double>())
    {
        data_writer.write(0, array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::NODE_BY_ELEMENT, double>())
    {
        data_writer.write(0, array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::QUAD, double>())
    {
        data_writer.write(0, array_id++, d.data().data(), d.data().rows());
    }

    if (!data.empty())
    {
        data_writer.close();
    }
}

static void write_serial_multi_partition(const Mesh& mesh,
                                         const MeshAssociatedData& data,
                                         const std::string file_stem)
{
    fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer;
    auto mesh_info = get_fea_io_mesh_info(mesh);
    mesh_writer.initialize(file_stem, mesh_info);

    // Create an adjacency list to quickly get the neighbors of elements while finding the
    // layer of ghost elements.
    ElementAdjacencyList adjacency_list;
    adjacency_list.update(mesh);

    auto group_infos = get_group_infos(mesh);
    fea_io::GroupWriter group_writer;
    if (group_infos.size() > 0)
    {
        group_writer.initialize(file_stem, group_infos);
    }

    for (auto partition_id : utility::range(mesh.n_partitions()))
    {
        // Get the nodes and elements that should appear in the partition corresponding to
        // partition_id, including ghost elements.
        auto node_and_element_indices =
          extract_partition_indices(mesh, partition_id, adjacency_list);
        auto& node_indices    = node_and_element_indices.first;
        auto& element_indices = node_and_element_indices.second;

        // For a serial write, we do not need to queue up all of the partitions for writing, since
        // the file does not have to be preallocated for serial writing.
        auto partition =
          get_fea_io_partition(mesh, partition_id, node_indices, element_indices, mesh_info);
        mesh_writer.write_partition(partition);

        // Note: this could be done with more memory but less processing if a map was used instead
        // of set.  If it becomes a wall time bottleneck, reconsider.
        for (auto& group_info : group_infos)
        {
            // Note: nodes that are not in the partition get the max limit
            std::vector<LocalIndex> mesh_to_partition_node_map(
              mesh.nodes().size(),
              std::numeric_limits<LocalIndex>::max());
            for (auto i : utility::indices(node_indices))
            {
                mesh_to_partition_node_map[node_indices[i]] = static_cast<LocalIndex>(i);
            }
            std::vector<LocalIndex> mesh_to_partition_element_map(
              mesh.elements().size(),
              std::numeric_limits<LocalIndex>::max());
            for (auto i : utility::indices(element_indices))
            {
                mesh_to_partition_element_map[element_indices[i]] = static_cast<LocalIndex>(i);
            }

            if (group_info.type == fea_io::GroupType::Node)
            {
                auto& nodes_in_group = mesh.nodes().group(group_info.name).node_indices;
                std::vector<LocalIndex> node_indices;
                node_indices.reserve(node_indices.size());
                for (auto& index : nodes_in_group)
                {
                    if (mesh_to_partition_node_map[index] != std::numeric_limits<LocalIndex>::max())
                    {
                        node_indices.push_back(mesh_to_partition_node_map[index]);
                    }
                }
                group_writer.write_list(group_info.name,
                                        partition_id,
                                        static_cast<LocalIndex>(node_indices.size()),
                                        node_indices.data());
            }
            else if (group_info.type == fea_io::GroupType::Element)
            {
                auto& elements_in_group = mesh.elements().group(group_info.name).element_indices;
                std::vector<LocalIndex> element_indices;
                element_indices.reserve(elements_in_group.size());
                for (auto& index : elements_in_group)
                {
                    if (mesh_to_partition_element_map[index] !=
                        std::numeric_limits<LocalIndex>::max())
                    {
                        element_indices.push_back(mesh_to_partition_element_map[index]);
                    }
                }
                group_writer.write_list(group_info.name,
                                        partition_id,
                                        static_cast<LocalIndex>(element_indices.size()),
                                        element_indices.data());
            }
        }
    }

    if (group_infos.size() > 0)
    {
        group_writer.close();
    }

    if (!data.empty())
    {
        throw _UTIL_MSG_EXCEPT(
          "Writing data for centralized partitioned meshes is not implemented yet.");
    }
}

static void write_distributed_multi_partition(const Mesh& mesh,
                                              const MeshAssociatedData& data,
                                              const std::string file_stem)
{
    fea_io::MeshParWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer;
    auto mesh_info = get_fea_io_mesh_info(mesh);
    mesh_writer.initialize(file_stem, mesh_info);

    auto partition_info = get_fea_io_partition_info(mesh, mesh_info);
    mesh_writer.queue_partition_for_write(partition_info);

    mesh_writer.start_parallel_write_phase();

    // Get the nodes and elements that should appear in the partition corresponding to
    // partition_id, which is all of them.
    std::vector<LocalIndex> node_indices(mesh.nodes().size());
    std::iota(node_indices.begin(), node_indices.end(), 0);
    std::vector<LocalIndex> element_indices(mesh.nodes().size());
    std::iota(element_indices.begin(), element_indices.end(), 0);

    // For a serial write, we do not need to queue up all of the partitions for writing, since
    // the file does not have to be preallocated for serial writing.
    auto partition =
      get_fea_io_partition(mesh, mesh.partition_id(), node_indices, element_indices, mesh_info);

    mesh_writer.write_partition(partition);

    auto group_infos = get_group_infos(mesh);
    if (group_infos.size() > 0)
    {
        throw _UTIL_MSG_EXCEPT("Writing groups in parallel is not implemented yet.");
    }

    mesh_writer.close();


    fea_io::DataWriter<GlobalIndex, LocalIndex, FloatType> data_writer;
    initialize_data_writer(data, file_stem, mesh_info, partition_info, data_writer);

    int array_id = 0;
    for (auto& d : data.get<DataFieldType::NODE, double>())
    {
        data_writer.write(mesh.partition_id(), array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::ELEMENT, double>())
    {
        data_writer.write(mesh.partition_id(), array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::NODE_BY_COMPONENT, double>())
    {
        data_writer.write(mesh.partition_id(), array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::NODE_BY_ELEMENT, double>())
    {
        data_writer.write(mesh.partition_id(), array_id++, d.data().data(), d.data().rows());
    }
    for (auto& d : data.get<DataFieldType::QUAD, double>())
    {
        data_writer.write(mesh.partition_id(), array_id++, d.data().data(), d.data().rows());
    }

    if (!data.empty())
    {
        data_writer.close();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Header defines

Mesh read_fea_io_mesh(const std::string& file_stem, const PartitionIndex partition_id)
{
    fea_io::MeshReader<GlobalIndex, LocalIndex, double> mesh_reader(file_stem);
    auto mesh_info = mesh_reader.get_mesh_info();
    if (mesh_info.connectivity_scheme != fea_io::ConnectivityScheme::FlexMesh)
    {
        throw _UTIL_MSG_EXCEPT("The mesh was written with a different connectivity scheme.  A "
                               "conversion has not been written.");
    }

    Mesh mesh;
    auto partition = mesh_reader.read_partition(partition_id);

    if (mesh_reader.get_mesh_info().total_num_partitions > 1)
    {
        mesh.n_partitions() = mesh_reader.get_mesh_info().total_num_partitions;
        mesh.partition_id() = partition_id;
    }

    for (auto i : utility::indices(partition.node_global_ids))
    {
        const auto offset = mesh_info.dimensionality * i;
        Point p           = { 0, 0, 0 };
        switch (mesh_info.dimensionality)
        {
            case 3: p[2] = partition.node_coords[offset + 2];
            case 2: p[1] = partition.node_coords[offset + 1];
            case 1: p[0] = partition.node_coords[offset];
        }
        mesh.nodes().add(p, partition.node_global_ids[i], partition.node_parent_partitions[i]);
    }

    for (auto i : utility::indices(partition.element_global_ids))
    {
        auto& connectivity = partition.element_connectivities[i];

        mesh.elements().add(convert_element_type(partition.element_types[i]),
                            connectivity,
                            partition.element_global_ids[i],
                            partition.element_parent_partitions[i]);
    }

    for (auto& subcomponent : partition.subcomponents)
    {
        auto find_iter =
          std::find_if(std::begin(mesh_info.components),
                       std::end(mesh_info.components),
                       [id = subcomponent.component_id](const fea_io::ComponentInfo<LocalIndex>& c)
                       { return c.id == id; });

        if (find_iter == std::end(mesh_info.components))
        {
            throw _UTIL_MSG_EXCEPT("A subcomponent referenced a component id that did not exist in "
                                   "the mesh information.");
        }
        if (mesh_info.components[0].lcs_parameterization != fea_io::LCSParameterization::Quaternion)
        {
            throw _UTIL_MSG_EXCEPT(
              "Local coordinate systems were parameterized with something other than Quaternions. "
              "A conversion needs to be added to read_partition.");
        }

        // Create a group for the elements in the component
        auto& component_group = mesh.elements().group(find_iter->name);
        for (auto local_elem_id : subcomponent.elements_in_subcomponent)
        {
            mesh.elements()[local_elem_id].component_id() = subcomponent.component_id;
            component_group.element_indices.push_back(local_elem_id);
        }

        // Assign local coordinate systems for elements in the subcomponent
        LCSAssignmentType component_lcs_assignment =
          convert_lcs_assignment((*find_iter).lcs_assignment);
        std::size_t curr_lcs_index = 0;
        for (auto local_elem_id : subcomponent.elements_in_subcomponent)
        {
            auto element = mesh.elements()[local_elem_id];
            auto& lcs    = element.local_coord_sys();
            element.local_coord_sys().initialize(
              component_lcs_assignment,
              static_cast<ConnIndex>(element.connectivity().size()));
            for (auto i = 0; i < lcs.size(); i++)
            {
                auto& w = subcomponent.lcs_parameters_in_subcomponent[curr_lcs_index++];
                auto& x = subcomponent.lcs_parameters_in_subcomponent[curr_lcs_index++];
                auto& y = subcomponent.lcs_parameters_in_subcomponent[curr_lcs_index++];
                auto& z = subcomponent.lcs_parameters_in_subcomponent[curr_lcs_index++];
                auto q  = simple_math::Quaternion(w, x, y, z);
                lcs.set(i, q);
            }
        }
    }

    return mesh;
}

Mesh read_fea_io_mesh_from_meta(const std::string& meta_filepath, const PartitionIndex partition_id)
{
    fea_io::MetaReader meta_reader;
    meta_reader.open(meta_filepath);
    auto& analysis_meta_info = meta_reader.get_analysis_meta();
    auto mesh_relative_path =
      analysis_meta_info.input_info.levels_of_detail[0].mesh_objects[0].mesh_path;
    auto mesh_path = (fs::path(meta_filepath).parent_path() / mesh_relative_path).generic_string();

    Mesh mesh = read_fea_io_mesh(mesh_path, partition_id);

    for (auto group_relative_path :
         analysis_meta_info.input_info.levels_of_detail[0].mesh_objects[0].group_paths)
    {
        fea_io::GroupReader group_reader(
          (fs::path(meta_filepath).parent_path() / group_relative_path).string());

        auto node_group_names = group_reader.get_node_group_names();
        for (auto& name : node_group_names)
        {
            auto& node_group = mesh.nodes().group(name);

            std::vector<LocalIndex> indices;
            group_reader.read_list<LocalIndex>(name, partition_id, indices);
            node_group.node_indices = std::move(indices);
        }

        auto element_group_names = group_reader.get_element_group_names();
        for (auto& name : element_group_names)
        {
            auto& element_group = mesh.elements().group(name);

            std::vector<LocalIndex> indices;
            group_reader.read_list<LocalIndex>(name, partition_id, indices);
            element_group.element_indices = std::move(indices);
        }
    }

    return mesh;
}

void write_fea_io(const Mesh& mesh,
                  const std::string& base_directory,
                  const std::string& mesh_name,
                  const std::string meta_file_name)
{
    MeshAssociatedData data;
    write_fea_io(mesh, data, base_directory, mesh_name, meta_file_name);
}

void write_fea_io(const Mesh& mesh,
                  const MeshAssociatedData& data,
                  const std::string& base_directory,
                  const std::string& mesh_name,
                  const std::string meta_file_name)
{
    if (mesh.contains_all_partitions())
    { // Centralized mesh case
        bool has_groups = mesh.nodes().groups().size() > 0 || mesh.elements().groups().size();
        auto file_stem =
          write_meta_file(base_directory, mesh_name, meta_file_name, has_groups, !data.empty());
        if (mesh.n_partitions() == 1)
        {
            write_serial_single_partition(mesh, data, file_stem);
        }
        else
        {
            write_serial_multi_partition(mesh, data, file_stem);
        }
    }
    else
    { // Distributed mesh case
        // This method treats the mesh as a single partition.
        auto file_stem =
          write_meta_file(base_directory, mesh_name, meta_file_name, false, !data.empty());
        write_distributed_multi_partition(mesh, data, file_stem);
    }
}
}

#endif