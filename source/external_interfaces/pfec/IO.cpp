#include "mesh/external_interfaces/pfec/IO.h"

#include "mesh/Defines.h"
#include "mesh/relations/ElementToComponentMap.h"

//#include "pfec/PFECNativeLibInterface.h"
#include "lib/pfecDataReader.h"
#include "lib/pfecGroupReader.h"
#include "lib/pfecMeshReader.h"
#include "lib/pfecMetaReader.h"

namespace mesh
{
namespace detail
{
ElementType convert_element_type(const PFECNativeInterface::ElemGeomType& geomType)
{
    // 1D
    if (geomType == PFECNativeInterface::L1D2N)
        return ElementType::L1D2N;
    else if (geomType == PFECNativeInterface::L1D3N)
        return ElementType::L1D3N;
    // 2D
    else if (geomType == PFECNativeInterface::T2D3N)
        return ElementType::T2D3N;
    else if (geomType == PFECNativeInterface::Q2D4N)
        return ElementType::Q2D4N;
    else if (geomType == PFECNativeInterface::T2D6N)
        return ElementType::T2D6N;
    else if (geomType == PFECNativeInterface::Q2D8N)
        return ElementType::Q2D8N;
    else if (geomType == PFECNativeInterface::C2D4N)
        return ElementType::C2D4N;
    else if (geomType == PFECNativeInterface::C2D6N)
        return ElementType::C2D6N;
    // 3D
    else if (geomType == PFECNativeInterface::T3D4N)
        return ElementType::T3D4N;
    else if (geomType == PFECNativeInterface::P3D5N)
        return ElementType::P3D5N;
    else if (geomType == PFECNativeInterface::W3D6N)
        return ElementType::W3D6N;
    else if (geomType == PFECNativeInterface::H3D8N)
        return ElementType::H3D8N;
    else if (geomType == PFECNativeInterface::T3D10N)
        return ElementType::T3D10N;
    else if (geomType == PFECNativeInterface::P3D13N)
        return ElementType::P3D13N;
    else if (geomType == PFECNativeInterface::W3D15N)
        return ElementType::W3D15N;
    else if (geomType == PFECNativeInterface::H3D20N)
        return ElementType::H3D20N;
    else if (geomType == PFECNativeInterface::C3D8N)
        return ElementType::C3D8N;
    else if (geomType == PFECNativeInterface::C3D16N)
        return ElementType::C3D16N;
    else
        return ElementType::UNKNGEOM;
}

Mesh read_pfec_mesh_kernel(PFECNativeInterface::MetaReader& meta_reader,
                           const PartitionIndex partition,
                           const int time_step_id)
{
    Mesh mesh;

    // Read the mesh for the given time step
    auto mesh_reader = meta_reader.get_mesh_reader(time_step_id);
    {
        std::vector<PFECNativeInterface::PFEC_IND> global_numbers;
        std::vector<PFECNativeInterface::PFEC_INT> partitions;
        std::vector<std::vector<double>> coordinates;
        mesh_reader->read_part_nodal(partition, global_numbers, partitions, coordinates);
        for (auto i = 0; i < global_numbers.size(); i++)
        {
            assert(coordinates[i].size() == 3);
            Point p = { coordinates[i][0], coordinates[i][1], coordinates[i][2] };
            mesh.nodes().add(p,
                             static_cast<GlobalIndex>(global_numbers[i]),
                             static_cast<PartitionIndex>(partitions[i]));
        }
    }
    // TODO make an adapter for reading to avoid the overhead here
    {
        std::vector<PFECNativeInterface::PFEC_IND> global_numbers;
        std::vector<PFECNativeInterface::PFEC_INT> partitions;
        std::vector<PFECNativeInterface::PFEC_INT> domains;
        std::vector<PFECNativeInterface::ElemGeomType> element_types;
        std::vector<std::vector<PFECNativeInterface::PFEC_INT>> local_connectivities;
        mesh_reader->read_part_element(partition,
                                       global_numbers,
                                       partitions,
                                       domains,
                                       element_types,
                                       local_connectivities);
        std::vector<LocalIndex> connectivity;
        for (auto i = 0; i < global_numbers.size(); i++)
        {
            connectivity.resize(local_connectivities[i].size());
            for (auto j : utility::indices(local_connectivities[i]))
            {
                connectivity[j] = local_connectivities[i][j];
            }
            mesh.elements().add(detail::convert_element_type(element_types[i]),
                                connectivity,
                                static_cast<GlobalIndex>(global_numbers[i]),
                                static_cast<PartitionIndex>(partitions[i]),
                                static_cast<ComponentIndex>(domains[i]));
        }

        // Create a group for the domains
        std::set<PFECNativeInterface::PFEC_INT> domains_set;
        for (auto& d : domains)
        {
            domains_set.insert(d);
        }
        for (auto& domain : domains_set)
        {
            auto& group = mesh.elements().group("Component" + std::to_string(domain));
            for (auto i : utility::indices(domains))
            {
                if (domains[i] == domain)
                {
                    group.element_indices.push_back(static_cast<LocalIndex>(i));
                }
            }
        }
    }
    mesh_reader->close();

    // Read the orientation data
    {
        auto orient_reader = meta_reader.get_orient_reader(time_step_id);
        if (orient_reader)
        {
            std::vector<std::vector<std::array<double, 4>>> quaterions;
            orient_reader->read_part_orientation_data(partition, 0, quaterions);
            for (auto& element : mesh.elements())
            {
                auto local_index             = element.local_index();
                auto n_quaterions_in_element = quaterions[local_index].size();
                auto& lcs                    = element.local_coord_sys();
                if (n_quaterions_in_element > 0)
                {
                    if (n_quaterions_in_element == 1)
                    {
                        lcs.initialize(LCSAssignmentType::Element,
                                       static_cast<ConnIndex>(element.connectivity().size()));
                        auto& w = quaterions[local_index][0][0];
                        auto& x = quaterions[local_index][0][1];
                        auto& y = quaterions[local_index][0][2];
                        auto& z = quaterions[local_index][0][3];
                        auto q  = simple_math::Quaternion(w, x, y, z);
                        lcs.set(0, q.inverse());
                    }
                    else if (n_quaterions_in_element > 1)
                    {
                        // Node-by-element is not easy to convert to a unique nodal value for the
                        // subset of nodes in the component. So, for now, just average them to be
                        // by element.
                        lcs.initialize(LCSAssignmentType::Element,
                                       static_cast<ConnIndex>(element.connectivity().size()));
                        double w = 0.;
                        double x = 0.;
                        double y = 0.;
                        double z = 0.;
                        for (auto j = 0; j < n_quaterions_in_element; ++j)
                        {
                            w += quaterions[local_index][j][0];
                            x += quaterions[local_index][j][1];
                            y += quaterions[local_index][j][2];
                            z += quaterions[local_index][j][3];
                        }
                        auto q = simple_math::Quaternion(w, x, y, z);
                        q.normalize();
                        lcs.set(0, q.inverse());
                    }
                }
            }
            orient_reader->close();
        }
    }

    // Read the group information
    {
        auto group_paths = meta_reader.get_group_file_paths(time_step_id);

        for (std::string group_path : group_paths)
        {
            PFECNativeInterface::GroupReader group_reader;
            group_reader.open(group_path);
            std::vector<std::string> node_group_names;
            std::vector<std::string> element_group_names;
            group_reader.get_group_info(node_group_names, element_group_names);

            // Read node groups for this file
            for (std::string node_group_name : node_group_names)
            {
                auto& group = mesh.nodes().group(node_group_name);

                std::vector<PFECNativeInterface::PFEC_IND> indices;
                bool uses_local_numbering;
                group_reader.get_part_node_group(partition,
                                                 node_group_name,
                                                 indices,
                                                 uses_local_numbering);
                if (uses_local_numbering)
                {
                    for (auto index : indices)
                    {
                        group.node_indices.push_back(static_cast<LocalIndex>(index));
                    }
                }
                else
                {
                    throw _UTIL_MSG_EXCEPT("Reading group info that uses global numbering has not "
                                           "been implemented yet.");
                }
            }

            // Read element groups for this file
            for (std::string element_group_name : element_group_names)
            {
                auto& group = mesh.elements().group(element_group_name);

                std::vector<PFECNativeInterface::PFEC_IND> indices;
                bool uses_local_numbering;
                group_reader.get_part_elem_group(partition,
                                                 element_group_name,
                                                 indices,
                                                 uses_local_numbering);
                if (uses_local_numbering)
                {
                    for (auto index : indices)
                    {
                        group.element_indices.push_back(static_cast<LocalIndex>(index));
                    }
                }
                else
                {
                    throw _UTIL_MSG_EXCEPT("Reading group info that uses global numbering has not "
                                           "been implemented yet.");
                }
            }

            group_reader.close();
        }
    }

    return mesh;
}


MeshAssociatedData read_pfec_data_kernel(PFECNativeInterface::DataReader& data_reader,
                                         const PartitionIndex partition)
{
    MeshAssociatedData data;
    auto data_infos = data_reader.get_data_info();
    for (auto& data_info : data_infos)
    {
        simple_math::GeneralMatrix extracted_data;

        if (data_info.type == PFECNativeInterface::DataType::Nodal)
        {
            std::vector<std::vector<double>> pfec_data;
            data_reader.read_part_unique_data(partition, data_info.id, pfec_data);

            extracted_data =
              simple_math::GeneralMatrix::Zero(pfec_data.size(), data_info.num_comps);
            for (auto i : utility::indices(pfec_data))
            {
                assert(pfec_data[i].size() == extracted_data.cols());
                for (auto j : utility::indices(pfec_data[i]))
                {
                    extracted_data(i, j) = pfec_data[i][j];
                }
            }

            data.add<DataFieldType::NODE, double>(
              DataField<double>::copy(data_info.name,
                                      data_info.comp_names,
                                      std::move(extracted_data)));
        }
        else if (data_info.type == PFECNativeInterface::DataType::Element)
        {
            std::vector<std::vector<double>> pfec_data;
            data_reader.read_part_unique_data(partition, data_info.id, pfec_data);

            extracted_data =
              simple_math::GeneralMatrix::Zero(pfec_data.size(), data_info.num_comps);
            for (auto i : utility::indices(pfec_data))
            {
                assert(pfec_data[i].size() == extracted_data.cols());
                for (auto j : utility::indices(pfec_data[i]))
                {
                    extracted_data(i, j) = pfec_data[i][j];
                }
            }

            data.add<DataFieldType::ELEMENT, double>(
              DataField<double>::copy(data_info.name,
                                      data_info.comp_names,
                                      std::move(extracted_data)));
        }
        else if (data_info.type == PFECNativeInterface::DataType::NodalByElement)
        {
            std::vector<std::vector<std::vector<double>>> pfec_data;
            data_reader.read_part_by_elem_data(partition, data_info.id, pfec_data);

            if (data.node_by_element_offsets.size() == 0)
            {
                data.node_by_element_offsets.resize(pfec_data.size());
                data.node_by_element_offsets[0] = 0;
                for (auto i : utility::indices(pfec_data))
                {
                    data.node_by_element_offsets[i + 1] =
                      data.node_by_element_offsets[i] + pfec_data[i].size();
                }
            }

            extracted_data =
              simple_math::GeneralMatrix::Zero(data.node_by_element_offsets.back() - 1,
                                               data_info.num_comps);
            std::size_t curr_row = 0;
            for (auto i : utility::indices(pfec_data))
            {
                for (auto j : utility::indices(pfec_data[i]))
                {
                    for (auto k : utility::indices(pfec_data[i][j]))
                    {
                        extracted_data(curr_row, j) = pfec_data[i][j][k];
                    }
                    curr_row++;
                }
            }

            data.add<DataFieldType::NODE_BY_ELEMENT, double>(
              DataField<double>::copy(data_info.name,
                                      data_info.comp_names,
                                      std::move(extracted_data)));
        }
        else if (data_info.type == PFECNativeInterface::DataType::QuadPoint)
        {
            std::vector<std::vector<std::vector<double>>> pfec_data;
            data_reader.read_part_by_elem_data(partition, data_info.id, pfec_data);

            if (data.quad_point_offsets.size() == 0)
            {
                data.quad_point_offsets.resize(pfec_data.size());
                data.quad_point_offsets[0] = 0;
                for (auto i : utility::indices(pfec_data))
                {
                    data.quad_point_offsets[i + 1] =
                      data.quad_point_offsets[i] + pfec_data[i].size();
                }
            }

            extracted_data = simple_math::GeneralMatrix::Zero(data.quad_point_offsets.back() - 1,
                                                              data_info.num_comps);
            std::size_t curr_row = 0;
            for (auto i : utility::indices(pfec_data))
            {
                for (auto j : utility::indices(pfec_data[i]))
                {
                    for (auto k : utility::indices(pfec_data[i][j]))
                    {
                        extracted_data(curr_row, j) = pfec_data[i][j][k];
                    }
                    curr_row++;
                }
            }

            data.add<DataFieldType::QUAD, double>(
              DataField<double>::copy(data_info.name,
                                      data_info.comp_names,
                                      std::move(extracted_data)));
        }
    }

    return data;
}

std::pair<Mesh, MeshAssociatedData> read_pfec_all_parts_kernel(
  PFECNativeInterface::MetaReader& meta_reader,
  const int time_step_id,
  const bool convert_node_data_to_element_data)
{
    Mesh mesh;

    // Read the mesh for the given time step
    auto mesh_reader = meta_reader.get_mesh_reader(time_step_id);
    auto partitions  = mesh_reader->get_parts_in_file();

    auto n_total_nodes    = mesh_reader->get_num_global_nodes();
    auto n_total_elements = mesh_reader->get_num_global_elems();
    auto orient_reader    = meta_reader.get_orient_reader(time_step_id);
    auto group_paths      = meta_reader.get_group_file_paths(time_step_id);
    std::vector<PFECNativeInterface::GroupReader> group_readers;
    for (std::string group_path : group_paths)
    {
        group_readers.push_back(PFECNativeInterface::GroupReader());
        group_readers.back().open(group_path);
    }

    // Node data
    std::vector<PFECNativeInterface::PFEC_IND> n_global_numbers(n_total_nodes);
    std::iota(n_global_numbers.begin(), n_global_numbers.end(), 0);
    std::vector<PFECNativeInterface::PFEC_INT> n_partitions(n_total_nodes);
    std::vector<std::vector<double>> n_coordinates(n_total_nodes);

    // Element data
    std::vector<PFECNativeInterface::PFEC_IND> e_global_numbers(n_total_elements);
    std::iota(e_global_numbers.begin(), e_global_numbers.end(), 0);
    std::vector<PFECNativeInterface::PFEC_INT> e_partitions(n_total_elements);
    std::vector<PFECNativeInterface::PFEC_INT> e_domains(n_total_elements);
    std::vector<PFECNativeInterface::ElemGeomType> e_element_types(n_total_elements);
    std::vector<std::vector<PFECNativeInterface::PFEC_INT>> e_connectivities(n_total_elements);
    std::vector<std::vector<std::array<double, 4>>> e_quaterions(n_total_elements);

    std::vector<std::vector<PFECNativeInterface::PFEC_IND>> parts_n_local_to_global_map;
    parts_n_local_to_global_map.reserve(partitions.size());
    std::vector<std::vector<PFECNativeInterface::PFEC_IND>> parts_e_local_to_global_map;
    parts_e_local_to_global_map.reserve(partitions.size());
    std::vector<std::vector<std::vector<PFECNativeInterface::PFEC_INT>>>
      parts_e_local_connectivities;
    parts_e_local_connectivities.reserve(partitions.size());

    for (auto& partition : partitions)
    {
        std::vector<PFECNativeInterface::PFEC_IND> part_n_global_numbers;
        std::vector<PFECNativeInterface::PFEC_INT> part_n_partitions;
        std::vector<std::vector<double>> part_n_coordinates;

        mesh_reader->read_part_nodal(partition,
                                     part_n_global_numbers,
                                     part_n_partitions,
                                     part_n_coordinates);

        parts_n_local_to_global_map.push_back(part_n_global_numbers);

        for (auto i = 0; i < part_n_global_numbers.size(); i++)
        {
            assert(part_n_coordinates[i].size() == 3);

            auto& global_id          = part_n_global_numbers[i];
            n_partitions[global_id]  = part_n_partitions[i];
            n_coordinates[global_id] = part_n_coordinates[i];
        }

        std::vector<PFECNativeInterface::PFEC_IND> part_e_global_numbers;
        std::vector<PFECNativeInterface::PFEC_INT> part_e_partitions;
        std::vector<PFECNativeInterface::PFEC_INT> part_e_domains;
        std::vector<PFECNativeInterface::ElemGeomType> part_e_element_types;
        std::vector<std::vector<PFECNativeInterface::PFEC_INT>> part_e_local_connectivities;

        mesh_reader->read_part_element(partition,
                                       part_e_global_numbers,
                                       part_e_partitions,
                                       part_e_domains,
                                       part_e_element_types,
                                       part_e_local_connectivities);

        parts_e_local_to_global_map.push_back(part_e_global_numbers);
        parts_e_local_connectivities.push_back(part_e_local_connectivities);

        for (auto i = 0; i < part_e_global_numbers.size(); i++)
        {
            auto& global_id            = part_e_global_numbers[i];
            e_partitions[global_id]    = part_e_partitions[i];
            e_domains[global_id]       = part_e_domains[i];
            e_element_types[global_id] = part_e_element_types[i];
            // Remap the connectivity to global node indices
            for (auto& c : part_e_local_connectivities[i])
            {
                c = part_n_global_numbers[c];
            }
            e_connectivities[global_id] = part_e_local_connectivities[i];
        }

        if (orient_reader)
        {
            std::vector<std::vector<std::array<double, 4>>> part_e_quaterions;
            orient_reader->read_part_orientation_data(partition, 0, part_e_quaterions);
            for (auto i = 0; i < part_e_global_numbers.size(); i++)
            {
                auto& global_id         = part_e_global_numbers[i];
                e_quaterions[global_id] = part_e_quaterions[i];
            }
        }

        for (auto& group_reader : group_readers)
        {
            std::vector<std::string> node_group_names;
            std::vector<std::string> element_group_names;
            group_reader.get_group_info(node_group_names, element_group_names);

            // Read node groups for this file
            for (auto& node_group_name : node_group_names)
            {
                auto& group = mesh.nodes().group(node_group_name);

                std::vector<PFECNativeInterface::PFEC_IND> local_indices;
                bool uses_local_numbering;
                group_reader.get_part_node_group(partition,
                                                 node_group_name,
                                                 local_indices,
                                                 uses_local_numbering);
                if (uses_local_numbering)
                {
                    for (auto local_index : local_indices)
                    {
                        group.node_indices.push_back(part_n_global_numbers[local_index]);
                    }
                }
                else
                {
                    throw _UTIL_MSG_EXCEPT("Reading group info that uses global numbering has not "
                                           "been implemented yet.");
                }
            }

            // Read element groups for this file
            for (std::string element_group_name : element_group_names)
            {
                auto& group = mesh.elements().group(element_group_name);

                std::vector<PFECNativeInterface::PFEC_IND> local_indices;
                bool uses_local_numbering;
                group_reader.get_part_elem_group(partition,
                                                 element_group_name,
                                                 local_indices,
                                                 uses_local_numbering);
                if (uses_local_numbering)
                {
                    for (auto local_index : local_indices)
                    {
                        group.element_indices.push_back(part_e_global_numbers[local_index]);
                    }
                }
                else
                {
                    throw _UTIL_MSG_EXCEPT("Reading group info that uses global numbering has not "
                                           "been implemented yet.");
                }
            }
        }
    }

    // Add nodes to the mesh
    // Note: an index map must be created because the mesh may only be a subset (clipped down), and
    // thus not every global index will be in this mesh.
    std::vector<GlobalIndex> n_global_to_new_local(n_total_nodes, NOT_SET_GlobalIndex);
    GlobalIndex next_new_local_index = 0;
    for (auto i : utility::range(n_total_nodes))
    {
        if (!n_coordinates[i].empty())
        {
            n_global_to_new_local[i] = next_new_local_index++;
            Point p = { n_coordinates[i][0], n_coordinates[i][1], n_coordinates[i][2] };
            mesh.nodes().add(p,
                             static_cast<GlobalIndex>(n_global_numbers[i]),
                             static_cast<PartitionIndex>(n_partitions[i]));
        }
    }
    // Permute the node indices in the groups
    for (auto& group : mesh.nodes().groups())
    {
        for (auto& n : group.node_indices)
        {
            n = n_global_to_new_local[n];
        }
    }

    // Add elements to the mesh
    std::vector<GlobalIndex> e_global_to_new_local(n_total_elements, NOT_SET_GlobalIndex);
    next_new_local_index = 0;
    std::vector<LocalIndex> connectivity;
    for (auto i : utility::range(n_total_elements))
    {
        if (!e_connectivities[i].empty())
        {
            e_global_to_new_local[i] = next_new_local_index++;

            connectivity.resize(e_connectivities[i].size());
            for (auto j : utility::indices(e_connectivities[i]))
            {
                connectivity[j] = n_global_to_new_local[e_connectivities[i][j]];
            }
            mesh.elements().add(detail::convert_element_type(e_element_types[i]),
                                connectivity,
                                static_cast<GlobalIndex>(e_global_numbers[i]),
                                static_cast<PartitionIndex>(e_partitions[i]),
                                static_cast<ComponentIndex>(e_domains[i]));

            auto n_quaterions_in_element = e_quaterions[i].size();
            if (n_quaterions_in_element > 0)
            {
                auto& lcs = mesh.elements().back().local_coord_sys();

                if (n_quaterions_in_element == 1)
                {
                    lcs.initialize(LCSAssignmentType::Element,
                                   static_cast<ConnIndex>(connectivity.size()));
                }
                else if (n_quaterions_in_element > 1)
                {
                    lcs.initialize(LCSAssignmentType::Node,
                                   static_cast<ConnIndex>(connectivity.size()));
                }

                for (auto j = 0; j < n_quaterions_in_element; ++j)
                {
                    auto& w = e_quaterions[i][j][0];
                    auto& x = e_quaterions[i][j][1];
                    auto& y = e_quaterions[i][j][2];
                    auto& z = e_quaterions[i][j][3];
                    auto q  = simple_math::Quaternion(w, x, y, z);
                    lcs.set(j, q.inverse());
                }
            }
        }
    }
    // Permute the node indices in the groups
    for (auto& group : mesh.elements().groups())
    {
        for (auto& n : group.element_indices)
        {
            n = e_global_to_new_local[n];
        }
    }


    // Create a group for the domains
    for (auto& domain_id : mesh_reader->get_domain_ids())
    {
        mesh.elements().group("Component" + std::to_string(domain_id));
    }
    std::map<int, ElementGroup*> domain_to_group;
    for (auto& domain_id : mesh_reader->get_domain_ids())
    {
        domain_to_group[domain_id] =
          &mesh.elements().group("Component" + std::to_string(domain_id));
    }
    for (auto e : mesh.elements())
    {
        domain_to_group[e.component_id()]->element_indices.push_back(e.local_index());
    }

    mesh_reader->close();
    if (orient_reader)
    {
        orient_reader->close();
    }
    for (auto& group_reader : group_readers)
    {
        group_reader.close();
    }


    // Note: to read NodalByElement data, it is complicated.  For now, just average the
    // NodalByElement data to create equivalent element data (which is coarser of course).
    MeshAssociatedData data;
    auto data_reader = meta_reader.get_data_reader(time_step_id);
    if (data_reader)
    {
        auto data_infos = data_reader->get_data_info();
        for (auto& data_info : data_infos)
        {
            if (data_info.type == PFECNativeInterface::DataType::Element)
            {
                data.add<DataFieldType::ELEMENT, double>(DataField<double>::copy(
                  data_info.name,
                  data_info.comp_names,
                  simple_math::GeneralMatrix::Zero(mesh.elements().size(), data_info.num_comps)));
            }
        }
        for (auto& data_info : data_infos)
        {
            if (data_info.type == PFECNativeInterface::DataType::Nodal)
            {
                if (!convert_node_data_to_element_data)
                {
                    data.add<DataFieldType::NODE, double>(DataField<double>::copy(
                      data_info.name,
                      data_info.comp_names,
                      simple_math::GeneralMatrix::Zero(mesh.nodes().size(), data_info.num_comps)));
                }
                else
                {
                    data.add<DataFieldType::ELEMENT, double>(DataField<double>::copy(
                      data_info.name,
                      data_info.comp_names,
                      simple_math::GeneralMatrix::Zero(mesh.elements().size(),
                                                       data_info.num_comps)));
                }
            }
        }
        for (auto& data_info : data_infos)
        {
            if (data_info.type == PFECNativeInterface::DataType::NodalByElement)
            {
                data.add<DataFieldType::ELEMENT, double>(DataField<double>::copy(
                  data_info.name,
                  data_info.comp_names,
                  simple_math::GeneralMatrix::Zero(mesh.elements().size(), data_info.num_comps)));
            }
            else if (data_info.type == PFECNativeInterface::DataType::QuadPoint)
            {
                throw _UTIL_MSG_EXCEPT("QuadPoint not implemented when reading all parts.");
            }
        }

        std::set<LocalIndex> debug_visited_elems;
        for (auto i : utility::indices(partitions))
        {
            auto& partition             = partitions[i];
            auto& n_local_to_global_map = parts_n_local_to_global_map[i];
            auto& e_local_to_global_map = parts_e_local_to_global_map[i];
            auto part_data              = read_pfec_data_kernel(*data_reader, partition);

            const std::vector<DataField<double>>& part_element_datas =
              part_data.get<DataFieldType::ELEMENT, double>();
            auto& total_element_datas = data.get<DataFieldType::ELEMENT, double>();
            for (auto j : utility::indices(part_element_datas))
            {
                auto& part_element_data = part_element_datas[j].data();
                auto total_element_data = total_element_datas[j].data_mutable();

                for (auto part_e_local_id : utility::indices(e_local_to_global_map))
                {
                    auto& new_local_id =
                      e_global_to_new_local[e_local_to_global_map[part_e_local_id]];

                    total_element_data.row(new_local_id) = part_element_data.row(part_e_local_id);
                    debug_visited_elems.insert(new_local_id);
                }
            }

            const std::vector<DataField<double>>& part_node_datas =
              part_data.get<DataFieldType::NODE, double>();
            auto& total_node_datas = data.get<DataFieldType::NODE, double>();
            for (auto j : utility::indices(part_node_datas))
            {
                auto& part_node_data = part_node_datas[j].data();

                if (convert_node_data_to_element_data)
                {
                    auto total_element_data =
                      total_element_datas[j + part_element_datas.size()].data_mutable();
                    auto& part_e_local_connectivities = parts_e_local_connectivities[i];

                    for (auto part_e_local_id : utility::indices(e_local_to_global_map))
                    {
                        auto& new_local_id =
                          e_global_to_new_local[e_local_to_global_map[part_e_local_id]];
                        auto& local_conn = part_e_local_connectivities[part_e_local_id];

                        total_element_data.row(new_local_id).setZero();
                        for (auto& c : local_conn)
                        {

                            total_element_data.row(new_local_id) += part_node_data.row(c);
                        }
                        total_element_data.row(new_local_id) /= local_conn.size();
                    }
                }
                else
                {
                    auto total_node_data = total_node_datas[j].data_mutable();

                    for (auto part_n_local_id : utility::indices(n_local_to_global_map))
                    {
                        auto& new_local_id =
                          n_global_to_new_local[n_local_to_global_map[part_n_local_id]];

                        total_node_data.row(new_local_id) = part_node_data.row(part_n_local_id);
                    }
                }
            }

            int offset = 0;
            if (convert_node_data_to_element_data)
            {
                offset += part_node_datas.size() + part_element_datas.size();
            }
            else
            {
                offset += part_element_datas.size();
            }

            const std::vector<DataField<double>>& part_node_by_element_datas =
              part_data.get<DataFieldType::NODE_BY_ELEMENT, double>();
            for (auto j : utility::indices(part_node_by_element_datas))
            {
                auto& part_node_by_element_data = part_node_by_element_datas[j].data();
                auto& part_offsets              = part_data.node_by_element_offsets;
                auto total_element_data         = total_element_datas[j + offset].data_mutable();

                for (auto part_e_local_id : utility::indices(e_local_to_global_map))
                {
                    auto& new_local_id =
                      e_global_to_new_local[e_local_to_global_map[part_e_local_id]];

                    total_element_data.row(new_local_id).setZero();
                    for (auto k : utility::range(part_offsets[part_e_local_id],
                                                 part_offsets[part_e_local_id + 1]))
                    {

                        total_element_data.row(new_local_id) += part_node_by_element_data.row(k);
                    }
                    total_element_data.row(new_local_id) /=
                      part_offsets[part_e_local_id + 1] - part_offsets[part_e_local_id];
                }
            }
        }

        data_reader->close();
    }

    return { mesh, std::move(data) };
}

} // end namespace detail

std::pair<Mesh, MeshAssociatedData> read_pfec_all_parts(
  const std::string& meta_file_path,
  const int time_step_id,
  const bool convert_node_data_to_element_data)
{
    PFECNativeInterface::MetaReader meta_reader;
    meta_reader.read_meta_file(meta_file_path);

    return detail::read_pfec_all_parts_kernel(meta_reader,
                                              time_step_id,
                                              convert_node_data_to_element_data);
}

std::pair<Mesh, MeshAssociatedData> read_pfec(const std::string& meta_file_path,
                                              const PartitionIndex partition,
                                              const int time_step_id)
{
    return { read_pfec_mesh(meta_file_path, partition, time_step_id),
             read_pfec_data(meta_file_path, partition, time_step_id) };
}

Mesh read_pfec_mesh(const std::string& meta_file_path,
                    const PartitionIndex partition,
                    const int time_step_id)
{
    PFECNativeInterface::MetaReader meta_reader;
    meta_reader.read_meta_file(meta_file_path);

    return detail::read_pfec_mesh_kernel(meta_reader, partition, time_step_id);
}

MeshAssociatedData read_pfec_data(const std::string& meta_file_path,
                                  const PartitionIndex partition,
                                  const int time_step_id)
{
    PFECNativeInterface::MetaReader meta_reader;
    meta_reader.read_meta_file(meta_file_path);
    auto data_reader = meta_reader.get_data_reader(time_step_id);
    return detail::read_pfec_data_kernel(*data_reader, partition);
}



void convert_pfec_to_fea_io(const std::string& pfec_meta_file_path,
                            const std::string& fea_io_path_stem)
{
    PFECNativeInterface::MetaReader meta_reader;
    meta_reader.read_meta_file(pfec_meta_file_path);

    auto n_time_steps = meta_reader.get_num_time_steps();

    // Read the mesh for the given time step
    auto mesh_reader = meta_reader.get_mesh_reader(0);

    throw _UTIL_MSG_EXCEPT("TODO: Implementation not complete.");
}
}