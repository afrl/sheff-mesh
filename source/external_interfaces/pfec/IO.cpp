#include "mesh/external_interfaces/pfec/IO.h"

//#include "pfec/PFECNativeLibInterface.h"
#include "lib/pfecDataReader.h"
#include "lib/pfecGroupReader.h"
#include "lib/pfecMeshReader.h"
#include "lib/pfecMetaReader.h"

namespace mesh
{
static ElementType convert_element_type(const PFECNativeInterface::ElemGeomType& geomType)
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

Mesh read_pfec(const std::string& meta_file_path, const PartitionIndex partition)
{
    Mesh mesh;

    PFECNativeInterface::MetaReader meta_reader;
    meta_reader.read_meta_file(meta_file_path);

    // Read the mesh for the given time step
    auto mesh_reader = meta_reader.get_mesh_reader(0);
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
            mesh.elements().add(convert_element_type(element_types[i]),
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
            auto& group = mesh.elements().group("Domain" + std::to_string(domain));
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
        auto orient_reader = meta_reader.get_orient_reader(0);
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
                    }
                    else if (n_quaterions_in_element > 1)
                    {
                        lcs.initialize(LCSAssignmentType::Node,
                                       static_cast<ConnIndex>(element.connectivity().size()));
                    }

                    for (auto j = 0; j < n_quaterions_in_element; ++j)
                    {
                        auto& w = quaterions[local_index][j][0];
                        auto& x = quaterions[local_index][j][1];
                        auto& y = quaterions[local_index][j][2];
                        auto& z = quaterions[local_index][j][3];
                        auto q  = simple_math::Quaternion(w, x, y, z);
                        lcs.set(j, q.inverse());
                    }
                }
            }
            orient_reader->close();
        }
    }

    // Read the group information
    {
        auto group_paths = meta_reader.get_group_file_paths(0);

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
}