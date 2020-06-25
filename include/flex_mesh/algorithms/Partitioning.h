#pragma once

#include "utility/Exceptions.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/Range.h"

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/mesh/MeshPartition.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/maps/AdjacencyLists.h"

#include "flex_mesh/algorithms/Search.h"

#include "second_party_metis/metis.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace flex_mesh
{
    /// Partitions the mesh using METIS.
    /// \param mesh - the mesh to partition
    /// \param n_parts - the number of partitions to break the mesh up into.
    template<typename MeshType>
    static void partition_metis(
        MeshType& mesh,
        const PartitionIndex n_parts);

    /// Extracts the mesh for a particular partition from the total mesh.
    /// \param mesh - the mesh to extract the partition mesh from (consumes this mesh)
    /// \param partition - the partition to extract.  If partition > # partitions in the mesh,
    ///	    then an empty mesh is returned.
    /// \important Consumes input parameter mesh.
    template<typename MeshType>
    static MeshPartition<typename MeshType::NodeTypename, typename  MeshType::ElementTypename>
    extract_partition_mesh(
        MeshType& mesh,
        const PartitionIndex partition);

    /// Returns the local indices of the nodes and elements respectively for the given partition.
    /// The lists of elements includes owned elements and a layer of ghost elements (if
    /// included_ghost is true).  The list of nodes includes all owned nodes and some nodes that
    /// are not owned but appear in the connectivity of owned (or ghost if included) elements.
    template<typename MeshType>
    static std::pair<std::vector<flex_mesh::LocalIndex>, std::vector<flex_mesh::LocalIndex>>
    extract_partition_indices(
        const MeshType& mesh,
        const PartitionIndex partition,
        const ElementAdjacencyList<MeshType>& adjacency_list,
        const bool include_ghost = true);

    //////////////////////////
    // Template definitions //
    //////////////////////////

    template<typename MeshType>
    void partition_metis(
        MeshType& mesh, 
        const PartitionIndex n_parts)
    {
        if (n_parts > 1)
        {
            // METIS input parameters
            idx_t n_elements = mesh.elements().size();
            idx_t n_nodes = mesh.nodes().size();
            std::vector<idx_t> eptr(n_elements + 1);
            std::vector<idx_t> eind;
            for (auto& e_ptr : mesh.elements())
            {
                eptr[e_ptr->get_global_number()] = eind.size();
                for (auto& n_ptr : *e_ptr)
                {
                    eind.push_back(n_ptr->get_global_number());
                }
            }
            eptr.back() = eind.size();
            idx_t n_partitions = n_parts;

            // METIS output parameters
            idx_t objective_function_value;
            std::vector<idx_t> element_partitions(n_elements);
            std::vector<idx_t> node_partitions(n_nodes);

            METIS_PartMeshNodal(&n_elements,
                &n_nodes,
                eptr.data(),
                eind.data(),
                nullptr, // vertex weights, treat as all equal
                nullptr, // vertex volumes, treat as all equal
                &n_partitions,
                nullptr, // target partition weights, treat all as equal
                nullptr, // options, use defaults
                &objective_function_value,
                element_partitions.data(),
                node_partitions.data());

            // Set node partitions
            for (auto& n_ptr : mesh.nodes())
            {
                n_ptr->set_parent_partition(
                    node_partitions[n_ptr->get_global_number()]);
            }
            // Set element partitions
            for (auto& e_ptr : mesh.elements())
            {
                e_ptr->set_parent_partition(
                    element_partitions[e_ptr->get_global_number()]);
            }

            mesh.nodes().sort_by_partition();
            mesh.elements().sort_by_partition();

            // Note: we should not change the global numbers of the nodes or elements
        }
        else
        {
            for (auto& n_ptr : mesh.nodes())
            {
                n_ptr->set_parent_partition(0);
            }
            for (auto& e_ptr : mesh.elements())
            {
                e_ptr->set_parent_partition(0);
            }
        }
    }

    template<typename MeshType>
    MeshPartition<typename MeshType::NodeTypename, typename MeshType::ElementTypename>
    extract_partition_mesh(MeshType& mesh, const PartitionIndex partition)
    {
        using NodeType = typename MeshType::NodeTypename;
        using ElementType = typename MeshType::ElementTypename;
        MeshPartition<NodeType, ElementType> part_mesh;
        part_mesh.move_to_self(mesh);

        // Determine the indices of elements belonging to the respective partition
        // or is adjacent to an element in the partition.
        std::set<flex_mesh::LocalIndex> elements_to_keep;
        flex_mesh::ElementAdjacencyList<MeshType> adjacency;
        adjacency.initialize(part_mesh);
        for (auto& e_ptr : part_mesh.elements())
        {
            if (e_ptr->get_parent_partition() == partition)
            {
                elements_to_keep.insert(e_ptr->get_local_number());
                for (auto neighbor_ptr : adjacency.get_neighbors(e_ptr->get_local_number()))
                {
                    elements_to_keep.insert(neighbor_ptr->get_local_number());
                }
            }
        }

        // Determine which nodes should be kept
        std::set<flex_mesh::LocalIndex> nodes_to_keep;
        for (auto& e_to_keep_index : elements_to_keep)
        {
            for (auto& n_ptr : *part_mesh.elements()[e_to_keep_index])
            {
                nodes_to_keep.insert(n_ptr->get_local_number());
            }
        }

        part_mesh.set_partition_id(partition);
        part_mesh.set_n_total_nodes(part_mesh.nodes().size());
        part_mesh.set_n_total_elements(part_mesh.elements().size());

        part_mesh.nodes().remove_all_except(nodes_to_keep);
        part_mesh.elements().remove_all_except(elements_to_keep);

        part_mesh.nodes().sort_by_partition();
        part_mesh.elements().sort_by_partition_owned_first(partition);

        return part_mesh;
    }

    template<typename MeshType>
    std::pair<std::vector<flex_mesh::LocalIndex>, std::vector<flex_mesh::LocalIndex>>
    extract_partition_indices(
        const MeshType& mesh,
        const PartitionIndex partition,
        ElementAdjacencyList<MeshType>& adjacency_list,
        const bool include_ghost = true)
    {
        // Loop over the elements
        std::set<flex_mesh::LocalIndex> nodes_to_include_set;
        std::vector<flex_mesh::LocalIndex> elements_to_include;
        for (auto& e_ptr : mesh.elements())
        {
            if (e_ptr->get_parent_partition() == partition)
            {
                elements_to_include.push_back(e_ptr->get_local_number());
                for (auto& n_ptr : *e_ptr)
                {
                    nodes_to_include_set.insert(n_ptr->get_local_number());
                }
            }
        }

        if (include_ghost)
        {
            adjacency_list.skip_in_neighbor_search(elements_to_include);

            std::set<flex_mesh::LocalIndex> ghost_layer_elements;
            for (auto& e_local_id : elements_to_include)
            {
                for (auto& neighbor_ptr : adjacency_list.get_neighbors(e_local_id))
                {
                    ghost_layer_elements.insert(neighbor_ptr->get_local_number());
                }
            }

            for (auto& e_local_id : ghost_layer_elements)
            {
                elements_to_include.push_back(e_local_id);
                for (auto& n_ptr : *mesh.elements()[e_local_id])
                {
                    nodes_to_include_set.insert(n_ptr->get_local_number());
                }
            }
        }

        return {
            std::vector<flex_mesh::LocalIndex>(nodes_to_include_set.begin(), nodes_to_include_set.end()),
            elements_to_include
        };
    }
}