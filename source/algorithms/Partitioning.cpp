#include "mesh/algorithms/Partitioning.h"

#include "mesh/relations/NodeToElementMap.h"
#include "second_party_metis/metis.h"

#include <set>

namespace mesh
{
void partition_metis(Mesh& mesh, const PartitionIndex n_parts)
{
    _UTIL_ASSERT(
      mesh.partition_id() == std::numeric_limits<PartitionIndex>::max(),
      "The mesh instance must contain the entire global mesh for the partition algorithm.");

    if (n_parts > 1)
    {
        // METIS input parameters
        idx_t n_elements = mesh.elements().size();
        idx_t n_nodes    = mesh.nodes().size();
        std::vector<idx_t> eptr(n_elements + 1);
        std::vector<idx_t> eind;
        for (auto& element : mesh.elements())
        {
            eptr[element.global_index()] = eind.size();
            for (auto& node : element.nodes())
            {
                eind.push_back(node.global_index());
            }
        }
        eptr.back()        = eind.size();
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


        // Set element partitions
        for (auto& element : mesh.elements())
        {
            element.parent_partition() =
              static_cast<PartitionIndex>(element_partitions[element.global_index()]);
        }

        // Loop over the nodes and ensure that the node is connected to at least one element with
        // the same partition. For some reason, this is not a requirement of the algorithm in METIS.
        NodeToElementMap node_to_element_map;
        node_to_element_map.update(mesh);
        for (auto n : mesh.nodes())
        {
            const auto& n_part_id = static_cast<PartitionIndex>(node_partitions[n.global_index()]);
            bool found_matching_part_id   = false;
            PartitionIndex last_e_part_id = n_part_id;
            for (auto e : node_to_element_map.get_elements(n.local_index()))
            {
                const auto& e_part_id =
                  static_cast<PartitionIndex>(element_partitions[e.global_index()]);
                found_matching_part_id |= n_part_id == e_part_id;
                last_e_part_id = e_part_id;
            }
            if (!found_matching_part_id) {
                n.parent_partition() = last_e_part_id;
            }
            else
            {
                n.parent_partition() = n_part_id;
            }
        }

        mesh.n_partitions() = n_parts;

        // Note: we should not change the global indices of the nodes or elements
    }
    else
    {
        for (auto& node : mesh.nodes())
        {
            node.parent_partition() = 0;
        }
        for (auto& element : mesh.elements())
        {
            element.parent_partition() = 0;
        }

        mesh.n_partitions() = n_parts;
    }
}

Mesh extract_partition_mesh(const Mesh& mesh, const PartitionIndex partition)
{
    Mesh mesh_part;
    mesh_part.n_partitions() = mesh.n_partitions();
    if (mesh_part.n_partitions() == 1)
    {
        mesh_part.partition_id() = mesh.partition_id();
    }
    else
    {
        mesh_part.partition_id() = partition;
    }

    ElementAdjacencyList adjacency;
    adjacency.update(mesh);

    // Determine the indices of elements belonging to the respective partition
    // or is adjacent to an element in the partition.
    std::set<LocalIndex> elements_to_keep;
    for (auto& element : mesh.elements())
    {
        if (element.parent_partition() == partition)
        {
            elements_to_keep.insert(element.local_index());
            for (auto& neighbor : adjacency.get_adjacent_elements(element.local_index()))
            {
                elements_to_keep.insert(neighbor.local_index());
            }
        }
    }

    mesh_part.copy(mesh, elements_to_keep);

    mesh_part.sort();

    return std::move(mesh_part);
}

std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> extract_partition_indices(
  const Mesh& mesh,
  const PartitionIndex partition,
  ElementAdjacencyList& adjacency_list,
  const bool include_ghost)
{
    // Loop over the elements
    std::set<LocalIndex> nodes_to_include_set;
    std::vector<LocalIndex> elements_to_include;
    for (auto& element : mesh.elements())
    {
        if (element.parent_partition() == partition)
        {
            elements_to_include.push_back(element.local_index());
            for (auto& node : element.nodes())
            {
                nodes_to_include_set.insert(node.local_index());
            }
        }
    }

    if (include_ghost)
    {
        adjacency_list.skip_in_neighbor_search(elements_to_include);

        std::set<LocalIndex> ghost_layer_elements;
        for (auto& e_local_id : elements_to_include)
        {
            for (auto& neighbor_index :
                 adjacency_list.get_adjacent_element_indices_with_skips(e_local_id))
            {
                ghost_layer_elements.insert(neighbor_index);
            }
        }

        for (auto& element_index : ghost_layer_elements)
        {
            elements_to_include.push_back(element_index);
            for (auto& node : mesh.elements()[element_index].nodes())
            {
                nodes_to_include_set.insert(node.local_index());
            }
        }
    }

    return { std::vector<LocalIndex>(nodes_to_include_set.begin(), nodes_to_include_set.end()),
             elements_to_include };
}
}
