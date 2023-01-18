#pragma once

#include "mesh/Mesh.h"
#include "mesh/relations/ElementAdjacencyList.h"

namespace mesh
{

/// @brief Partitions the mesh using METIS.
/// @param mesh - the mesh to partition
/// @param n_parts - the number of partitions to break the mesh up into.
void partition_metis(Mesh& mesh, const PartitionIndex n_parts);

/// @brief Extracts the mesh for a particular partition from the total mesh.
/// @param mesh - the mesh to extract the partition mesh from (consumes this mesh)
/// @param partition - the partition to extract.  If partition > # partitions in the mesh,
///	    then an empty mesh is returned.
/// @attention Consumes input parameter mesh.
Mesh extract_partition_mesh(const Mesh& mesh, const PartitionIndex partition);

/// @brief Returns the local indices of the nodes and elements respectively for the given partition.
/// The lists of elements includes owned elements and a layer of ghost elements (if
/// included_ghost is true).  The list of nodes includes all owned nodes and some nodes that
/// are not owned but appear in the connectivity of owned (or ghost if included) elements.
/// @param mesh - mesh to extract from
/// @param partition - 0-based index of the target partition
/// @param adjacency_list - note: changes skip flags if include_ghost is set to true.
/// @param include_ghost - indicates whether to include indices of ghost nodes/elements (owned by another partition)
std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> extract_partition_indices(
  const Mesh& mesh,
  const PartitionIndex partition,
  ElementAdjacencyList& adjacency_list,
  const bool include_ghost = true);
}