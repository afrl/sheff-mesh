#pragma once

#include "ElementAccessor.h"
#include "mesh/Defines.h"

namespace mesh
{
/// Given an element and a hash of the connectivity indices for which the face is being searched,
/// find the index of the face in the element that contains all of the connectivity entries.
/// @returns A pair indicating if it was successful and the index of the face relative to the
/// element.  Note: if the bool is false, then the index is meaningless.
std::pair<bool, ConnIndex> find_face_index(const ElementConstAccessor element,
                                           const std::size_t conn_hash);

/// Given an element and a hash of the connectivity indices for which the edge is being searched,
/// find the index of the edge in the element that contains all of the connectivity entries.
/// @returns A pair indicating if it was successful and the index of the edge relative to the
/// element.  Note: if the bool is false, then the index is meaningless.
std::pair<bool, ConnIndex> find_edge_index(const ElementConstAccessor element,
                                           const std::size_t conn_hash);

/// Given an element and a list of connectivity indices, this finds the element topology component
/// that contains all of the connectivity entries (order insensitive).
///
/// This is particularly useful if you have some list of connectivity entries and you want to find
/// if a face, edge, or just a node match the connectivity.
ElementTopologyComponent find_topology_component(const ElementConstAccessor element,
                                                 const std::vector<ConnIndex>& connectivity);

/// Very similiar to find_topology_component but based on shared nodes.
///
/// @param element - the element to search for the topology component
/// @param shared_connectivities - the list of shared connectivity entries between two elements. One
///        of which should be the element passed as a parameter.
/// @param in_first_element - bool indicating if the given element is the first in the list of
///        shared connectivities or the second.
ElementTopologyComponent find_shared_topology_component(
  const ElementConstAccessor element,
  const bool in_first_element,
  const std::vector<std::pair<ConnIndex, ConnIndex>>& shared_connectivities);

std::vector<LocalIndex> find_shared_nodes(
  const ElementConstAccessor element_1,
  const ElementConstAccessor element_2);

std::vector<std::pair<ConnIndex, ConnIndex>> find_shared_node_positions(
  const ElementConstAccessor element_1,
  const ElementConstAccessor element_2);

/// Returns a mean coordinate of the nodes which should be interior to the element in most cases.
Point calc_mean_coordinate(const ElementConstAccessor element);

}