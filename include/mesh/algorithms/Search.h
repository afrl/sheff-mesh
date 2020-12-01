#pragma once

#include "mesh/Mesh.h"
#include "mesh/surface/Surface.h"

#include <utility>

namespace mesh
{
/// @brief Finds extent of the mesh (most negative and positive point forming the bounding box of
/// the mesh).
/// @returns A pair of points with most negative point being first and most positive being second.
std::pair<Point, Point> find_extent(const Mesh& mesh);

/// @brief Returns the extent of multiple extents.
std::pair<Point, Point> find_extent(const std::vector<std::pair<Point, Point>>& extents);

/// @brief TODO write doc string
/// @param mesh 
/// @param elements_to_search 
MeshSurface extract_boundary(Mesh& mesh, const std::vector<LocalIndex>& elements_to_search);

/// A helper struct for find_faces_from_node_set
struct ElementFacePair
{
    LocalIndex local_element_index;
    ConnIndex face_index;
};

/// Find all faces such that all nodes in the face appears in the given node set.
std::vector<ElementFacePair> find_faces_from_node_set(
  const Mesh& mesh,
  const std::set<LocalIndex>& node_set);

}