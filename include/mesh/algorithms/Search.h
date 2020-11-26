#pragma once

#include "mesh/Mesh.h"

#include <utility>

namespace mesh
{
/// @brief Finds extent of the mesh (most negative and positive point forming the bounding box of
/// the mesh).
/// @returns A pair of points with most negative point being first and most positive being second.
std::pair<Point, Point> find_extent(const Mesh& mesh);

/// @brief Returns the extent of multiple extents.
std::pair<Point, Point> find_extent(const std::vector<std::pair<Point, Point>>& extents);
}