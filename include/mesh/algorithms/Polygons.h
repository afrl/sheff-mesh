#pragma once

#include "mesh/Mesh.h"
#include "simple_math/EigenDefines.h"
#include "utility/DynamicBitset.h"

#include <unordered_map>

namespace mesh
{

struct SegmentStats
{
    double mean;
    double std_dev;
};

std::vector<std::pair<ConnIndex, ConnIndex>> polygon_segments(
  const utility::Slice<const LocalIndex, ConnIndex>& connectivity);

std::vector<std::pair<ConnIndex, ConnIndex>> polygon_segments(
  const utility::Slice<LocalIndex, ConnIndex>& connectivity);


// Implements Newell's method to calculate a polygons normal vector
simple_math::Vector3D calc_polygon_unit_normal(const ElementConstAccessor& polygon_element);

SegmentStats calc_segment_length_stats(const Mesh& m);
SegmentStats calc_segment_length_stats(const ElementConstAccessor& polygon_element);

struct CollapseSegmentsResult
{
    std::vector<std::pair<LocalIndex, LocalIndex>> node_combinations;
    std::vector<LocalIndex> elements_removed;
};

// TODO implement collapse algorithm that respects features
CollapseSegmentsResult collapse_short_polygon_segments(
  Mesh& mesh,
  const double min_edge_length,
  const std::vector<std::pair<LocalIndex, LocalIndex>>& feature_edges,
  const std::vector<LocalIndex>& edge_id_to_feature_ids_map);

/// @brief Refines segments of polygons that have a length larger than the given maximum.
/// @returns A map from the hash of the original segment to the new connectivity of the segment
/// (with the new nodes).
std::unordered_map<std::size_t, std::vector<LocalIndex>> refine_all_polygon_segments(
  Mesh& mesh,
  double length_max);

// Returns the indices of the new triangles
std::vector<LocalIndex> triangulate_polygon(Mesh& total_mesh,
                                            const ElementConstAccessor& polygon_element,
                                            const double max_area     = -1.0,
                                            const bool remove_polygon = false);

std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> triangulate_all_polygons(
  Mesh& total_mesh,
  double max_tri_area = 1.0);

}