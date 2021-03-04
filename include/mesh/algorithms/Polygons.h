#pragma once

#include "mesh/Mesh.h"
#include "simple_math/EigenDefines.h"

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

void triangulate_polygon(Mesh& total_mesh,
                         const ElementConstAccessor& polygon_element,
                         double max_area = std::numeric_limits<double>::max());

}