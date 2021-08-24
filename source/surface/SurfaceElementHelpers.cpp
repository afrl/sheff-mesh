#include "mesh/surface/SurfaceElementHelpers.h"

#include "mesh/algorithms/Polygons.h" // polygon_segments
#include "simple_math/Utility.h"

namespace mesh
{
Point calc_mean_coordinate_surface_element(const SurfaceElementAccessor& element)
{
    using namespace simple_math::point_operators;

    Point mean                = { 0, 0, 0 };
    const FloatType num_nodes = static_cast<FloatType>(element.connectivity().size());
    for (auto& node : element.nodes())
    {
        mean = mean + node.point();
    }
    mean = (1.0 / num_nodes) * mean;
    return mean;
}

simple_math::Vector3D calc_outward_normal(const SurfaceElementAccessor& surf_element)
{
    // Implements Newell's method to calculate a polygons normal vector
    // Assumes that nodes are orderd counter-clockwise around face, which is true for the types
    // below, but algorithm may need tweaked for other types.
    _UTIL_ASSERT(
      surf_element.type() == ElementType::T2D3N || surf_element.type() == ElementType::T2D6N ||
        surf_element.type() == ElementType::Q2D4N || surf_element.type() == ElementType::Q2D8N,
      "Algorithm used to calculate the outward normal not implemented for element type.");

    simple_math::Vector3D normal = simple_math::Vector3D::Zero();
    for (auto segment : polygon_segments(surf_element.connectivity()))
    {
        auto p0 = surf_element.node(segment.first).coord();
        auto p1 = surf_element.node(segment.second).coord();
        normal += p0.cross(p1);
    }
    normal.normalize();
    return normal;
}
}