#include "mesh/surface/SurfaceElementHelpers.h"

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
}