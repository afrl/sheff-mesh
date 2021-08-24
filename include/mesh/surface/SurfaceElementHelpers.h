#pragma once

#include "mesh/surface/SurfaceElementAccessor.h"

namespace mesh
{
/// Returns a mean coordinate of the nodes which should be interior to the element in most cases.
Point calc_mean_coordinate_surface_element(const SurfaceElementAccessor& surf_element);

simple_math::Vector3D calc_outward_normal(const SurfaceElementAccessor& surf_element);
}