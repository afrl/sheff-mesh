#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Search.h"
#include "mesh/algorithms/Transformation.h"
#include "simple_math/Utility.h"

using namespace mesh;

int main()
{
    Mesh linear_grid = make_uniform_linear_grid({ 1, 1, 1 }, { 2, 3, 3 }, { 2, 2, 2 });

    scale(linear_grid, 0.5);
    auto extent = find_extent(linear_grid);
    assert(simple_math::point_as_vector(extent.first).isApprox(simple_math::Vector3D(0.5, 0.5, 0.5)));
    assert(
      simple_math::point_as_vector(extent.second).isApprox(simple_math::Vector3D(1, 1.5, 1.5)));

    translate(linear_grid, { 1, 2, 1 });
    extent = find_extent(linear_grid);
    assert(
      simple_math::point_as_vector(extent.first).isApprox(simple_math::Vector3D(1.5, 2.5, 1.5)));
    assert(
      simple_math::point_as_vector(extent.second).isApprox(simple_math::Vector3D(2, 3.5, 2.5)));

    write_vtk(linear_grid, "linear_grid.vtk", false);

    return 0;
}