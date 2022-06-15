#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/algorithms/IO.h"

using namespace mesh;

int main()
{
    Mesh grid = make_uniform_linear_tet_grid({ 0, 0, 0 }, {9, 1, 1}, { 9, 1, 1 });

    assert(grid.nodes().size() == 40);
    assert(grid.elements().size() == 6 * 9);

    convert_to_quadratic(grid);

    write_vtk(grid, "quadratic_tet_grid.vtk", false);

	return 0;
}