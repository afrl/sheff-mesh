#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"

using namespace mesh;

int main()
{
    Mesh linear_grid = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

    write_vtk(linear_grid, "linear_grid.vtk", false);

	assert(linear_grid.nodes().size() == 15);
    assert(linear_grid.elements().size() == 8);

	return 0;
}