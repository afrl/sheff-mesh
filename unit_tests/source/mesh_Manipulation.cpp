#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/algorithms/IO.h"

using namespace mesh;

int main()
{
    Mesh grid = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

	assert(grid.nodes().size() == 15);
    assert(grid.elements().size() == 8);

    convert_to_quadratic(grid);

    write_vtk(grid, "grid.vtk", false);
    
    assert(grid.nodes().size() == 37);
    assert(grid.elements().size() == 8);

	return 0;
}