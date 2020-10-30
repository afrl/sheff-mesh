#include "mesh/algorithms/Creation.h"

using namespace mesh;

int main()
{
    Mesh linear_grid = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

	assert(linear_grid.nodes().size() == 15);
    assert(linear_grid.elements().size() == 8);

    Mesh quadratic_grid = make_quadratic_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

    assert(quadratic_grid.nodes().size() == 33);
    assert(quadratic_grid.elements().size() == 8);

	return 0;
}