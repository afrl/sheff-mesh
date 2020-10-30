#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"

using namespace mesh;

int main()
{
    Mesh mesh = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

	write_vtk(mesh, "grid.vtk");

	return 0;
}