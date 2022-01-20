#include "unit_test_helper.h"
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace mesh;

int main()
{
    auto fiber_matrix_mesh = read_pfec(get_test_data_filepath("FiberMatrixRVE.pfecx"));

    write_vtk(fiber_matrix_mesh.first, "FiberMatrixRVE.vtk");

    return 0;
}