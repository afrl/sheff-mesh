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
    auto mesh = read_abaqus(get_test_data_filepath("BarWithHole.inp"));

    assert(mesh.nodes().size() == 5724);
    assert(mesh.elements().size() == 3596);

    write_vtk(mesh, "BarWithHole.vtk");

    auto mesh2 = read_abaqus(get_test_data_filepath("BoxAndCylinder.inp"));
    write_vtk(mesh2, "BoxAndCylinder.vtk");

    return 0;
}