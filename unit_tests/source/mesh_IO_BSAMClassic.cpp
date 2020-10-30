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
    auto mesh = read_bsam_classic(get_test_data_filepath("SquarePanel_Elast.in"));

    assert(mesh.nodes().size() == 512 * 4);
    assert(mesh.elements().size() == 225 * 4);
    assert(mesh.nodes().groups().size() == 7);
    assert(mesh.elements().groups().size() == 4);

    write_vtk(mesh, "SquarePanel_Elast.vtk");

    return 0;
}