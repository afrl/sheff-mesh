#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace mesh;

int main()
{
    Mesh mesh = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2, 3, 4 });

    MeshAssociatedData data;
    data.add<DataFieldType::NODE, double>(
      DataField<double>::copy("random_data",
                              { "x", "y" },
                              simple_math::DynamicDoubleMatrix::Random(mesh.nodes().size(), 2)));



    write_vtk(mesh, "grid.vtk", data);

    return 0;
}