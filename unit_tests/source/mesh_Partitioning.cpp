#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Partitioning.h"

using namespace mesh;

int main()
{
    std::vector<double> coordinates(5);
    std::iota(coordinates.begin(), coordinates.end(), 0.0);

    Mesh mesh = make_linear_grid(coordinates, coordinates, coordinates);

    partition_metis(mesh, 3);
    write_vtk(mesh, "total_mesh.vtk");

    auto mesh_part0 = extract_partition_mesh(mesh, 0);
    write_vtk(mesh_part0, "partition_0.vtk");

    return 0;
}