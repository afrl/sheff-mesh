#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Partitioning.h"
#include "utility/MPIUtils.h"

#include <iostream>

using namespace mesh;

int main(int argc, char** argv)
{
    utility::MPIUtils::init(&argc, &argv);

    std::vector<double> coordinates(5);
    std::iota(coordinates.begin(), coordinates.end(), 0.0);

    Mesh mesh = make_linear_grid(coordinates, coordinates, coordinates);

    // Test read and write with no meta file
    {
        write_fea_io(mesh, ".", "grid_mesh_no_meta");
        auto grid_mesh = read_fea_io_mesh("./grid_mesh_no_meta", 0);
        assert(grid_mesh.nodes().size() == mesh.nodes().size());
        assert(grid_mesh.elements().size() == mesh.elements().size());
    }

    // Test read and write with a meta file
    {
        write_fea_io(mesh, ".", "grid_mesh", "grid_mesh");
        auto grid_mesh = read_fea_io_mesh_from_meta("./grid_mesh.x.json", 0);
        assert(grid_mesh.nodes().size() == mesh.nodes().size());
        assert(grid_mesh.elements().size() == mesh.elements().size());
    }

    partition_metis(mesh, 3);

    // Note: you cannot check that each partition based on previous results, since
    // METIS is not deterministic across platforms.

    // Test read and write of partitioned mesh with no meta file
    {
        write_fea_io(mesh, ".", "grid_mesh_partitioned_no_meta");
        auto grid_mesh = read_fea_io_mesh("./grid_mesh_partitioned_no_meta", 0);
        assert(grid_mesh.nodes().size() < mesh.nodes().size());
        assert(grid_mesh.elements().size() < mesh.elements().size());
    }

    // Test read and write of partitioned mesh with no meta file
    {
        write_fea_io(mesh, ".", "grid_mesh_partitioned", "grid_mesh_partitioned");
        auto grid_mesh = read_fea_io_mesh_from_meta("./grid_mesh_partitioned.x.json", 0);
        assert(grid_mesh.nodes().size() < mesh.nodes().size());
        assert(grid_mesh.elements().size() < mesh.elements().size());
    }

    return 0;
}