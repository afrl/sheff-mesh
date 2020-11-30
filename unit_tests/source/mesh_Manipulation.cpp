#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Manipulation.h"

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

    Mesh tri;
    tri.nodes().add({ 0, 0, 0 });
    tri.nodes().add({ 1, 0, 0 });
    tri.nodes().add({ 0, 1, 0 });
    tri.elements().add(ElementType::T2D3N, { 0, 1, 2 });

    convert_to_quadratic(tri);

    assert(tri.elements()[0].type() == ElementType::T2D6N);
    assert(tri.elements()[0].connectivity()[0] == 0);
    assert(tri.elements()[0].connectivity()[2] == 1);
    assert(tri.elements()[0].connectivity()[4] == 2);

    write_vtk(tri, "quadratic_tri.vtk", false);

    Mesh tet;
    tet.nodes().add({ 0, 0, 0 });
    tet.nodes().add({ 1, 0, 0 });
    tet.nodes().add({ 0, 1, 0 });
    tet.nodes().add({ 0, 0, 1 });
    tet.elements().add(ElementType::T3D4N, { 0, 1, 2, 3 });
    write_vtk(tet, "linear_tet.vtk", false);

    convert_to_quadratic(tet);

    assert(tri.elements()[0].type() == ElementType::T3D10N);
    assert(tri.elements()[0].connectivity()[0] == 0);
    assert(tri.elements()[0].connectivity()[1] == 1);
    assert(tri.elements()[0].connectivity()[2] == 2);
    assert(tri.elements()[0].connectivity()[3] == 3);

    write_vtk(tet, "quadratic_tet.vtk", false);

    return 0;
}