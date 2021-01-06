#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Search.h"

using namespace mesh;

int main()
{
    Mesh mixed_grid = make_mixed_element_grid({ 0, 0, 0 }, { 1, 1, 1 }, { 3, 3, 3 });
    auto n_blocks   = 3 * 3 * 3;

    write_vtk(mixed_grid, "mixed_grid.vtk", false);

    assert(mixed_grid.elements().size() == 92 * n_blocks);

    // Make sure there is the correct amount of each element type
    auto type_counts = count_element_types(mixed_grid);
    assert(type_counts[ElementType::P3D5N] == 50 * n_blocks);
    assert(type_counts[ElementType::T3D4N] == 20 * n_blocks);
    assert(type_counts[ElementType::W3D6N] == 10 * n_blocks);
    assert(type_counts[ElementType::H3D8N] == 12 * n_blocks);

    // Make sure there are no extra unmatched faces internally
    auto surface = extract_boundary(mixed_grid);
    write_vtk(surface.mesh(), "mixed_grid_surface.vtk", false);
    assert(surface.size() == 3 * 18 * 9);

    return 0;
}