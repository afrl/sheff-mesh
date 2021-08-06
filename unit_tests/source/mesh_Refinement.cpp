#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/algorithms/Refinement.h"
#include "mesh/algorithms/Search.h"
#include "mesh/algorithms/UnitElements.h"

using namespace mesh;

int main()
{
    Mesh linear_block = get_unit_element(ElementType::T3D4N);
    assert(extract_boundary(linear_block).surface.size() == 4);

    refine_uniform(linear_block);
    assert(extract_boundary(linear_block).surface.size() == 4 * 4);

    refine_uniform(linear_block);
    assert(extract_boundary(linear_block).surface.size() == 4 * 4 * 4);


    Mesh quadratic_block = get_unit_element(ElementType::T3D10N);
    assert(extract_boundary(quadratic_block).surface.size() == 4);

    refine_uniform(quadratic_block);
    assert(extract_boundary(quadratic_block).surface.size() == 4 * 4);

    refine_uniform(quadratic_block);
    assert(extract_boundary(quadratic_block).surface.size() == 4 * 4 * 4);

    return 0;
}