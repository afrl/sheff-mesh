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


    Mesh quadratic_block                         = get_unit_element(ElementType::T3D10N);
    quadratic_block.elements()[0].component_id() = 1;
    auto& group                                   = quadratic_block.elements().group("test");
    group.element_indices.push_back(0);
    assert(extract_boundary(quadratic_block).surface.size() == 4);

    refine_uniform(quadratic_block);
    assert(extract_boundary(quadratic_block).surface.size() == 4 * 4);
    for (auto e : quadratic_block.elements())
    {
        assert(e.component_id() == 1);
        assert(std::find(group.element_indices.begin(),
                         group.element_indices.end(),
                         e.local_index()) != group.element_indices.end());
    }

    refine_uniform(quadratic_block);
    assert(extract_boundary(quadratic_block).surface.size() == 4 * 4 * 4);
    for (auto e : quadratic_block.elements())
    {
        assert(e.component_id() == 1);
        assert(std::find(group.element_indices.begin(),
                         group.element_indices.end(),
                         e.local_index()) != group.element_indices.end());
    }

    return 0;
}