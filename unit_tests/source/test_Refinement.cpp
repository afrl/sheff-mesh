#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/algorithms/Refinement.h"
#include "mesh/algorithms/Search.h"
#include "mesh/algorithms/UnitElements.h"

using namespace mesh;

struct RefinementCase
{
    ElementType element_type;
    int n_faces;
    int refine_factor;
};

int main()
{
    std::vector<RefinementCase> refine_cases = { { ElementType::T3D4N, 4, 4 },
                                                 { ElementType::T3D10N, 4, 4 },
                                                 { ElementType::H3D8N, 6, 4 },
                                                 { ElementType::H3D20N, 6, 4 } };
    for (auto& refine_case : refine_cases)
    {
        Mesh unit_element_mesh                         = get_unit_element(refine_case.element_type);
        unit_element_mesh.elements()[0].component_id() = 1;
        auto& group                                    = unit_element_mesh.elements().group("test");
        group.element_indices.push_back(0);
        write_vtk(unit_element_mesh, to_string(refine_case.element_type) + "_0.vtk", false);
        assert(extract_boundary(unit_element_mesh).surface.size() == refine_case.n_faces);

        refine_uniform(unit_element_mesh);
        write_vtk(unit_element_mesh, to_string(refine_case.element_type) + "_1.vtk", false);
        assert(extract_boundary(unit_element_mesh).surface.size() ==
               refine_case.n_faces * refine_case.refine_factor);
        for (auto e : unit_element_mesh.elements())
        {
            assert(e.component_id() == 1);
            assert(std::find(group.element_indices.begin(),
                             group.element_indices.end(),
                             e.local_index()) != group.element_indices.end());
        }

        refine_uniform(unit_element_mesh);
        write_vtk(unit_element_mesh, to_string(refine_case.element_type) + "_2.vtk", false);
        assert(extract_boundary(unit_element_mesh).surface.size() ==
               refine_case.n_faces * refine_case.refine_factor * refine_case.refine_factor);
        for (auto e : unit_element_mesh.elements())
        {
            assert(e.component_id() == 1);
            assert(std::find(group.element_indices.begin(),
                             group.element_indices.end(),
                             e.local_index()) != group.element_indices.end());
        }
    }

    return 0;
}