#include "mesh/Mesh.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/IsoparametricMappingUtils.h"
#include "simple_math/Utility.h"

using namespace mesh;

int main()
{
    // TODO implement tests for other parts of search algorithms

    {
        Mesh mesh;
        mesh.nodes().add({ 3.889, 0.041, 10.223 });
        mesh.nodes().add({ 3.991, 0.038, 10.223 });
        mesh.nodes().add({ 3.937, 0.199, 10.224 });
        mesh.nodes().add({ 3.811, 0.186, 10.224 });
        mesh.nodes().add({ 3.895, 0.041, 10.446 });
        mesh.nodes().add({ 3.996, 0.039, 10.446 });
        mesh.nodes().add({ 3.943, 0.202, 10.446 });
        mesh.nodes().add({ 3.818, 0.187, 10.446 });

        mesh.elements().add(ElementType::H3D8N, { 0, 1, 2, 3, 4, 5, 6, 7 });
        auto element = mesh.elements().back();

        auto xis         = simple_math::linspace(-0.1, 1.1, 11);
        auto node_coords = get_nodal_global_coords(element);

        for (auto i : utility::range(xis.size()))
        {
            for (auto j : utility::range(xis.size()))
            {
                for (auto k : utility::range(xis.size()))
                {
                    auto global_point =
                      parametric_to_global(element, { xis[i], xis[j], xis[j] }, node_coords);
                    auto result = global_to_parametric_vec(element, global_point, node_coords);
                    assert(result.first == true);
                    auto global_point2 =
                      parametric_to_global_vec(element, result.second, node_coords);
                    assert((global_point2 - global_point).norm() < 1e-6);
                }
            }
        }
    }

    return 0;
}