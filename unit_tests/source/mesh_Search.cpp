#include "mesh/Mesh.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Search.h"
#include "simple_math/Utility.h"

using namespace mesh;

int main()
{
    // TODO implement tests for other parts of search algorithms

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

    write_vtk(mesh, "mesh.vtk", false);

    auto p = find_parametric_point(mesh.elements()[0],
                                   { 3.9152717844331821, 0.19828546882018369, 10.399085879880049 });

    assert((simple_math::point_as_vector(p) -
            simple_math::Vector3D(0.57735026918826038, 0.57735027312199760, 1.0000000014369588))
             .norm() < 1e-10);

    return 0;
}