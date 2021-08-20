#include "mesh/algorithms/UnitElements.h"
#include "mesh/face/Face.h"

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace mesh;

int main()
{
    Mesh mesh = get_unit_element(mesh::ElementType::H3D8N);

    ElementConstAccessor element = mesh.elements()[0];
    FaceContainer faces(element);

    // Checks begin/end iterators for range-based loop and face types
    for (auto& face : faces)
    {
        assert(face.type() == FaceType::Q3D4N);

        face.node(0).point(); // make sure node accessor doesn't throw

        for (auto& node : face.nodes())
        {
            node.point(); // make sure node accessor doesn't throw
        }

        // Ensure face lies on coordinate plane
        assert(std::all_of(face.nodes().begin(),
                            face.nodes().end(),
                            [](NodeConstAccessor n) -> bool { return n.x() == 1.0; }) ||
               std::all_of(face.nodes().begin(),
                            face.nodes().end(),
                           [](NodeConstAccessor n) -> bool { return n.x() == -1.0; }) ||
               std::all_of(face.nodes().begin(),
                            face.nodes().end(),
                           [](NodeConstAccessor n) -> bool { return n.y() == 1.0; }) ||
               std::all_of(face.nodes().begin(),
                            face.nodes().end(),
                           [](NodeConstAccessor n) -> bool { return n.y() == -1.0; }) ||
               std::all_of(face.nodes().begin(),
                            face.nodes().end(),
                           [](NodeConstAccessor n) -> bool { return n.z() == 1.0; }) ||
               std::all_of(face.nodes().begin(),
                            face.nodes().end(),
                            [](NodeConstAccessor n) -> bool { return n.z() == -1.0; }));
    }

    // Test the helpers.  Note: contains_point calls:
    // - calc_mean_coordinate
    // - calc_signed_norm_projection
    // So we can weakly test all of those at once
    assert(contains_point(element, {0, 0, 0}));
    assert(contains_point(element, { 1, 0, 0 }));
    assert(contains_point(element, { 1, 1, 1 }));
    assert(!contains_point(element, { 2, 2, 2 }));
    assert(!contains_point(element, { 2, 0, 0 }));

    return 0;
}