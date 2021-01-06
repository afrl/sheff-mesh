#include "mesh/algorithms/Transformation.h"

#include "simple_math/Utility.h"

namespace mesh
{
Mesh& translate(Mesh& mesh, const simple_math::Point& translation_vector)
{
    auto t = simple_math::point_as_vector(translation_vector);
    for (auto& node : mesh.nodes())
    {
        node.coord() += t;
    }
    return mesh;
}

Mesh& scale(Mesh& mesh, const double& scale_factor, const simple_math::Point& origin)
{
    auto o = simple_math::point_as_vector(origin);
    for (auto& node : mesh.nodes())
    {
        node.coord() = scale_factor * (node.coord() - o) + o;
    }
    return mesh;
}
Mesh& scale(Mesh& mesh,
            const std::array<double, 3>& scale_factors,
            const simple_math::Point& origin)
{
    for (auto& node : mesh.nodes())
    {
        for (auto i : utility::range(3))
        {
            node.point(i) = scale_factors[i] * (node.point(i) - origin[i]) + origin[i];
        }
    }
    return mesh;
}
}