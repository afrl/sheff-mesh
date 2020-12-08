#include "mesh/face/FaceHelpers.h"

#include "mesh/face/FaceContainer.h"
#include "mesh/face/FaceIterator.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

#include <set>

using namespace simple_math::point_operators;

namespace mesh
{
Point calc_mean_coordinate(const FaceConstAccessor& face)
{
    Point mean                = { 0, 0, 0 };
    const FloatType num_nodes = static_cast<FloatType>(face.connectivity().size());
    for (auto& node : face.nodes())
    {
        mean = mean + node.point();
    }
    mean = (1.0 / num_nodes) * mean;
    return mean;
}

double calc_signed_norm_projection(const FaceConstAccessor& face, const Point& x)
{
    // This is kind of a hack but it works for all the elements at the moment
    // but for the types we have the normal vector can be found by crossing these two
    // in-plane vectors with each other:
    // vec 1: the first node in the connectivity to the second node
    // vec 2: the first node in the connectivity to the last node
    // I cannot imagine a situation where this would not work, but revaulate if it ever occurs.
    const Point centroid       = calc_mean_coordinate(face);
    const Vector3D centroidVec = simple_math::point_as_vector(centroid);
    const Vector3D vec1        = face.node(1).coord() - face.node(0).coord();
    const Vector3D vec2 =
      face.node(static_cast<ConnIndex>(face.connectivity().size()) - 1).coord() -
      face.node(0).coord();
    const Vector3D normalVec = vec1.cross(vec2).normalized();
    // Center the point relative to the centroid so we get an accurate measure of direction
    const Vector3D pointCentered = simple_math::point_as_vector(x) - centroidVec;
    double projection            = pointCentered.dot(normalVec);
    return projection;
}

bool contains_point(const ElementConstAccessor& element, const Point& x)
{
    FaceContainer faces(element);
    return contains_point(faces, x);
}

bool contains_point(const FaceContainer& faces, const Point& x)
{
    for (auto& face : faces)
    {
        if (calc_signed_norm_projection(face, x) > 0)
        {
            return false;
        }
    }
    return true;
}
}
