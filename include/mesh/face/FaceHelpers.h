#pragma once

#include "FaceAccessor.h"

namespace mesh
{
/// @brief Returns a mean coordinate of the nodes which should be interior to the face in most
/// cases.
Point calc_mean_coordinate(const FaceConstAccessor& face);

/// @brief Calculate the signed normal projection of a given point on the given face.
///
/// @return Signed distance from face with positive indicating the point is outside
/// of element, and negative indicating the point is inside the element.
double calc_signed_norm_projection(const FaceConstAccessor& face, const Point& x);

/// @brief Determines if an element contains a given point based on projection on the faces.
bool contains_point(const ElementConstAccessor& element, const Point& x);

/// @brief Equivalent to other contains_point() function, but requires the FaceContainer of the
/// element to be checked.
bool contains_point(const FaceContainer& faces, const Point& x);
}