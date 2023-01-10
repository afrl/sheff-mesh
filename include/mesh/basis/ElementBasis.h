#pragma once

#include "BasisFunctions.h"
#include "mesh/element/ElementAccessor.h"

namespace mesh
{
/// @brief Calculates the jacobian matrix at the given parametric point.
///
/// The jacobian matrix is in the form:
/// \code
/// J = [ dx/dxi    dy/dxi    dz/dxi   ]
/// 	[ dx/deta   dy/deta   dz/deta  ]
/// 	[ dx/dzeta  dy/dzeta  dz/dzeta ]
/// \endcode
/// @return Returns a 3x3 matrix regardless of the dimensionality of the element, though if the
/// element is 2D on the first two rows and columns will be filled and if the element is 1D then
/// just the (0,0) entry will be filled.
simple_math::Matrix3D calc_jacobian(const ElementConstAccessor element,
                                    const simple_math::Point& parametric_coordinate);
}