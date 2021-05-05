#include "simple_math/quadrature/QuadratureDefines.h"
#include "mesh/common/FaceTypes.h"
#include "mesh/common/ElementTypes.h"

namespace mesh
{
/// @brief Selects the appropriate Gauss quadrature scheme to use for a solid element.
///
/// Note: the integration order is selected such that a constant integrand is exactly integrated for the element.
simple_math::QuadratureScheme select_gauss_quadrature_scheme(const ElementType& element_type);

/// @brief Selects the appropriate Gauss quadrature scheme to use for a face.
///
/// Note: the integration order is selected such that a constant integrand is exactly integrated for the face.
simple_math::QuadratureScheme select_gauss_quadrature_scheme(const FaceType& face_type);
}