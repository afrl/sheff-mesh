#include "mesh/common/QuadratureHelpers.h"

namespace mesh
{
simple_math::QuadratureScheme select_gauss_quadrature_scheme(
  const mesh::ElementType& element_type)
{
    // Set the default quadrature type
    // For now, pyramids will use standard integration, but it should be specialized in
    // the long run see the pyramid element pdf in docs.
    switch (element_type)
    {
        case mesh::ElementType::DEF1D:
        case mesh::ElementType::L1D2N: return { simple_math::QuadratureType::Gauss1D, 1 };
        case mesh::ElementType::L1D3N: return { simple_math::QuadratureType::Gauss1D, 2 };
        case mesh::ElementType::Q2D4N: return { simple_math::QuadratureType::Gauss2D, 1 };
        case mesh::ElementType::Q2D8N: return { simple_math::QuadratureType::Gauss2D, 2 };
        case mesh::ElementType::T2D3N: return { simple_math::QuadratureType::Gauss2DTriangle, 1 };
        case mesh::ElementType::T2D6N: return { simple_math::QuadratureType::Gauss2DTriangle, 2 };
        case mesh::ElementType::H3D8N: return { simple_math::QuadratureType::Gauss3D, 1 };
        case mesh::ElementType::H3D20N: return { simple_math::QuadratureType::Gauss3D, 2 };
        case mesh::ElementType::P3D5N: return { simple_math::QuadratureType::Gauss3D, 1 };
        case mesh::ElementType::P3D13N: return { simple_math::QuadratureType::Gauss3D, 2 };
        case mesh::ElementType::T3D4N: return { simple_math::QuadratureType::Gauss3DTet, 1 };
        case mesh::ElementType::T3D10N: return { simple_math::QuadratureType::Gauss3DTet, 2 };
        case mesh::ElementType::W3D6N: return { simple_math::QuadratureType::Gauss3DWedge, 1 };
        case mesh::ElementType::W3D15N: return { simple_math::QuadratureType::Gauss3DWedge, 2 };
        default: throw _UTIL_MSG_EXCEPT("Unknown Element Type.");
    }
}

simple_math::QuadratureScheme select_gauss_quadrature_scheme(const mesh::FaceType& face_type)
{
    switch (face_type)
    {
        case mesh::FaceType::L2D2N: return { simple_math::QuadratureType::Gauss1D, 1 };
        case mesh::FaceType::L2D3N: return { simple_math::QuadratureType::Gauss1D, 2 };
        case mesh::FaceType::Q3D4N: return { simple_math::QuadratureType::Gauss2D, 1 };
        case mesh::FaceType::Q3D8N: return { simple_math::QuadratureType::Gauss2D, 2 };
        case mesh::FaceType::T3D3N: return { simple_math::QuadratureType::Gauss2DTriangle, 1 };
        case mesh::FaceType::T3D6N: return { simple_math::QuadratureType::Gauss2DTriangle, 2 };
        default: throw _UTIL_MSG_EXCEPT("Unknown Element Type.");
    }
}
}