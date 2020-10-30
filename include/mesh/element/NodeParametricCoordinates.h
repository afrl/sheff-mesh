#pragma once

#include "mesh/common/ElementTypes.h"
#include "simple_math/EigenDefines.h"
#include "simple_math/quadrature/Quadrature.h"
#include "utility/Exceptions.h"

namespace mesh
{

/// Provides the nodal parametric coordinate given the element type and quadrature type.
/// This relies on the other static functions in NodeParametricCoordinates.h that give
/// the tables for each tyep of element. (definition at end of header)
static Eigen::Map<const simple_math::GeneralMatrix> get_nodal_parametric_points(
  const mesh::ElementType& element_type,
  const simple_math::QuadratureType& quad_type);


//////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions

/// A helper to check if two extents match.  This function exists to allow initialization lists to
/// be used conveniently in the methods below.
static bool extent_matches(const std::array<double, 6>& extent1,
                           const std::array<double, 6>& extent2)
{
    return extent1 == extent2;
}


//////////////////////////////////////////////////////////////////////////////////////////////
// 1D elements

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::L1D2N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_L1D2N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, 0, 0, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[2] = { -1.0, 1.0 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::L1D3N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_L1D3N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, 0, 0, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[3] = { -1.0, 0.0, 1.0 };
    return points;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Linear 2D elements

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::T2D3N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_T2D3N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { 0, 1, 0, 0, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[6] = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::Q2D4N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_Q2D4N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[8] = { -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::C2D4N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_C2D4N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[8] = { -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0 };
    return points;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Quadratic 2D elements

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::T2D6N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_T2D6N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { 0, 1, 0, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[12] = { 0.0, 0.0, 0.5, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 1.0, 0.0, 0.5 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::Q2D8N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_Q2D8N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[16] = { -1.0, -1.0, 0.0, -1.0, 1.0,  -1.0, 1.0,  0.0,
                                 1.0,  1.0,  0.0, 1.0,  -1.0, 1.0,  -1.0, 0.0 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::C2D6N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_C2D6N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[16] = { -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
    return points;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Linear 3D elements

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::T3D4N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_T3D4N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { 0, 1, 0, 1, 0, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::P3D5N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_P3D5N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, -1, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[15] = { -1, -1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 0, 0, 1 };
    return points;
}


/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::W3D6N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_W3D6N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { 0, 1, 0, 1, -1, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[18] = { 0, 0, -1, 1, 0, -1, 0, 1, -1, 0, 0, 1, 1, 0, 1, 0, 1, 1 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::H3D8N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_H3D8N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, -1, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[24] = { -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1,
                                 1,  -1, -1, 1,  1, -1, 1,  1, 1, 1,  -1, 1 };
    return points;
}

// Note: missing Q3D4N

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::C3D8N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_C3D8N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[24] = { -1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0,
                                 -1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0 };
    return points;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Quadratic 3D elements

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::H3D8N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_T3D10N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { 0, 1, 0, 1, 0, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[30] = { 0,   0,   0, 1, 0,   0, 0, 1, 0,   0,   0, 1,   0.5, 0,   0,
                                 0.5, 0.5, 0, 0, 0.5, 0, 0, 0, 0.5, 0.5, 0, 0.5, 0,   0.5, 0.5 };
    return points;
}


/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::P3D13N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_P3D13N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, -1, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[39] = { -1, -1,   -1,   -1, 1,   -1,   1, 1,   -1,  1, -1,   -1,  0,
                                 0,  1,    0,    -1, -1,  1,    0, -1,  0,   1, -1,   -1,  0,
                                 -1, -0.5, -0.5, 0,  0.5, -0.5, 0, 0.5, 0.5, 0, -0.5, 0.5, 0 };
    return points;
}

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::W3D15N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_W3D15N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { 0, 1, 0, 1, -1, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[45] = { 0,   0,   -1, 0.5, 0, -1, 1,   0,   -1, 0.5, 0.5, -1, 0, 1,   -1,
                                 0,   0.5, -1, 0,   0, 0,  1,   0,   0,  0,   1,   0,  0, 0,   1,
                                 0.5, 0,   1,  1,   0, 1,  0.5, 0.5, 1,  0,   1,   1,  0, 0.5, 1 };
    return points;
}


/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::H3D20N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_H3D20N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, -1, 1 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[60] = { -1, -1, -1, -1, 0,  -1, -1, 1,  -1, -1, 1,  0,  -1, 1,  1,
                                 -1, 0,  1,  -1, -1, 1,  -1, -1, 0,  0,  -1, -1, 0,  1,  -1,
                                 0,  1,  1,  0,  -1, 1,  1,  -1, -1, 1,  0,  -1, 1,  1,  -1,
                                 1,  1,  0,  1,  1,  1,  1,  0,  1,  1,  -1, 1,  1,  -1, 0 };
    return points;
}

// Note: missing Q3D8N

/// This return the parametric coordinates of the nodes for elements of type
/// mesh::ElementType::C3D16N as a single array, which can be mapped to
/// a matrix with row-major ordering.
/// Note: an error will be thrown if the coordinates are not hard-coded for
///       the given type of quadrature.
static const double* get_parametric_points_C3D16N(const std::array<double, 6>& parametric_extent)
{
    _UTIL_ASSERT(extent_matches(parametric_extent, { -1, 1, -1, 1, 0, 0 }),
                 "Parametric points not implemented for given type of quadrature.  \
                      Maybe update NodeParametricCoordinates.h for a new type.");

    static double points[48] = { -1, -1, 0,  0, -1, 0,  1, -1, 0,  1,  0,  0, 1,  1,  0, 0,
                                 1,  0,  -1, 1, 0,  -1, 0, 0,  -1, -1, 0,  0, -1, 0,  1, -1,
                                 0,  1,  0,  0, 1,  1,  0, 0,  1,  0,  -1, 1, 0,  -1, 0, 0 };
    return points;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions

static Eigen::Map<const simple_math::GeneralMatrix> get_nodal_parametric_points(
  const mesh::ElementType& element_type,
  const simple_math::QuadratureType& quad_type)
{
    auto parametric_extent = simple_math::get_extent(quad_type);
    switch (element_type)
    {
        // 1D types
        case mesh::ElementType::L1D2N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_L1D2N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              1);                                         // # cols
        case mesh::ElementType::L1D3N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_L1D3N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              1);                                         // # cols
        // Linear 2D types
        case mesh::ElementType::T2D3N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_T2D3N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              2);                                         // # cols
        case mesh::ElementType::Q2D4N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_Q2D4N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              2);                                         // # cols
        case mesh::ElementType::C2D4N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_C2D4N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              2);                                         // # cols
        // Quadratic 2D types
        case mesh::ElementType::T2D6N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_T2D6N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              2);                                         // # cols
        case mesh::ElementType::Q2D8N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_Q2D8N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              2);                                         // # cols
        case mesh::ElementType::C2D6N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_C2D6N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              2);                                         // # cols
        // Linear 2D types
        case mesh::ElementType::T3D4N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_T3D4N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::P3D5N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_P3D5N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::W3D6N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_W3D6N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::H3D8N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_H3D8N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::C3D8N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_C3D8N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        // Quadratic 2D types
        case mesh::ElementType::T3D10N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_T3D10N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::P3D13N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_P3D13N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::W3D15N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_W3D15N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::H3D20N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_H3D20N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        case mesh::ElementType::C3D16N:
            return Eigen::Map<const simple_math::GeneralMatrix>(
              get_parametric_points_C3D16N(parametric_extent),
              mesh::get_num_nodes_for_type(element_type), // # rows
              3);                                         // # cols
        default:
            throw _UTIL_MSG_EXCEPT("No table implemented for parametric coordinates of the nodes \
                                    for the given element type.");
    }
}
}