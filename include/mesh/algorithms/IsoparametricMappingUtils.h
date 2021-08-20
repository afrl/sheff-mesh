#pragma once

#include "mesh/Mesh.h"

#include <utility>

namespace mesh
{
/// @brief Given an element, returns the global coordinates of the nodes for the element.
/// @returns A Nx3 matrix (where N is the number of nodes in the element) with each row
///          in the matrix ordered by the element's connectivity.
simple_math::GeneralMatrix get_nodal_global_coords(const ElementConstAccessor& element);

/// @brief Given an element and parametric point, returns the corresponding global coordinate.
///
/// Note: if you are trying to map fixed quadrature points, then this will be inefficient, since
///       the BasisFunctionCache is not used.  Only use this for arbitrary parametric points.
/// Note: if this is called more than once for an element, use the other definition that takes
///       nodal_global_coords.
simple_math::Vector3D parametric_to_global(const ElementConstAccessor& element,
                                           const simple_math::Point& parametric_point);

/// @brief Given an element and parametric point, returns the corresponding global coordinate.
///
/// Note: if you are trying to map fixed quadrature points, then this will be inefficient, since
///       the BasisFunctionCache is not used.  Only use this for arbitrary parametric points.
simple_math::Vector3D parametric_to_global(const ElementConstAccessor& element,
                                           const simple_math::Point& parametric_point,
                                           const simple_math::GeneralMatrix& nodal_global_coords);

/// @brief Given an element and parametric point, returns the corresponding global coordinate.
///
/// Note: if you are trying to map fixed quadrature points, then this will be inefficient, since
///       the BasisFunctionCache is not used.  Only use this for arbitrary parametric points.
/// Note: if this is called more than once for an element, use the other definition that takes
///       nodal_global_coords.
simple_math::Vector3D parametric_to_global_vec(const ElementConstAccessor& element,
                                               const simple_math::Vector3D& parametric_point);

/// @brief Given an element and parametric point, returns the corresponding global coordinate.
///
/// Note: if you are trying to map fixed quadrature points, then this will be inefficient, since
///       the BasisFunctionCache is not used.  Only use this for arbitrary parametric points.
simple_math::Vector3D parametric_to_global_vec(
  const ElementConstAccessor& element,
  const simple_math::Vector3D& parametric_point,
  const simple_math::GeneralMatrix& nodal_global_coords);

/// @brief Given an element and a global point, returns the corresponding parametric coordinate
///        in the element.
///
/// For irregular elements (non-parallel faces), the search element may fail if the global point
/// is outside the element.  Consequently, check the bool in the returned pair before trusting
/// the resulting parametric point.
///
/// Note: if this is called more than once for an element, use the other definition that takes
///       nodal_global_coords.
///
/// @returns Pair of (first) whether the search was successful and (second) the parametric point.
std::pair<bool, simple_math::Vector3D> global_to_parametric(const ElementConstAccessor& element,
                                                            const simple_math::Point& global_point);

/// @brief Given an element and a global point, returns the corresponding parametric coordinate
///        in the element.
///
/// For irregular elements (non-parallel faces), the search element may fail if the global point
/// is outside the element.  Consequently, check the bool in the returned pair before trusting
/// the resulting parametric point.
///
/// @returns Pair of (first) whether the search was successful and (second) the parametric point.
std::pair<bool, simple_math::Vector3D> global_to_parametric(
  const ElementConstAccessor& element,
  const simple_math::Point& global_point,
  const simple_math::GeneralMatrix& nodal_global_coords);

/// @brief Given an element and a global point, returns the corresponding parametric coordinate
///        in the element.
///
/// For irregular elements (non-parallel faces), the search element may fail if the global point
/// is outside the element.  Consequently, check the bool in the returned pair before trusting
/// the resulting parametric point.
///
/// Note: if this is called more than once for an element, use the other definition that takes
///       nodal_global_coords.
///
/// @returns Pair of (first) whether the search was successful and (second) the parametric point.
std::pair<bool, simple_math::Vector3D> global_to_parametric_vec(
  const ElementConstAccessor& element,
  const simple_math::Vector3D& global_point);

/// @brief Given an element and a global point, returns the corresponding parametric coordinate
///        in the element.
///
/// For irregular elements (non-parallel faces), the search element may fail if the global point
/// is outside the element.  Consequently, check the bool in the returned pair before trusting
/// the resulting parametric point.
///
/// @returns Pair of (first) whether the search was successful and (second) the parametric point.
std::pair<bool, simple_math::Vector3D> global_to_parametric_vec(
  const ElementConstAccessor& element,
  const simple_math::Vector3D& global_point,
  const simple_math::GeneralMatrix& nodal_global_coords);
}