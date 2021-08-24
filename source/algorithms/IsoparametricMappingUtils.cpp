#include "mesh/algorithms/IsoparametricMappingUtils.h"

#include "mesh/basis/BasisFunctions.h"
#include "simple_math/Utility.h"

namespace mesh
{
//////////////////////////////
// Mapping for an element

simple_math::GeneralMatrix get_nodal_global_coords(const ElementConstAccessor& element)
{
    simple_math::GeneralMatrix node_coords(element.connectivity().size(), 3);
    for (auto i : utility::range(element.connectivity().size()))
    {
        node_coords.row(i) = element.node(i).coord();
    }
    return node_coords;
}

simple_math::Vector3D parametric_to_global(const ElementConstAccessor& element,
                                           const simple_math::Point& parametric_point)
{
    auto nodal_global_coords = get_nodal_global_coords(element);
    return parametric_to_global(element, parametric_point, nodal_global_coords);
}

simple_math::Vector3D parametric_to_global(const ElementConstAccessor& element,
                                           const simple_math::Point& parametric_point,
                                           const simple_math::GeneralMatrix& nodal_global_coords)
{
    simple_math::GeneralVector basis_values;
    BasisFunctions::values(element.type(), parametric_point, basis_values);
    return nodal_global_coords.transpose() * basis_values;
}

simple_math::Vector3D parametric_to_global_vec(const ElementConstAccessor& element,
                                               const simple_math::Vector3D& parametric_point)
{
    auto nodal_global_coords = get_nodal_global_coords(element);
    return parametric_to_global_vec(element, parametric_point, nodal_global_coords);
}

simple_math::Vector3D parametric_to_global_vec(
  const ElementConstAccessor& element,
  const simple_math::Vector3D& parametric_point,
  const simple_math::GeneralMatrix& nodal_global_coords)
{
    simple_math::GeneralVector basis_values;
    BasisFunctions::values(element.type(), simple_math::to_point(parametric_point), basis_values);
    return nodal_global_coords.transpose() * basis_values;
}

std::pair<bool, simple_math::Vector3D> global_to_parametric(const ElementConstAccessor& element,
                                                            const simple_math::Point& global_point)
{
    auto nodal_global_coords = get_nodal_global_coords(element);
    return global_to_parametric_vec(element,
                                    simple_math::point_as_vector(global_point),
                                    nodal_global_coords);
}

std::pair<bool, simple_math::Vector3D> global_to_parametric(
  const ElementConstAccessor& element,
  const simple_math::Point& global_point,
  const simple_math::GeneralMatrix& nodal_global_coords)
{
    return global_to_parametric_vec(element,
                                    simple_math::point_as_vector(global_point),
                                    nodal_global_coords);
}

std::pair<bool, simple_math::Vector3D> global_to_parametric_vec(
  const ElementConstAccessor& element,
  const simple_math::Vector3D& global_point)
{
    auto nodal_global_coords = get_nodal_global_coords(element);
    return global_to_parametric_vec(element, global_point, nodal_global_coords);
}

std::pair<bool, simple_math::Vector3D> global_to_parametric_vec(
  const ElementConstAccessor& element,
  const simple_math::Vector3D& global_coord,
  const simple_math::GeneralMatrix& nodal_global_coords)
{
    simple_math::GeneralVector basis_values;
    simple_math::GeneralMatrix basis_derivatives;
    simple_math::Vector3D curr_guess = simple_math::Vector3D::Zero();
    int curr_iteration               = 0;
    do
    {
        BasisFunctionsVec::values(element.type(), curr_guess, basis_values);
        auto residual = nodal_global_coords.transpose() * basis_values - global_coord;
        if (residual.norm() < 1e-8)
        {
            return { true, curr_guess };
        }

        BasisFunctionsVec::derivatives(element.type(), curr_guess, basis_derivatives);
        auto jacobian = nodal_global_coords.transpose() * basis_derivatives;

        curr_guess -= jacobian.inverse() * residual;
        curr_iteration++;
    } while (curr_iteration < 20);

    return { false, curr_guess };
}

//////////////////////////////////////
// Mappings for a surface element

simple_math::GeneralMatrix get_nodal_global_coords(const SurfaceElementAccessor& element)
{
    simple_math::GeneralMatrix node_coords(element.connectivity().size(), 3);
    for (auto i : utility::range(element.connectivity().size()))
    {
        node_coords.row(i) = element.node(i).coord();
    }
    return node_coords;
}

simple_math::Vector3D parametric_to_global(const SurfaceElementAccessor& element,
                                           const simple_math::Point& parametric_point)
{
    auto nodal_global_coords = get_nodal_global_coords(element);
    return parametric_to_global(element, parametric_point, nodal_global_coords);
}

simple_math::Vector3D parametric_to_global(const SurfaceElementAccessor& element,
                                           const simple_math::Point& parametric_point,
                                           const simple_math::GeneralMatrix& nodal_global_coords)
{
    simple_math::GeneralVector basis_values;
    BasisFunctions::values(element.type(), parametric_point, basis_values);
    return nodal_global_coords.transpose() * basis_values;
}

simple_math::GeneralMatrix quadrature_parametric_to_global(const SurfaceElementAccessor& element,
                                                           const simple_math::Quadrature quadrature)
{
    BasisFunctionsCache basis_cache;
    return quadrature_parametric_to_global(element, quadrature, basis_cache);
}

simple_math::GeneralMatrix quadrature_parametric_to_global(const SurfaceElementAccessor& element,
                                                           const simple_math::Quadrature quadrature,
                                                           BasisFunctionsCache& basis_cache)
{
    basis_cache.setup(element.type(), quadrature);
    return basis_cache.values() * quadrature.parametric_points();
}

}