#pragma once

#include "BasisFunctions.h"
#include "mesh/common/ElementTypes.h"
#include "simple_math/quadrature/Quadrature.h"
#include "utility/DynamicMatrix.h"

#include <map>
#include <vector>

namespace mesh
{
/// @brief A class that caches the basis function and basis function derivatives for different
/// combinations of element type and quadrature.
///
/// When one of these quantities need to be evaluated many times, this class will provide
/// significant speed up. However, if you are just trying to evaulate the basis functions or
/// derivatives at specific point, then use the BasisFunctionsBase class instead.
class BasisFunctionsCache
{
  public:
    BasisFunctionsCache()
      : current_basis_values_(nullptr)
      , current_basis_derivs_(nullptr)
    {
    }

    ~BasisFunctionsCache() {};

    /// Setup the basis functions and derivatives for a specific element and quadrature.
    /// Subsequent calls to values and derivatives will return the
    /// cached values for the last combination of element and quadrature passed to
    /// this method.  Note: the values are only actually calculated once per a combination
    /// of element type, quadrature type, and integration order.  Subsequent calls will
    /// use cached values.
    void setup(const ElementType& element_type, const simple_math::Quadrature& quadrature)
    {
        const auto quad_type = quadrature.type();

        const auto element_type_index =
          static_cast<std::underlying_type<ElementType>::type>(element_type);
        const auto quad_type_index =
          static_cast<std::underlying_type<simple_math::QuadratureType>::type>(quad_type);
        auto int_order_index = quadrature.order() - 1;
        bool force_update    = false;

        if (quadrature.is_adaptive())
        {
            int_order_index = 0;
            force_update    = true;
        }

        current_basis_values_ =
          &std_basis_values_[std::tie(element_type_index, quad_type_index, int_order_index)];
        current_basis_derivs_ =
          &std_basis_derivs_[std::tie(element_type_index, quad_type_index, int_order_index)];
        if (current_basis_values_->size() == 0 || force_update)
        {
            // using EigenRowType = Eigen::DenseBase<simple_math::GeneralMatrix>::RowXpr;

            simple_math::GeneralVector basis_values_for_quad(get_num_nodes_for_type(element_type));

            current_basis_values_->resize(quadrature.n_quad_points(),
                                          get_num_nodes_for_type(element_type));
            current_basis_derivs_->resize(quadrature.n_quad_points());

            for (simple_math::QuadratureIndex qp = 0; qp < quadrature.n_quad_points(); qp++)
            {
                Point parametric_coord = quadrature.parametric_point(qp);
                BasisFunctions::values(element_type, parametric_coord, basis_values_for_quad);
                current_basis_values_->row(qp) = basis_values_for_quad.transpose();
                BasisFunctions::derivatives(element_type,
                                                parametric_coord,
                                                (*current_basis_derivs_)[qp]);
            }
        }
    }

    /// Returns a matrix of the basis function values for the element and quadrature type
    /// indicated in the last setup() call. \see BasisFunctionsCache::setup
    /// The first index of the matrix (row) indicates the quadrature point index, while the
    /// the second index of the matrix (column) indicates the node index in the connectivity.
    /// To interpolate a quantity from the nodes to the quad points, use:
    /// \code [cpp]
    /// BasisFunctionsCache basis_cache;
    /// ... // later in the code, the cache can be used to interpolate:
    /// basis_cache.setup(...);
    /// auto quad_values = basis_cache.values() * nodal_values;
    /// \endcode
    /// where nodal_values is a matrix organized where the rows are the node indices and
    /// the columns are component of the quantity (1 if a scalar quantity).
    const simple_math::GeneralMatrix& values() const
    {
        return *current_basis_values_;
    }
    /// Returns a vector of matrices of the basis function derivatives with respect to the
    /// parametric coordinate system for the element and quadrature type indicated in the last
    /// setup() call. \see BasisFunctionsCache::setup
    /// The index of the vector indicates the quadrature point index.
    /// For each matrix, the first index (row) indicates the node index in the connectivity, while
    ///  the second index (column) indicates the axis that the derivative is with respect to.
    const std::vector<simple_math::GeneralMatrix>& derivatives() const
    {
        return *current_basis_derivs_;
    }

  private:
    using CombinationKey = std::tuple<std::underlying_type<mesh::ElementType>::type,
                                      std::underlying_type<simple_math::QuadratureType>::type,
                                      simple_math::QuadratureOrder>;

    /// Stores the basis functions for at each quadrature point for a range of combinations of
    /// element type, quadrature type, and integration order.  Each entry stores a matrix that
    /// has each basis function value for each quadrature point.  The first index is the quadrature
    /// point index ordered by the Quadrature class and the second index is the basis function index
    /// ordered by the connectivity of the element type.
    /// For all_basis_values_, the first index is the element type, the second index is the
    /// quadrature type, and the third index is (integration_order-1).
    std::map<CombinationKey, simple_math::GeneralMatrix> std_basis_values_;
    /// The basis function values for the current combination of element type, quadrature type, and
    /// integration order given to setup()
    simple_math::GeneralMatrix* current_basis_values_;

    /// Very similar to all_basis_values_.
    /// For entry, the first index (of the vector) is the matrix of the derivatives of the basis
    /// function values for a quadrature point.  The matrix is organized as the first index is the
    /// node/basis index and the second index is the parametric axis that the derivative is taken
    /// with respect to.
    std::map<CombinationKey, std::vector<simple_math::GeneralMatrix>> std_basis_derivs_;

    /// The basis function derivatives for the current combination of element type, quadrature type,
    /// and integration order given to setup()
    std::vector<simple_math::GeneralMatrix>* current_basis_derivs_;
};
}