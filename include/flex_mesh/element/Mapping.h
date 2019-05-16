#pragma once

#include "simple_math/quadrature/Quadrature.h"
#include "utility/DynamicMatrix.h"

#include "Element.h"
#include "ShapeFunctions.h"

#include <vector>

namespace flex_mesh
{
	/// A class that caches the shape function and shape function derivatives
	/// for different combinations of element type and quadrature.  When one of these quantities
	/// need to be evaluated many times, this class will provide significant speed up.
	/// However, if you are just trying to evaulate the shape functions or derivatives at specific
	/// point, then use the Element::get_shape_function_values and
	/// Element::get_shape_function_derivatives methods instead.
	class Mapping
	{
	public:
		Mapping()
			: current_shape_values_(all_shape_values_[0][0][0]),
			current_shape_derivs_(all_shape_derivs_[0][0][0])
		{
		}

		~Mapping() {};

		/// Setup the mapping for a specific element and quadrature.
		/// Subsequent calls to get_shape_values and get_shape_derivatives will return the
		/// cached values for the last combination of element and quadrature passed to
		/// this method.  Note: the values are only actually calculated once per a combination
		/// of element type, quadrature type, and integration order.  Subsequent calls will
		/// use cached values.
		void setup(const ElementType& element_type, const simple_math::Quadrature& quadrature)
		{
			const auto element_type_index = static_cast<std::underlying_type<ElementType>::type>(element_type);
			const auto quad_type_index = static_cast<std::underlying_type<simple_math::QuadratureType>::type>(quadrature.get_type());
			const auto int_order_index = quadrature.get_integration_order() - 1;

			current_shape_values_ = all_shape_values_[element_type_index][quad_type_index][int_order_index];
			if (current_shape_values_.size() == 0)
			{
				//using EigenRowType = Eigen::DenseBase<simple_math::GeneralMatrix>::RowXpr;
				using ShapeFunctionsType = ShapeFunctions<Point, simple_math::GeneralVector, simple_math::GeneralMatrix>;

				simple_math::GeneralVector shape_values_for_quad(get_num_nodes_for_type(element_type));

				current_shape_values_.resize(quadrature.n_quad_points(), get_num_nodes_for_type(element_type));
				current_shape_derivs_.resize(quadrature.n_quad_points());

				for (unsigned int qp = 0; qp < quadrature.n_quad_points(); qp++)
				{
					Point parametric_coord = quadrature.parametric_point(qp);
					ShapeFunctionsType::values(element_type, parametric_coord, shape_values_for_quad);
					current_shape_values_.row(qp) = shape_values_for_quad.transpose();
					ShapeFunctionsType::derivatives(element_type, parametric_coord, current_shape_derivs_[qp]);
				}
			}
		}

		/// Returns a matrix of the shape function values for the element and quadrature type
		/// indicated in the last setup() call. \see Mapping::setup
		/// The first index of the matrix (row) indicates the quadrature point index, while the
		/// the second index of the matrix (column) indicates the node index in the connectivity.
		/// To get a quantity mapped from the nodes to the quad points, use:
		/// \code [cpp]
		/// Mapping mapping;
		/// mapping.setup(...);
		/// auto quad_values = mapping.get_shape_values() * nodal_values;
		/// \endcode
		/// where nodal_values is a matrix organized where the rows are the node indices and
		/// the columns are component of the quantity (1 if a scalar quantity).
		const simple_math::GeneralMatrix& get_shape_values() const { return current_shape_values_; }
		/// Returns a vector of matrices of the shape function derivatives with respect to the
		/// parametric coordinate system for the element and quadrature type indicated in the last
		/// setup() call. \see Mapping::setup
		/// The index of the vector indicates the quadrature point index.
		/// For each matrix, the first index (row) indicates the node index in the connectivity, while
		///  the second index (column) indicates the axis that the derivative is with respect to.
		const std::vector<simple_math::GeneralMatrix>& get_shape_derivatives() const { return current_shape_derivs_; }

	private:
		// The numeric limits of the underlying types will be used to create the array
		static constexpr long long MAX_QUADRATURE_TYPE_VALUE = static_cast<std::underlying_type<simple_math::QuadratureType>::type>(simple_math::QuadratureType::UNKNOWN);
		static constexpr long long MAX_INTEGRATION_ORDER = simple_math::MAX_GAUSS_ORDER;

		/// Stores the shape functions for at each quadrature point for a range of combinations of
		/// element type, quadrature type, and integration order.  Each entry stores a matrix that
		/// has each shape function value for each quadrature point.  The first index is the quadrature
		/// point index ordered by the Quadrature class and the second index is the shape function index
		/// ordered by the connectivity of the element type.
		/// For all_shape_values_, the first index is the element type, the second index is the quadrature
		/// type, and the third index is (integration_order-1).
		simple_math::GeneralMatrix all_shape_values_ [MAX_ELEMENT_TYPE_VALUE][MAX_QUADRATURE_TYPE_VALUE][MAX_INTEGRATION_ORDER];
		/// The shape function values for the current combination of element type, quadrature type, and
		/// integration order given to setup()
		simple_math::GeneralMatrix& current_shape_values_;

		/// Very similar to all_shape_values_.
		/// For entry, the first index (of the vector) is the matrix of the derivatives of the shape function
		/// values for a quadrature point.  The matrix is organized as the first index is the node/basis index
		/// and the second index is the parametric axis that the derivative is taken with respect to.
		std::vector<simple_math::GeneralMatrix> all_shape_derivs_[MAX_ELEMENT_TYPE_VALUE][MAX_QUADRATURE_TYPE_VALUE][MAX_INTEGRATION_ORDER];
		/// The shape function derivatives for the current combination of element type, quadrature type, and
		/// integration order given to setup()
		std::vector<simple_math::GeneralMatrix>& current_shape_derivs_;
	};
}