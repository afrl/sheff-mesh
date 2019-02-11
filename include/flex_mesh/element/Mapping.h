#pragma once

#include "simple_math/quadrature/Quadrature.h"
#include "utility/DynamicMatrix.h"

#include "ShapeFunctions.h"

#include <vector>

namespace flex_mesh
{
	class Mapping
	{
		void setup(const ElementType element_type, const simple_math::Quadrature& quadrature)
		{
			const auto element_type_index = static_cast<std::underlying_type<ElementType>::type>(element_type);
			const auto quad_type_index = static_cast<std::underlying_type<simple_math::QuadratureType>::type>(quadrature.get_type);
			const auto int_order_index = quadrature.get_integration_order() - 1;

			current_shape_values_ = all_shape_values_[element_type_index][quad_type_index][int_order_index];
			if (current_shape_values_.size() == 0)
			{
				using EigenRowType = Eigen::DenseBase<simple_math::GeneralMatrix>::RowXpr;
				for (auto qp = 0; qp < quadrature.n_quad_points(); qp++)
				{
					ShapeFunctions<Point, EigenRowType, simple_math::GeneralMatrix>::values(element_type, quadrature.parametric_point(qp), current_shape_values_.row(qp));
					// TODO fill derivatives
				}
			}
		}

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

	};
}