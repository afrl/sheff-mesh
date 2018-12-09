#pragma once


#include "flex_mesh/common/Utility.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace flex_mesh
{

	template<typename NodeType>
	struct Vertex
	{
		Connectivity<NodeType, uint8_t> neighbors;

		/// Returns the neighbor at a given index, but this method allows indices to be negative
		/// or larger than the number of neighbors to more easily facilitate algorithms with neighbors.
		/// For example, -1 will return the last node and the neighbors.size() index will return
		/// the first node.  The valid range of indices are [-neighbors.size()+1, 2*neighbors.size()-1].
		/// In other words, if there are 3 neighbors, then the valid range is [-2, 5]
		NodeType* get_neighbor(const int16_t i) const {
			assert(i >= -neighbors.size() + 1);
			assert(i <= 2 * neighbors.size() - 1);
			if (i < 0)
				return neighbors[i + neighbors.size()];
			else if (i >= neighbors.size())
				return neighbors[i - neighbors.size()];
			else
				return neighbors[static_cast<uint8_t>(i)];
		}
	};

	template<typename NodeType>
	struct ShapeFunctionHelper
	{
		/// A helper function for computing shape functions, see Eqn. 7 in
		/// http://graphics.stanford.edu/~wicke/publications/FEM/wbg-afmocp-07.pdf
		double compute_c(const NodeType& a, const NodeType& b, const Vector3D& x);
		/// A helper function for computing shape functions, returns the volume a tetrahedral
		/// defined by three nodes and a point.
		double compute_v(const NodeType& a, const NodeType& b, const NodeType& c, const Vector3D& x);

		/// Refer to Eqn. 22 in Appendix A of:
		/// http://graphics.stanford.edu/~wicke/publications/FEM/wbg-afmocp-07.pdf
		double compute_grad_c(const NodeType& a, const NodeType& b, const Vector3D& x);
		/// Refer to Eqn. 21 in Appendix A of:
		/// http://graphics.stanford.edu/~wicke/publications/FEM/wbg-afmocp-07.pdf
		double compute_grad_v(const NodeType& a, const NodeType& b, const NodeType& c, const Vector3D& x);

		/// Evaluates x / sin(x) using a series expansion with error O(x^6),
		/// avoids numerical problems near x = 0, which should result in 1.
		/// Refer to https://www.wolframalpha.com/input/?i=x%2Fsin(x)
		inline double eval_x_over_sinx(const double x);

		/// Approximates arccosine with an absolute error <= 6.7e-5 radians
		/// This call requires 20-30% less time to evaluate than the built in acos
		static double approx_acos(double x);
		/// Approximates arccosine with an absolute error <= 1.7e-3 radians
		/// This call requires 50-70% less time to evaluate than the built in acos
		static double approx_acos2(double x);

		// Consider hashing combinations move compute methods to a class
		// and track combinations internally
		std::map<std::size_t, double> c_hash_map;
		std::map<std::size_t, double> V_hash_map;
	};


	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::compute_c(const NodeType & a, const NodeType & b, const Vector3D & x)
	{
		//const auto hash = VectorPtrHasher<const NodeType>::hash({ &a, &b });
		//const auto find_iter = this->c_hash_map.find(hash);
		//if (find_iter != this->c_hash_map.end())
		//	return find_iter->second;

		const auto xa = a.as_eigen_vec3d();
		const auto xb = b.as_eigen_vec3d();
		const auto xa_x = xa - x;
		const auto xb_x = xb - x;

		const double c1 = (xa_x.cross(xb_x)).norm() / 6.0;
		const double c2_inner = (xa_x).dot(xb_x) / (xa_x.norm() * xb_x.norm());
		const double c2 = approx_acos(c2_inner);
		const double c = c1 * c2;

		//this->c_hash_map[hash] = c;

		return c;
	}

	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::compute_v(const NodeType & a, const NodeType & b, const NodeType & c, const Vector3D & x)
	{
		const auto xa = a.as_eigen_vec3d();
		const auto xb = b.as_eigen_vec3d();
		const auto xc = c.as_eigen_vec3d();

		return fabs(((xa - x).dot((xb - x).cross(xc - x)))) / 6.0;
	}

	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::compute_grad_c(const NodeType & a, const NodeType & b, const Vector3D & x)
	{
		const auto xa = a.as_eigen_vec3d();
		const auto xb = b.as_eigen_vec3d();
		const auto da = x - xa;
		const auto db = x - xb;
		const auto da_hat = da.normalized();
		const auto db_hat = db.normalized();

		const double grad_c = (da_hat.dot(db_hat) * (da_hat * db.norm() + db_hat * da.norm())
			- da - db + eval_x_over_sinx(approx_acos2(da_hat.dot(db_hat))) * ((da_hat.cross(db_hat).cross(xb - xa)))) / 6.0;

		return grad_c;
	}

	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::compute_grad_v(const NodeType & a, const NodeType & b, const NodeType & c, const Vector3D & x)
	{
		const auto xa = a.as_eigen_vec3d();
		const auto xb = b.as_eigen_vec3d();
		const auto xc = c.as_eigen_vec3d();

		return ((xc - xa).cross(xb - xa)).norm() / 6.0;
	}

	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::eval_x_over_sinx(const double x)
	{
		return 1.0 + x * x / 6.0 + 7.0 * x *x *x * x / 360.0;
	}

	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::approx_acos(double x)
	{
		const double negate = double(x < 0);
		x = abs(x);
		double ret = -0.0187293;
		ret = ret * x;
		ret = ret + 0.0742610;
		ret = ret * x;
		ret = ret - 0.2121144;
		ret = ret * x;
		ret = ret + 1.5707288;
		ret = ret * sqrt(1.0 - x);
		ret = ret - 2 * negate * ret;
		return negate * 3.14159265358979 + ret;
	}

	template<typename NodeType>
	inline double ShapeFunctionHelper<NodeType>::approx_acos2(double x)
	{
		const double a = -0.939115566365855;
		const double b = 0.9217841528914573;
		const double c = -1.2845906244690837;
		const double d = 0.295624144969963174;
		return M_PI / 2 + (a*x + b * x*x*x) / (1 + c * x*x + d * x*x*x*x);
	}

}