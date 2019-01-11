#pragma once

#define EIGEN_MAX_CPP_VER 11
#include "Eigen/Eigen"

#include <cstdint>
#include <numeric>
#include <array>
#include <memory>

namespace flex_mesh
{
	typedef std::uint32_t LocalIndex;
	typedef std::uint32_t GlobalIndex;
	typedef std::uint16_t PartitionIndex;
	typedef std::uint16_t ComponentIndex;
	typedef double FloatType;
	typedef std::array<FloatType, 3> Point;
	typedef Eigen::Matrix<FloatType, 3, 1> Vector3D;
	/// A vector with a fixed size that is larger than the connectivity expected for any element.
	typedef Eigen::Matrix<FloatType, Eigen::Dynamic, 1, 0, 120> ElementVector;
	/// A vector with a fixed size that is larger than the connectivity expected for any element
	/// and 3 columns for spatial derivatives.
	typedef Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic, 0, 120, 3> ElementMatrix3D;
};
