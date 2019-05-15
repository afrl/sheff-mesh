#pragma once

#include "simple_math/EigenDefines.h"

#include <cstdint>
#include <numeric>
#include <array>
#include <memory>

namespace flex_mesh
{
    using ConnIndex = uint8_t;
	using LocalIndex = std::uint32_t;
	using GlobalIndex = std::uint32_t;
	using PartitionIndex = std::uint16_t;
	using ComponentIndex = std::uint16_t;
	using FloatType = simple_math::FloatType;
	using Point = simple_math::Point;
	using Vector3D = simple_math::Vector3D;
	using NodeVector = simple_math::GeneralVector;
	using NodeX3Matrix = simple_math::GeneralMatrix;
};
