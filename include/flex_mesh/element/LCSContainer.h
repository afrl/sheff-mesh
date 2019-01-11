#pragma once

#include "simple_math/EigenDefines.h"
#include "utility/DynamicArray.h"

#include "flex_mesh/common/Defines.h"

#include <assert.h>
#include <vector>

namespace flex_mesh
{
	/// An enum describing what level the local coordinate system will be
	/// assigned for a component.  This could be:
	/// - None: no LCS for the entire component
	/// - Component: one LCS for the entire component
	/// - Element: one LCS for each element in the component
	/// - Node: one LCS for each node in the component
	enum class LCSAssignment : char
	{
		None,
		Component,
		Element,
		Node
	};

	using Quaternion = simple_math::Quaternion;

	template<typename SizeType>
	class LCSContainer : public utility::DynamicArray<Quaternion, SizeType>
	{
	public:
		using Parent = utility::DynamicArray<Quaternion, SizeType>;

		LCSContainer() : Parent() {};
		template<typename QuaternionContainer>
		LCSContainer(const QuaternionContainer& quaternions) : Parent(quaternions) {};
		LCSContainer(const std::initializer_list<Quaternion>& quaternions) : Parent(quaternions) {};
		LCSContainer(const utility::DynamicArray<Quaternion, SizeType >& arg) : Parent(arg) {};

		/// Assign the array to a new list of nodes
		template<typename OrientationContainer>
		LCSContainer& operator = (const OrientationContainer &nodes) {
			Parent::operator=(nodes);
			return *this;
		};
	};
}