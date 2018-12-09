#pragma once

#include "Eigen/Eigen"

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

	typedef Eigen::Quaternion<FloatType> Quaternion;

	template<typename SizeType>
	class Orientations : public DynamicArray<Quaternion, SizeType>
	{
	public:
		using Parent = DynamicArray<Quaternion, SizeType>;

		Orientations() : Parent() {};
		template<typename QuaternionContainer>
		Orientations(const QuaternionContainer& quaternions) : Parent(quaternions) {};
		Orientations(const std::initializer_list<Quaternion>& quaternions) : Parent(quaternions) {};
		Orientations(const DynamicArray<Quaternion, SizeType >& arg) : Parent(arg) {};

		/// Assign the array to a new list of nodes
		template<typename OrientationContainer>
		Orientations& operator = (const OrientationContainer &nodes) {
			Parent::operator=(nodes);
			return *this;
		};
	};
}