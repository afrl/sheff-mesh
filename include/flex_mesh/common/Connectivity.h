#pragma once

#include "utility/DynamicArray.h"

namespace flex_mesh
{
	/// A dynamic array node pointers used to describe  connectivity of edges, faces, and elements.
	/// Note: if the connectivity is pointers to constant nodes, use Connectivity<const Node, SizeType>
	/// Note: specify SizeType such that the type can hold the longest size expected.
	template<typename NodeType, typename SizeType>
	class Connectivity : public utility::DynamicArray<NodeType *, SizeType>
	{
	public:
		using Parent = utility::DynamicArray<NodeType *, SizeType>;

		Connectivity() : Parent() {};
		template<typename NodePtrContainer>
		Connectivity(const NodePtrContainer& nodes) : Parent(nodes) {};
		Connectivity(const std::initializer_list<NodeType *>& nodes) : Parent(nodes) {};
		Connectivity(const Parent& arg) : Parent(arg) {};

		/// Assign the array to a new list of nodes
		template<typename NodePtrContainer>
		Connectivity& operator = (const NodePtrContainer &nodes) {
			Parent::operator=(nodes);
			return *this;
		};
	};
};