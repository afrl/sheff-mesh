#pragma once

#include "Node.h"

namespace flex_mesh
{
	/// Add a local number and parent partition to the node class
	/// Pass the base node type as the template argument.
	/// This is tempalted to allow customizeable inheritances for types
	/// needed by the user of this library, similar to combination of Element,
	/// DistribElement, and MeshGenElement types.  For now, there is only one
	/// extension of the Node class, though a MeshGenNode class is expected.
	template<typename BaseType>
	class DistribNode : public BaseType
	{
		static_assert(
			std::is_base_of<BaseType, Node>::value,
			"BaseType must be a derived type of flex_mesh::Node"
		);

	public:
		using Parent = BaseType;

		// Inherit constructors from Node
		using BaseType::BaseType;
		DistribNode(const LocalIndex local_num, const GlobalIndex global_num, const Point& coords, const PartitionIndex parent_partition = 0);
		DistribNode(const DistribNode& other) = default;
		DistribNode(DistribNode&& other) = default;

		/////////////////
		// Get-setters //
		/////////////////

		const LocalIndex& get_local_number() const;
		void set_local_number(const LocalIndex& localNum);

		const PartitionIndex& get_parent_partition() const;
		void set_parent_partition(const PartitionIndex& partNum);

		///////////////////////////////
		// Static comparison methods //
		///////////////////////////////

		/// Nodes owned by lower partition numbers are less than those owned by higher partition numbers,
		/// and within the same partition, compares local number.  Useful for sorting by partition the number.
		static bool part_based_compare(const DistribNode& node1, const DistribNode& node2);

	protected:

		/// The local number of the node is the number within a mesh for a particular partition.
		/// If the mesh contains all the partitions, then this will match the global number
		/// (assuming the the nodes are ordered by global number).
		LocalIndex local_num_;
		/// It is assumed that only one partition will be responsible for each node.
		/// A mesh for a partition may have nodes belonging to other parititons depending
		/// element connectivity and ghost elements.  A FEA code using this library may or may
		/// not actually use this information.
		PartitionIndex parent_partition_;
	};

	template<typename BaseType>
	DistribNode<BaseType>::DistribNode(const LocalIndex local_num, const GlobalIndex global_num, const Point & coords, const PartitionIndex parent_partition)
		: BaseType(global_num, coords), local_num_(local_num), parent_partition_(parent_partition)
	{
	}

	template<typename BaseType>
	inline const LocalIndex & DistribNode<BaseType>::get_local_number() const
	{
		return local_num_;
	}

	template<typename BaseType>
	inline void DistribNode<BaseType>::set_local_number(const LocalIndex & localNum)
	{
		local_num_ = localNum;
	}

	template<typename BaseType>
	inline const PartitionIndex & DistribNode<BaseType>::get_parent_partition() const
	{
		return parent_partition_;
	}

	template<typename BaseType>
	inline void DistribNode<BaseType>::set_parent_partition(const PartitionIndex & partNum)
	{
		parent_partition_ = partNum;
	}

	template<typename BaseType>
	inline bool DistribNode<BaseType>::part_based_compare(const DistribNode & node1, const DistribNode & node2) {
		if (node1.parent_partition_ == node2.parent_partition_)
			return node1.local_num_ < node2.local_num_;
		else
			return node1.parent_partition_ < node2.parent_partition_;
	}
}