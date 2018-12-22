#pragma once

#include "Element.h"

namespace flex_mesh
{
	/// This class extends the Element type to include methods and properties
	/// to support distributed meshes.  Namely, there is the introduction of
	/// the parent partition and local number.  The base type is passed as a
	/// template parameter, allowing you mix and match of the different element
	/// classes.  For example, if you just wanted a distributed element without
	/// the extras that is helpful for mesh generation, just use:
	/// DistribElement<NodeType, Element<NodeType>>
	/// as your type.  However, if you want the mesh generation and distribute
	/// element functionality then simply use:
	/// MeshGenElement<NodeType, DistribElement<NodeType, Element<NodeType>>>
	/// as your type.  And so on.  It is highly recommended that you inherit
	/// from a concrete variation of the templates for a class you define for you
	/// application of this library.  The templates may be a verbose and a bit
	/// tricky, but it allows strong typing, with mix and match inheritance
	/// without the need for the woes of multiple inheritance.
	template<typename BaseType>
	class DistribElement : public BaseType
	{
	public:
		using NodeTypename = typename BaseType::NodeTypename;
		using NodePtr = NodeTypename * ;

		static_assert(
			std::is_base_of<BaseType, Element<NodeTypename>>::value,
			"BaseType must be a derived type of flex_mesh::Element<NodeType>"
			);

		//////////////////
		// Constructors //
		//////////////////

		DistribElement();
		DistribElement(const DistribElement&) = default;
		virtual ~DistribElement() {};

		/////////////////
		// Get-setters //
		/////////////////

		const LocalIndex& get_local_number() const;
		void set_local_number(const LocalIndex& local_num);

		const PartitionIndex& get_parent_partition() const;
		void set_parent_partition(const PartitionIndex& parent_partition);


	protected:
		/// The local number of the element is the number within a mesh for a particular partition.
		/// If the mesh contains all the partitions, then this will match the global number
		/// (assuming the the elements are ordered by global number).
		LocalIndex local_num_;
		/// It is assumed that only one partition will be responsible for each element.
		/// A mesh for a partition may have elements belonging to other parititons depending
		/// ghost elements.
		PartitionIndex parent_partition_;
	};

	/////////////////
	// Definitions //
	/////////////////

	template<typename BaseType>
	inline DistribElement<BaseType>::DistribElement()
		: BaseType(), local_num_(0), parent_partition_(0)
	{
	}

	template<typename BaseType>
	inline const LocalIndex & DistribElement<BaseType>::get_local_number() const
	{
		return local_num_;
	}

	template<typename BaseType>
	inline void DistribElement<BaseType>::set_local_number(const LocalIndex & local_num)
	{
		local_num_ = local_num;
	}

	template<typename BaseType>
	inline const PartitionIndex & DistribElement<BaseType>::get_parent_partition() const
	{
		return parent_partition_;
	}

	template<typename BaseType>
	inline void DistribElement<BaseType>::set_parent_partition(const PartitionIndex & parent_partition)
	{
		parent_partition_ = parent_partition;
	}
}