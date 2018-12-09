#pragma once

#include "flex_mesh/common/Connectivity.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/node/NodeUtilities.h"

namespace flex_mesh
{
	/// Note: to loop over the connectivity, use a ranged based-loop, such as:
	///       for (auto& node_ptr : edge)
	///         // do something with node_ptr
	template<typename NodeType>
	class Edge
	{
	public:
		typedef uint8_t ConnSize;

		/// Note: you can also use initializer lists here, such as:
		/// Edge<Node> e(EdgeType::E2N, { &n1, &n2 }); // n1 and n2 are Node objects
		Edge(const EdgeType& geom_type, const Connectivity<NodeType, ConnSize>& connectivity);
		virtual ~Edge() {};

		/// Get the type by const-ref
		const EdgeType& get_type() const;
		/// Returns an identifier that is unique to edges that have the same connectivity
		const std::size_t& get_conn_hash_id() const;
		/// Returns vector that is mean coordinate of the nodes defining this edge
		Point get_mean_coordinate() const;
		/// Returns the length of the edge
		FloatType get_length() const;

		/// Returns the node pointer at the beginning of the connectivity
		auto begin() const;
		/// Returns the node pointer at the end of the connectivity
		auto end() const;

		/// Swap a node in the connectivity for another and update the connectivity hash id
		void swap_nodes(const NodeType* old_node, NodeType* new_node);

		/// Compares based on mean coordinate of the two edges (takes into account tolerance)
		bool operator < (const Edge& other) const;
		/// Compares based on mean coordinate of the two edges (takes into account tolerance)
		bool operator > (const Edge& other) const;

	protected:
		/// The geometry/topolgoy type of edge (determines the number of nodes)
		EdgeType type_;
		/// The ordered list of nodes defining the edge.
		/// NOTE: this object does not own the nodes and will not delete them
		Connectivity<NodeType, ConnSize> connectivity_;
		/// An identifier that will match any other faces that have the exact same connectivity
		/// This is generated when the face is created
		std::size_t conn_hash_id_;

		/// Calculate the unique id for the connectivity based on the hash of the sorted connectivity
		/// pointers.  Any edge with the same nodes (in any order) will have the same hash.
		void update_conn_hash_id();
	};

	/////////////////
	// Definitions //
	/////////////////

	template<typename NodeType>
	inline Edge<NodeType>::Edge(const EdgeType & geom_type, const Connectivity<NodeType, ConnSize> & connectivity)
		: type_(geom_type), connectivity_(connectivity), conn_hash_id_(0)
	{
		assert(connectivity_.size() == get_num_nodes_for_type(type_));
		update_conn_hash_id();
	}

	template<typename NodeType>
	inline const EdgeType & Edge<NodeType>::get_type() const
	{
		return type_;
	}

	template<typename NodeType>
	inline const std::size_t & Edge<NodeType>::get_conn_hash_id() const
	{
		return conn_hash_id_;
	}

	template<typename NodeType>
	inline Point Edge<NodeType>::get_mean_coordinate() const
	{
		return get_mean_coordinate_of_nodes(connectivity_);
	}

	template<typename NodeType>
	inline FloatType Edge<NodeType>::get_length() const
	{
		if (type_ == EdgeType::E2N)
			return distance(*connectivity_[0], *connectivity_[1]);
		else if (type_ == EdgeType::E3N) // An approximation (should this be fixed?)
			return distance(*connectivity_[0], *connectivity_[2]);
		else
			throw _FM_MSG_EXCEPT("Unsupported edge type.");
	}

	template<typename NodeType>
	inline auto Edge<NodeType>::begin() const
	{
		return connectivity_.begin();
	}

	template<typename NodeType>
	inline auto Edge<NodeType>::end() const
	{
		return connectivity_.end();
	}

	template<typename NodeType>
	inline void Edge<NodeType>::swap_nodes(const NodeType* old_node, NodeType* new_node)
	{
		for (unsigned int i = 0; i < get_num_nodes_for_type(type_); ++i)
			if (connectivity_[i] == old_node)
				connectivity_[i] = new_node;
		update_conn_hash_id();
	}

	template<typename NodeType>
	inline bool Edge<NodeType>::operator<(const Edge & other) const
	{
		return compare_points(get_mean_coordinate(), other.get_mean_coordinate(), NodeType::get_compare_tolerance());
	}

	template<typename NodeType>
	inline bool Edge<NodeType>::operator>(const Edge & other) const
	{
		return compare_points(other.get_mean_coordinate(), get_mean_coordinate(), NodeType::get_compare_tolerance());
	}

	template<typename NodeType>
	inline void Edge<NodeType>::update_conn_hash_id()
	{
		std::vector<NodeType *> sorted_conn = connectivity_.copy_as_vector();
		std::sort(sorted_conn.begin(), sorted_conn.end());
		conn_hash_id_ = VectorPtrHasher<NodeType>::hash(sorted_conn);
	}

};