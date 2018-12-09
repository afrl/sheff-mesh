#pragma once

#include "DistribElement.h"

namespace flex_mesh
{
	/// This class extends the Element type to include methods and properties
	/// to support mesh generation and other complex algorithms.
	/// Namely, there is the introduction of faces, edges, and a hash identifier.
	/// The has identifier will be unique for any combination of a connectivity
	/// involving the same nodes, which is helpful for finding matching elements,
	/// faces, and edges.
	/// The base type is passed as a template parameter, allowing you mix and
	/// match of the different element classes.  For example, if you just wanted
	/// a mesh gen element without the extras provided from DistribElement just use:
	/// MeshGenElement<NodeType, Element<NodeType>>
	/// as your type.  However, if you want the mesh generation and distribute
	/// element functionality then simply use:
	/// MeshGenElement<NodeType, DistribElement<NodeType, Element<NodeType>>>
	/// as your type.  And so on.  It is highly recommended that you inherit
	/// from a concrete variation of the templates for a class you define for you
	/// application of this library.  The templates may be a verbose and a bit
	/// tricky, but it allows strong typing, with mix and match inheritance
	/// without the need for the woes of multiple inheritance.
	template<typename NodeType, typename BaseType>
	class MeshGenElement : public BaseType
	{
		static_assert(
			std::is_base_of<Element<NodeType>, BaseType>::value,
			"BaseType must be a derived type of flex_mesh::Element<NodeType>"
			);

	public:

		using Parent = BaseType;

		//////////////////
		// Constructors //
		//////////////////

		using Parent::Parent;

		/// Initializes the element topology object
		/// Note: orientations are not allocated by default
		void initialize(const ElementType type,
			const std::vector<NodeType *> & connectivity,
			const LCSAssignment lcs_assignment = LCSAssignment::None);

		///////////////
		// Accessors //
		///////////////

		/// Returns an identifier that is unique to elements that have the same connectivity
		const std::size_t& get_conn_hash_id() const;
		/// Returns vector that is mean coordinate of the nodes defining this face
		Point get_mean_coordinate() const;
		/// Return true if the element is linear, false otherwise
		bool is_linear() const;
		/// Return a vector of pointers for nodes on the element's corners (in order)
		std::vector<NodeType*> get_corner_nodes() const;

		//////////////////
		// Manipulators //
		//////////////////

		/// Given a map from old to new nodes, this will replace the pointers in
		/// the element and face connectivity containers.
		virtual void remap_element_nodes(const std::map<NodeType *, NodeType *> & old_to_new_nodes) override;
		/// Eliminates midnodes from the element, making it linear.  Also eleminates
		/// extra orientation data.
		void make_element_linear();
		/// Swap one node for another node in this element.
		/// Note: invalidates parent element map for the nodes in the mesh, if it was being tracked
		/// call update_parents via the mesh class.
		void swap_nodes(const NodeType* old_node, NodeType* new_node);

		///////////////////////
		// Face/Edge methods //
		///////////////////////

		/// Returns the faces of this element.
		const std::vector<std::unique_ptr<FaceType>>& get_faces();
		/// Return true if this element shares a face with otherElement
		bool shares_face(const MeshGenElement &otherElement) const;
		/// Returns vector with faces that are shared between this and another element
		std::vector<Face<NodeType> *> get_adjacent_faces(MeshGenElement &other);

		/// Returns a set of edges (this element owns them and they will cleaned up
		/// by the element with reset_faces_and_edges, initialize, or the destructor
		/// is called.
		const std::vector<std::unique_ptr<EdgeType>>& get_edges();
		/// Like get_edges but the caller gets ownership.  They will get cleaned up
		/// automatically when the returned vector falls out of scope (if no new shared
		/// pointers are created).
		std::vector<std::shared_ptr<EdgeType>> get_edges_shared();

		/// Deletes any face and edge pointers, but new ones will be created whenever
		/// faces or edges are asked for again.
		void reset_faces_and_edges();

		//////////////////
		// Orientations //
		//////////////////

		/// Mirror element orientations
		void mirror_orientations(const short& axis_index);

		/// Stretch element orientations
		void stretch_orientations(const short& axis_index,
			const double& factor);

	protected:
		/// An identifier that will match any other faces that have the exact same connectivity
		/// This is generated when the face is created
		std::size_t conn_hash_id_;
		/// The faces of the element.  Only generated when get_faces is called.
		/// It is more memory efficient to have a pointer than an empty vector when
		/// faces are not created.
		std::unique_ptr<std::vector<std::unique_ptr<FaceType>>> faces_;
		/// If the edges are asked for once, then the element will save them for access in the future.
		/// Note: they are not created until get_edges is called, so even for methods in this class,
		/// use the method to access them.
		/// It is more memory efficient to have a pointer than an empty vector when
		/// faces are not created.
		std::unique_ptr<std::vector<std::unique_ptr<EdgeType>>> edges_;

		/// Calculate the unique id for the connectivity based on the hash of the sorted connectivity
		/// pointers.  Any edge with the same nodes (in any order) will have the same hash.
		void update_conn_hash_id();
	};

	/////////////////
	// Definitions //
	/////////////////

	template<typename NodeType, typename BaseType>
	inline void MeshGenElement<NodeType, BaseType>::initialize(const ElementType type, const std::vector<NodeType*>& connectivity, const LCSAssignment lcs_assignment)
	{
		Parent::initialize(type, connectivity, lcs_assignment);
		update_conn_hash_id();
	}
	template<typename NodeType, typename BaseType>
	inline const std::size_t & MeshGenElement<NodeType, BaseType>::get_conn_hash_id() const
	{
		return conn_hash_id_;
	}

	template<typename NodeType, typename BaseType>
	inline Point MeshGenElement<NodeType, BaseType>::get_mean_coordinate() const
	{
		return get_mean_coordinate_of_nodes(Parent::connectivity_);
	}

	template<typename NodeType, typename BaseType>
	inline bool MeshGenElement<NodeType, BaseType>::is_linear() const
	{
		return is_type_linear(Parent::geom_type_);
	}
	template<typename NodeType, typename BaseType>
	inline std::vector<NodeType*> MeshGenElement<NodeType, BaseType>::get_corner_nodes() const
	{
		std::vector<NodeType*> corner_nodes;
		throw _FM_MSG_EXCEPT_C("Not implemented.");
		return corner_nodes;
	}

	template<typename NodeType, typename BaseType>
	inline void MeshGenElement<NodeType, BaseType>::remap_element_nodes(const std::map<NodeType*, NodeType*>& old_to_new_nodes)
	{
		Parent::remap_element_nodes(old_to_new_nodes);
		update_conn_hash_id();
		reset_faces_and_edges();
	}

	template<typename NodeType, typename BaseType>
	inline void MeshGenElement<NodeType, BaseType>::make_element_linear()
	{
		throw _FM_MSG_EXCEPT_C("Not implemented.");
	}

	template<typename NodeType, typename BaseType>
	inline void MeshGenElement<NodeType, BaseType>::swap_nodes(const NodeType * old_node, NodeType * new_node)
	{
		for (unsigned int i = 0; i < get_num_nodes_for_type(Parent::type_); ++i)
			if (Parent::connectivity_[i] == old_node)
				Parent::connectivity_[i] = new_node;
		update_conn_hash_id();
		reset_faces_and_edges();
	}

	template<typename NodeType, typename BaseType>
	inline void MeshGenElement<NodeType, BaseType>::update_conn_hash_id()
	{
		// Create connectivity hash
		std::vector<NodeType *> sorted_conn = Parent::connectivity_.copy_as_vector();
		std::sort(sorted_conn.begin(), sorted_conn.end());
		conn_hash_id_ = VectorPtrHasher<NodeType>::hash(sorted_conn);
	}

	template<typename NodeType, typename BaseType>
	inline void MeshGenElement<NodeType, BaseType>::reset_faces_and_edges()
	{
		faces_->clear();
		faces_ = nullptr;
		edges_->clear();
		edges_ = nullptr;
	}
}