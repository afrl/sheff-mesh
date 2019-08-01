#pragma once

#include "flex_mesh/common/Defines.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/edge/Edge.h"
#include "flex_mesh/face/Face.h"

#include "LCSContainer.h"
#include "ShapeFunctions.h"

#include <set>

namespace flex_mesh
{
	template<typename NodeType>
	class Element
	{
	public:
		using NodeTypename = NodeType;
		using FaceTypename = Face<NodeType>;
		using EdgeTypename = Edge<NodeType>;
		using NodePtr = NodeType*;
		using FacePtr = Face<NodeType>*;
		using EdgePtr = Edge<NodeType>*;
		using ConnSize = ConnIndex;

		//////////////////
		// Constructors //
		//////////////////

		/// Note: component ID and global number are 0 by default
		///       use set_global_number() and set_component_id()
		Element();
		Element(const GlobalIndex global_num,
			const ElementType type,
			const std::vector<NodeType *> & connectivity,
			const LCSAssignment lcs_assignment = LCSAssignment::None) noexcept;
		Element(const Element&) = default;
		// Move constructor cannot be defaulted for MSVC, TODO implement
		//Element(const Element&&) noexcept = default;
		virtual ~Element() {};

		/// Initializes the element topology object
		/// Note: lcs are not allocated by default
		virtual void initialize(const ElementType type,
			const std::vector<NodeType *> & connectivity,
			const LCSAssignment lcs_assignment = LCSAssignment::None);

		/////////////////
		// Get-setters //
		/////////////////

		const GlobalIndex& get_global_number() const;
		void set_global_number(const GlobalIndex& global_num);

		const LocalIndex& get_local_number() const;
		void set_local_number(const LocalIndex& local_num);

		const ComponentIndex& get_component_id() const;
		void set_component_id(const ComponentIndex& component_id);

		const PartitionIndex& get_parent_partition() const;
		void set_parent_partition(const PartitionIndex& parent_partition);

		/// Get the type by const-ref
		/// Note: to set this, call initialize(...)
		const ElementType& get_type() const;
		/// Returns the number of nodes in the connectivity of the element.
		/// Note: to set this, call initialize(...)
		inline ConnSize get_num_nodes() const;

		/// Get the connectivity as pointers to the nodes
		const Connectivity<NodeType, ConnSize>& get_conn() const;
		/// Get the connectivity in terms of the local node numbering for a partition
		std::vector<LocalIndex> get_local_conn() const;

		/// Returns the node pointer at the beginning of the connectivity
		auto begin() const;
		/// Returns the node pointer at the end of the connectivity
		auto end() const;

		/////////////////
		// Comparisons //
		/////////////////

		/// Compares based on the mean coordinate
		bool operator < (const Element &Other) const;
		/// Compares based on the mean coordinate
		bool operator > (const Element &Other) const;

		///////////////////////
		// Interrogate Nodes //
		///////////////////////

		/// Returns vector that is mean coordinate of the nodes defining this face
		Point get_mean_coordinate() const;
		/// Fills a vector (which it will resize) with the global numbers of the nodes
		void get_global_node_numbers(std::vector<GlobalIndex>& numbers) const;
		/// Return a set of the nodes that are shared between this and another element
		std::set<NodeType *> get_shared_nodes(const Element &other) const;
        /// Returns a vector of pairs indicating the positions of any shared nodes in the
        /// connectivity of each element.  For each pair in the vector, there is a shared
        /// node.  For each pair, the first index is the position in this element's connectivity
        /// and the second index is the position in the other element's connectivity.
        /// Order of the vector is not predictable.
        std::vector<std::pair<ConnSize, ConnSize>> get_shared_node_positions(const Element &other) const;
		/// Returns if this element contains the given node
		bool contains_node(const NodeType& node) const;
		/// Given a map from old to new nodes, this will replace the pointers in
		/// the element and face connectivity containers.
		virtual void remap_element_nodes(const std::map<NodeType *, NodeType *> & old_to_new_nodes);

		//////////////////
		// Orientations //
		//////////////////

		/// Allocates the local coordinate systems
		void allocate_lcs(const LCSAssignment lcs_assignment);
		/// Rotate all lcs in the element by the same amount
		/// axes are defined in the global coordinate system
		void rotate_lcs(const Vector3D & axis,
			const double & angle);
		/// Rotate the lcs about one axis.
		void rotate_lcs(const short & axis_index,
			const double & angle);
		/// Set all lcs in the element to an identity matrix
		void reset_lcs();
		/// Directly get the local coordinate systems
		/// \attention Be careful.  The lcs may only contain one value per an element
		LCSContainer<ConnSize>& get_lcs();
		/// Get the local coordinate system for a specific node
		Quaternion get_lcs(const ConnSize& i);
		/// Returns the type of local coordinate system assignment for the element
		LCSAssignment get_lcs_assignment() const;

		///////////////////
		// Interpolation //
		///////////////////

		/// Evaluates the shape functions at the given parametric coordinate
		/// Resizes and fills the Eigen matrix with the shape function values
		virtual void get_shape_function_values(
			const Point& parametric_coordinate,
			NodeVector& shape_function_values) const;
		/// Evaluates the shape functions at the given parametric coordinate
		NodeVector get_shape_function_values(const Point& parametric_coordinate) const;

		/// Evaluates the derivatives of shape functions at the given parametric coordinate in
		/// terms of the parametric coordinates (dPhi/dxi)
		/// Resizes and fills the Eigen matrix with the shape function derivatives
		virtual void get_shape_function_derivatives(
			const Point& parametric_coordinate,
			NodeX3Matrix& shape_function_derivs) const;
		/// Evaluates the derivatives of shape functions at the given parametric coordinate in
		/// terms of the parametric coordinates (dPhi/dxi)
		NodeX3Matrix get_shape_function_derivatives(
			const Point& parametric_coordinate) const;

		/// Resizes and fills the Eigen matrix with the nodal coordinates
		/// Will be arranged with each row being a node
		void get_node_coordinates(NodeX3Matrix& coordinates) const;

        //////////
        // Misc //
        //////////

        bool contains_coordinate(const Point& x) const;

	protected:
		GlobalIndex global_num_;
		LocalIndex local_num_;
		ComponentIndex component_id_;
		ElementType type_;
		Connectivity<NodeType, ConnSize> connectivity_;
		LCSContainer<ConnSize> lcs_;
		LCSAssignment lcs_assignment_;
		/// It is assumed that only one partition will be responsible for each element.
		/// A mesh for a partition may have elements belonging to other parititons depending ghost elements.
		PartitionIndex parent_partition_;
	};

	/////////////////
	// Definitions //
	/////////////////

	template<typename NodeType>
	inline Element<NodeType>::Element()
		: global_num_(0), local_num_(0), component_id_(0), type_(ElementType::UNKNGEOM), lcs_assignment_(LCSAssignment::None), parent_partition_(0)
	{
	}

	template<typename NodeType>
	inline Element<NodeType>::Element(const GlobalIndex global_num, const ElementType type, const std::vector<NodeType*>& connectivity, const LCSAssignment lcs_assignment) noexcept
		: global_num_(global_num), local_num_(global_num), component_id_(0), parent_partition_(0)
	{
		initialize(type, connectivity, lcs_assignment);
	}

	template<typename NodeType>
	inline void flex_mesh::Element<NodeType>::initialize(const ElementType type, const std::vector<NodeType*>& connectivity, const LCSAssignment lcs_assignment)
	{
		type_ = type;
		connectivity_ = connectivity;
		allocate_lcs(lcs_assignment);
	}

	template<typename NodeType>
	inline const GlobalIndex & flex_mesh::Element<NodeType>::get_global_number() const
	{
		return global_num_;
	}

	template<typename NodeType>
	inline void Element<NodeType>::set_global_number(const GlobalIndex & global_num)
	{
		global_num_ = global_num;
	}

	template<typename NodeType>
	inline const LocalIndex & Element<NodeType>::get_local_number() const
	{
		return local_num_;
	}

	template<typename NodeType>
	inline void Element<NodeType>::set_local_number(const LocalIndex & local_num)
	{
		local_num_ = local_num;
	}

	template<typename NodeType>
	inline const ComponentIndex & Element<NodeType>::get_component_id() const
	{
		return component_id_;
	}

	template<typename NodeType>
	inline void Element<NodeType>::set_component_id(const ComponentIndex & component_id)
	{
		component_id_ = component_id;
	}

	template<typename NodeType>
	inline const PartitionIndex & Element<NodeType>::get_parent_partition() const
	{
		return parent_partition_;
	}

	template<typename NodeType>
	inline void Element<NodeType>::set_parent_partition(const PartitionIndex & parent_partition)
	{
		parent_partition_ = parent_partition;
	}

	template<typename NodeType>
	inline const ElementType & Element<NodeType>::get_type() const
	{
		return type_;
	}

	template<typename NodeType>
	inline typename Element<NodeType>::ConnSize Element<NodeType>::get_num_nodes() const
	{
		return connectivity_.size();
	}

	template<typename NodeType>
	inline const Connectivity<NodeType, typename Element<NodeType>::ConnSize>& Element<NodeType>::get_conn() const
	{
		return connectivity_;
	}

	template<typename NodeType>
	inline std::vector<LocalIndex> Element<NodeType>::get_local_conn() const
	{
		std::vector<LocalIndex> local_indices;
		for (auto& n_ptr : *this)
			local_indices.push_back(n_ptr->get_local_number());
		return local_indices;
	}

	template<typename NodeType>
	inline auto Element<NodeType>::begin() const
	{
		return connectivity_.begin();
	}

	template<typename NodeType>
	inline auto Element<NodeType>::end() const
	{
		return connectivity_.end();
	}

	template<typename NodeType>
	inline bool Element<NodeType>::operator<(const Element & other) const
	{
		return compare_points(get_mean_coordinate(), other.get_mean_coordinate(), NodeType::get_compare_tolerance());
	}

	template<typename NodeType>
	inline bool Element<NodeType>::operator>(const Element & other) const
	{
		return compare_points(other.get_mean_coordinate(), get_mean_coordinate(), NodeType::get_compare_tolerance());
	}

	template<typename NodeType>
	inline Point Element<NodeType>::get_mean_coordinate() const
	{
		return get_mean_coordinate_of_nodes(connectivity_);
	}

	template<typename NodeType>
	inline void Element<NodeType>::get_global_node_numbers(std::vector<GlobalIndex>& numbers) const
	{
		numbers.resize(connectivity_.size());
		for (auto i = 0; i < connectivity_.size(); i++) {
			numbers[i] = connectivity_[i]->get_global_number();
		}
	}

	template<typename NodeType>
	inline std::set<NodeType*> Element<NodeType>::get_shared_nodes(const Element & other) const
	{
		std::set<NodeType*> common_nodes(connectivity_.begin(), connectivity_.end());
		for (auto& n_ptr : other) {
			auto find_iter = common_nodes.find(n_ptr);
			if (find_iter == common_nodes.end())
				common_nodes.remove(find_iter);
		}
		return common_nodes;
	}

    template<typename NodeType>
    inline std::vector<std::pair<typename Element<NodeType>::ConnSize, typename Element<NodeType>::ConnSize>> Element<NodeType>::get_shared_node_positions(const Element & other) const
    {
        std::vector<std::pair<ConnSize, ConnSize>> shared_node_pos;
        const auto & other_conn = other.get_conn();
        for (ConnSize other_pos = 0; other_pos < other_conn.size(); other_pos++)
        {
            for (ConnSize this_pos = 0; this_pos < connectivity_.size(); this_pos++)
            {
                if (other_conn[other_pos] == connectivity_[this_pos])
                {
                    shared_node_pos.push_back({ this_pos, other_pos });
                }
            }
        }
        return shared_node_pos;
    }

    template<typename NodeType>
	inline bool Element<NodeType>::contains_node(const NodeType & node) const
	{
		for (auto& n_ptr : *this)
			if (n_ptr == &node)
				return true;
		return false;
	}

	template<typename NodeType>
	inline void Element<NodeType>::remap_element_nodes(const std::map<NodeType*, NodeType*>& old_to_new_nodes)
	{
		for (auto i = 0; i < connectivity_.size(); i++) {
			auto find_iter = old_to_new_nodes.find(connectivity_[i]);
			if (find_iter != old_to_new_nodes.end())
				connectivity_[i] = find_iter->second;
		}
	}

	template<typename NodeType>
	inline void Element<NodeType>::allocate_lcs(const LCSAssignment lcs_assignment)
	{
		lcs_assignment_ = lcs_assignment;
		if (lcs_assignment == LCSAssignment::None)
			lcs_.clear();
		else if (lcs_assignment == LCSAssignment::Node) {
			lcs_ = std::vector<Quaternion>(connectivity_.size());
			for (auto& lcs : lcs_)
				lcs.setIdentity();
		}
		// There will only be one lcs per an element for LCSAssignment::Element
		// and LCSAssignment::Component
		else {
			lcs_ = std::vector<Quaternion>(1);
			for (auto& lcs : lcs_)
				lcs.setIdentity();
		}
	}

	template<typename NodeType>
	inline void Element<NodeType>::rotate_lcs(const Vector3D & axis, const double & angle)
	{
		Quaternion rotation(Eigen::AngleAxis<FloatType>(angle, axis));
		for (auto& orient : lcs_)
			orient *= rotation;
	}

	template<typename NodeType>
	inline void Element<NodeType>::rotate_lcs(const short & axis_index, const double & angle)
	{
		Vector3D axis; axis.setZero();
		axis(axis_index, 1.0);
		rotate_lcs(axis, angle);
	}

	template<typename NodeType>
	inline void Element<NodeType>::reset_lcs()
	{
		for (auto& orient : lcs_)
			orient.setIdentity();
	}

	template<typename NodeType>
	inline LCSContainer<typename Element<NodeType>::ConnSize> & Element<NodeType>::get_lcs()
	{
		return lcs_;
	}

	template<typename NodeType>
	inline Quaternion Element<NodeType>::get_lcs(const ConnSize & i)
	{
		if (lcs_.size() == 0)
			return Quaternion::Identity();
		// Handle the case if there is only one lcs for the whole element
		else if (lcs_.size() == 1)
			return lcs_[0];
		return lcs_[i];
	}

	template<typename NodeType>
	inline LCSAssignment Element<NodeType>::get_lcs_assignment() const
	{
		return lcs_assignment_;
	}

	template<typename NodeType>
	inline void flex_mesh::Element<NodeType>::get_shape_function_values(const Point & parametric_coordinate, NodeVector & shape_function_values) const
	{
		ShapeFunctions<Point, NodeVector, NodeX3Matrix>::values(type_, parametric_coordinate, shape_function_values);
	}

	template<typename NodeType>
	inline NodeVector Element<NodeType>::get_shape_function_values(const Point & parametric_coordinate) const
	{
		NodeVector shape_function_values;
		get_shape_function_values(parametric_coordinate, shape_function_values);
		return shape_function_values;
	}

	template<typename NodeType>
	inline void flex_mesh::Element<NodeType>::get_shape_function_derivatives(const Point & parametric_coordinate, NodeX3Matrix & shape_function_derivs) const
	{
		ShapeFunctions<Point, NodeVector, NodeX3Matrix>::derivatives(type_, parametric_coordinate, shape_function_derivs);
	}

	template<typename NodeType>
	inline NodeX3Matrix Element<NodeType>::get_shape_function_derivatives(const Point & parametric_coordinate) const
	{
		NodeX3Matrix derivatives;
		get_shape_function_derivatives(derivatives);
		return derivatives;
	}

	template<typename NodeType>
	inline void Element<NodeType>::get_node_coordinates(NodeX3Matrix & coordinates) const
	{
		coordinates.resize(connectivity_.size(), 3);
		for (ConnSize i = 0; i < connectivity_.size(); i++)
			for (auto j = 0; j < 3; j++)
				coordinates(i, j) = (*connectivity_[i])[j];
	}

    template<typename NodeType>
    inline bool Element<NodeType>::contains_coordinate(const Point& x) const
    {
        auto face_infos = get_face_info(type_);
        for (auto& face_info : face_infos)
        {
            Connectivity<NodeType, ConnSize> face_connectivity;
            for (auto node_index : face_info.connectivity)
            {
                face_connectivity.push_back(connectivity_[node_index]);
            }
            Face<NodeType> face(face_info.type, face_connectivity);
            if (face.calc_signed_norm_projection(x) > 0)
            {
                return false;
            }
        }
        return true;
    }

};
