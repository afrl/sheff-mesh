#pragma once

#include <array>
#include <forward_list>

#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/Utility.h"

#include "Eigen/Eigen"

namespace flex_mesh
{

	/// This class captures the geometrical/topological information of a node
	/// such as coordinates, numbering, parents, etc.
	class Node
	{

	public:
		//////////////////
		// Constructors //
		//////////////////

		/// Initializes every property to 0
		Node() noexcept;
		/// Creates a node with all variables initialized to 0 except the coordinates.
		Node(const Point& coord) noexcept;
		/// Creates a node with all variables initialized to 0 except the global number and coordinates.
		Node(const GlobalIndex global_num, const Point& coord) noexcept;
		/// Creates a node with all variables initialized, except the parent elements
		Node(const GlobalIndex global_num, const Point& coord, const PartitionIndex& parent_partition) noexcept;
		/// Copy constructor
		Node(const Node& other);
		/// Move constructor
		Node(Node&& other) noexcept = default;

		/////////////////
		// Get-setters //
		/////////////////

		const Point& get_coord() const;
		void set_coord(const Point& coord);
        const double& x() const;
        const double& y() const;
        const double& z() const;

		const LocalIndex& get_local_number() const;
		void set_local_number(const LocalIndex& local_num);

		const GlobalIndex& get_global_number() const;
		void set_global_number(const GlobalIndex& global_num);

		const PartitionIndex& get_parent_partition() const;
		void set_parent_partition(const PartitionIndex& parent_partition);

		/////////////////
		// Eigen casts //
		/////////////////

		/// Returns an Eigen::Map object that maps the coordinates to a 2D Eigen vector.
		/// It can be used as a normal Eigen vector, and it will modify the underlying
		/// coordinates if the vector is modified.
		Eigen::Map<Eigen::Vector2d> as_eigen_vec2d();
		/// Const version of as_eigen_vec2d().  Will not allow modification to the underlying
		/// coordinates of the node.
		Eigen::Map<const Eigen::Vector2d> as_eigen_vec2d() const;
		/// Returns an Eigen::Map object that maps the coordinates to a 3D Eigen vector.
		/// It can be used as a normal Eigen vector, and it will modify the underlying
		/// coordinates if the vector is modified.
		Eigen::Map<Eigen::Vector3d> as_eigen_vec3d();
		/// Const version of as_eigen_vec3d().  Will not allow modification to the underlying
		/// coordinates of the node.
		Eigen::Map<const Eigen::Vector3d> as_eigen_vec3d() const;

		///////////////
		// Operators //
		///////////////

		/// Tests if the addresses are the same
		bool operator == (const Node &other) const;
		/// Compares based on the coordinates of each node (first by x, then by y, then by z),
		/// but takes into account a tolerance.
		bool operator > (const Node &other) const;
		/// Compares based on the coordinates of each node (first by x, then by y, then by z),
		/// but takes into account a tolerance.
		bool operator < (const Node &other) const;
		/// Adds a vector to the coordinates of a node, effectively translating the node.
		/// Note: may invalidate existing KD trees, be sure to update them as needed
		Node& operator += (const Point &translation_vec);
		/// Assigns the coordinates of the node to a new set of coordinates.
		/// Note: may invalidate existing KD trees, be sure to update them as needed
		Node& operator = (const Point &coord);
		/// Const coordinate accessor.  Note: does not check if i is within 0-2.
		const FloatType& operator [] (unsigned int i) const;
		/// Non-const coordinate accessor.  Note: does not check if i is within 0-2.
		/// Note: may invalidate existing KD trees, be sure to update them as needed
		///       if the value of a coordinat eis modified.
		FloatType& operator [] (unsigned int i);

		///////////////////////////////
		// Static comparison methods //
		///////////////////////////////

		/// Compares nodes based on global number.  Useful for sorting by number
		static bool num_based_compare(const Node& node1, const Node& node2);

		/// Nodes owned by lower partition numbers are less than those owned by higher partition numbers,
		/// and within the same partition, compares local number.  Useful for sorting by partition the number.
		static bool part_based_compare(const Node& node1, const Node& node2);

		///////////////////////
		// Tolerance methods //
		///////////////////////

		/// Sets the tolerance used for > and < operators
		/// Note: Not thread-safe
		static void set_compare_tolerance(const double& tol);
		/// Sets the tolerance used for > and < operators within the node class
		static double get_compare_tolerance();

		////////////////////////////
		// Methods for mesh class //
		////////////////////////////

		/// Add a parent element to the vector of parents (which are internally stored as void*)
		template<typename ElementType>
		void add_parent_element(ElementType* elem);
		/// Remove a parent element to the vector of parents (which are internally stored as void*)
		template<typename ElementType>
		void remove_parent_element(ElementType* elem);
		/// Reset the partent list to contain nothing
		void reset_parents();
		/// Get the list of parent elements, though it must cast to the correct types
		std::vector<void*> const * const get_parent_elements() const;

	protected:
		/// The global number of the node, which must be an unique, 0-based index with all of the nodes
		/// forming a contiguous range from 0 to N-1, where N is the number of nodes.
		GlobalIndex global_num_;
		/// The local number of the node on the rank, which must be an unique, 0-based index with all
		/// of the nodes on the rank forming a contiguous range from 0 to N-1, where N is the number of nodes.
        LocalIndex local_num_;
		/// The coordinates of the node.  Assuemd to be 3D.
		Point coord_;
		/// It is assumed that only one partition will be responsible for each node.
		/// A mesh for a partition may have nodes belonging to other parititons depending
		/// element connectivity and ghost elements.  A FEA code using this library may or may
		/// not actually use this information.
		PartitionIndex parent_partition_;

		/// The mesh class will be responsible for using this property, but since the
		/// NodeTopology class cannot know what type of element will contain it, the parents
		/// will be stored as a vector of void pointers.  The mesh class will allocate and cast them.
		/// DO NOT mess with the pointers outside of the mesh class.
		/// AND NEVER delete these pointers.
		std::unique_ptr<std::vector<void*>> parent_elements_;
	};

	/////////////////
	// Definitions //
	/////////////////

	inline Node::Node() noexcept
		: global_num_(0), local_num_(0), coord_({ 0, 0, 0 }), parent_partition_(0)
	{
	}

	inline Node::Node(const Point& coord) noexcept
		: global_num_(0), local_num_(0), coord_(coord), parent_partition_(0)
	{
	}

	inline Node::Node(const GlobalIndex global_num, const Point & coord) noexcept
		: global_num_(global_num), local_num_(global_num), coord_(coord), parent_partition_(0)
	{
	}

	inline Node::Node(const GlobalIndex global_num, const Point & coord, const PartitionIndex & parent_partition) noexcept
		: global_num_(global_num), local_num_(global_num), coord_(coord), parent_partition_(parent_partition)
	{
	}

	inline Node::Node(const Node & other)
		: global_num_(other.global_num_), local_num_(other.local_num_), coord_(other.coord_), parent_partition_(other.parent_partition_)
	{
	}

	inline const Point& Node::get_coord() const
	{
		return coord_;
	}

	inline void Node::set_coord(const Point& coord)
	{
		coord_ = coord;
	}

    inline const double & Node::x() const
    {
        return coord_[0];
    }

    inline const double & Node::y() const
    {
        return coord_[1];
    }

    inline const double & Node::z() const
    {
        return coord_[2];
    }

	inline const LocalIndex & Node::get_local_number() const
	{
		return local_num_;
	}

	inline void Node::set_local_number(const LocalIndex & local_num)
	{
		local_num_ = local_num;
	}

	inline const GlobalIndex & Node::get_global_number() const
	{
		return global_num_;
	}

	inline void Node::set_global_number(const GlobalIndex & global_num)
	{
		global_num_ = global_num;
	}

	inline const PartitionIndex & Node::get_parent_partition() const
	{
		return parent_partition_;
	}

	inline void Node::set_parent_partition(const PartitionIndex & parent_partition)
	{
		parent_partition_ = parent_partition;
	}

	inline Eigen::Map<Eigen::Vector2d> Node::as_eigen_vec2d()
	{
		return Eigen::Map<Eigen::Vector2d>(&coord_[0], 2);
	}

	inline Eigen::Map<const Eigen::Vector2d> Node::as_eigen_vec2d() const
	{
		return Eigen::Map<const Eigen::Vector2d>(&coord_[0], 2);
	}

	inline Eigen::Map<Eigen::Vector3d> Node::as_eigen_vec3d()
	{
		return Eigen::Map<Eigen::Vector3d>(&coord_[0], 3);
	}

	inline Eigen::Map<const Eigen::Vector3d> Node::as_eigen_vec3d() const
	{
		return Eigen::Map<const Eigen::Vector3d>(&coord_[0], 3);
	}

	inline bool Node::operator==(const Node & other) const
	{
		return this == &other;
	}

	inline bool Node::operator>(const Node & other) const
	{
		return other < *this;
	}

	inline bool Node::operator<(const Node & other) const
	{
		return compare_points(this->coord_, other.coord_, get_compare_tolerance());
	}

	inline const FloatType & Node::operator[](unsigned int i) const
	{
		return coord_[i];
	}

	inline FloatType & Node::operator[](unsigned int i)
	{
		return coord_[i];
	}

	inline Node & Node::operator+=(const Point& translation_vec)
	{
		for (unsigned int i = 0; i < 3; i++)
			coord_[i] += translation_vec[i];
		return *this;
	}

	inline Node & Node::operator=(const Point& coord)
	{
		coord_ = coord;
		return *this;
	}

	inline bool Node::num_based_compare(const Node & node1, const Node & node2)
	{
		return node1.global_num_ < node2.global_num_;
	}

	inline bool Node::part_based_compare(const Node & node1, const Node & node2)
	{
		if (node1.parent_partition_ == node2.parent_partition_)
			return node1.global_num_ < node2.global_num_;
		else
			return node1.parent_partition_ < node2.parent_partition_;
	}

	inline void Node::set_compare_tolerance(const double & tol)
	{
		static double compare_tol_ = 1e-14;
		compare_tol_ = tol;
	}


	inline double Node::get_compare_tolerance()
	{
		static double compare_tol_;
		return compare_tol_;
	}

	inline void Node::reset_parents()
	{
		parent_elements_ = nullptr;
	}

	inline std::vector<void*> const * const Node::get_parent_elements() const
	{
		return parent_elements_.get();
	}

	template<typename ElementType>
	inline void Node::add_parent_element(ElementType * elem)
	{
		if (parent_elements_ == nullptr)
			parent_elements_ = std::make_unique<std::vector<void*>>();

		parent_elements_->push_back(static_cast<void*>(elem));
	}

	template<typename ElementType>
	inline void Node::remove_parent_element(ElementType * elem)
	{
		if (parent_elements_ == nullptr)
			return;

		for (auto i = 0; i < parent_elements_->size(); i++) {
			if ((*parent_elements_)[i] == elem) {
				(*parent_elements_)[i] = parent_elements_->back();
				parent_elements_->resize(parent_elements_->size() - 1);
			}
		}
	}
};
