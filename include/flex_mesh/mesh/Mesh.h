#pragma once

#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/Exceptions.h"
#include "flex_mesh/common/Traits.h"

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace flex_mesh
{
// Tests to see if elements contain distributed information
#define HAS_LOCAL_NUMBERING_DECLARE(T) typename T_ = T, typename std::enable_if<flex_mesh_traits::has_get_local_number<T_>::value>::type* = nullptr
#define HAS_LOCAL_NUMBERING_DEFINE(T) typename T_, typename std::enable_if<flex_mesh_traits::has_get_local_number<T_>::value>::type*
#define MISSING_LOCAL_NUMBERING_DECLARE(T) typename T_ = nullptr, typename std::enable_if<!flex_mesh_traits::has_get_local_number<T_>::value>::type* = nullptr
#define MISSING_LOCAL_NUMBERING_DEFINE(T) typename T_, typename std::enable_if<!flex_mesh_traits::has_get_local_number<T_>::value>::type*

	template<typename T>
	class MeshEntityContainer
	{
	public:
		typedef typename std::vector<std::unique_ptr<T>>::iterator Iterator;
		typedef typename std::vector<std::unique_ptr<T>>::const_iterator ConstIterator;

		MeshEntityContainer() {};
		~MeshEntityContainer() {};

		/// Returns the size of the container
		size_t size() const { return objects_.size(); }
		/// Returns the begin iterator of the container
		/// Note: accessing this non-const iterator will mark the container as changed,
		/// which will result in updates for other data structures (such as any KD trees,
		/// etc.) when a data structure that depends on this container is requested. Use
		/// the const version whenever possible.
		auto begin() {
			mark_as_changed();
			return objects_.begin();
		}
		/// Returns the begin iterator of the container
		auto begin() const { return objects_.begin(); }
		/// Returns the begin iterator of the container (forcing const)
		auto begin_const() const { return objects_.begin(); }
		/// Returns the end iterator of the container
		auto end() const { return objects_.end(); }
		/// Returns a raw pointer to the first object added
		T * const front() const { return  !objects_.empty() ? objects_.front().get() : nullptr; }
		/// Returns a raw pointer to the last object added
		T * const back() const { return !objects_.empty() ? objects_.back().get() : nullptr; }
		/// Add an object to the container
		/// Returns a raw pointer to that object
#ifndef SWIG
		T * add(std::unique_ptr<T> obj) {
			mark_as_changed();
			objects_.push_back(std::move(obj));
			return back();
		};
#endif
		/// Add an object to the container and have it call its default constructor
		/// Returns a raw pointer to that object
		T * add() {
			mark_as_changed();
			objects_.push_back(std::make_unique<T>());
			return back();
		};
		/// Remove an object from the container given an iterator
		/// Returns the next iterator
		Iterator remove(Iterator it) {
			mark_as_changed();
			*it = std::move(objects_.back());
			objects_.resize(objects_.size() - 1);
			return it++;
		};
		/// Remove a list of objects from the container given iterators
		void remove(std::list<Iterator> iter_list) {
			mark_as_changed();
			for (auto it : iter_list)
				this->remove(it);
		}
		/// Given a set of object pointers, this will remove all of them from the container.
		/// Note: if the container stores elements, it does not update element connectivities,
		/// that must be done in the algorithm that call this.
		void remove_nodes(std::set<T *> objects_to_remove) {
			mark_as_changed();
			for (auto it = begin(); nIt != end(); ) {
				auto setIt = objects_to_remove.find(it->get());
				if (setIt != objects_to_remove.end()) {
					objects_to_remove.erase(setIt);
					it = remove_node(it);
				}
				else
					it++;
			}
		};

		template<typename Node_, typename Element_>
		friend class Mesh;

	private:
		std::vector<std::unique_ptr<T>> objects_;
		bool changed_flag_for_kd_tree_;

		/// This is called whenever a method is invoke that could change the data structure
		/// or the contained elements.
		void mark_as_changed() {
			changed_flag_for_kd_tree_ = true;
		}
		/// Called by KDTree class
		void reset_kd_tree_changed_flag() {
			changed_flag_for_kd_tree_ = false;
		}
	};

	/// A base container for NodeTopology and ElementTopology objects
	/// that compose a mesh. This is meant to be very simple and only
	/// contain methods related to accessing and manipulating a mesh.
	/// All algorithms will be found in other classes that derive
	template<typename Node_, typename Element_>
	class Mesh
	{
	public:
		using NodeTypename = Node_;
		using ElementTypename = Element_;

		//////////////////////////
		// Constructors and     //
		// Mesh-Mesh operations //
		//////////////////////////

		Mesh() = default;
		/// A deepcopy constructor
		Mesh(const Mesh& other);
		/// Move constructor
		Mesh(Mesh&& other) noexcept;
		/// Add another to this mesh, making a deepcopy of nodes and elements
		void add_to_self(const Mesh & other);
		/// Move another mesh into this mesh
		/// \important Note: the other mesh will be empty after this
		void move_to_self(Mesh & other);
		/// Add another to this mesh, making a deepcopy of nodes and elements
		Mesh& operator+=(const Mesh& other) noexcept;
		/// Destroys the current mesh and copies the other mesh.
		Mesh& operator=(const Mesh& other) noexcept;
		/// Destroys the current mesh and moves the other mesh to this one.
		Mesh& operator=(Mesh&&) noexcept = default;

		/// Clears the mesh object of its contents
		void clear();

		/////////////////////
		// Node operations //
		/////////////////////

		/// Access element container by non-const reference
		const MeshEntityContainer<Node_>& nodes() const;
		/// Access element container by const reference
		MeshEntityContainer<Node_>& nodes();
        /// Add a node via the coordinate
		Node_ * add_node_by_coord(const Point& coord);
		/// Set the tolerance used for comparisons between node coordinates
		/// Note: tolerance is common to all Mesh instances if they are templated
		/// for the same Node_ type.
		void set_compare_tolerance(double tol);
		/// Get the tolerance used for comparisons between node coordinates
		/// Note: tolerance is common to all Mesh instances if they are templated
		/// for the same Node_ type.
		double get_compare_tolerance() const;

		////////////////////////
		// Element operations //
		////////////////////////

		/// Access element container by non-const reference
		const MeshEntityContainer<Element_>& elements() const;
		/// Access element container by const reference
		MeshEntityContainer<Element_>& elements();
        /// Reset element face/edges
		void reset_element_faces_and_edges();
		/// Returns the set of component IDs for all the elements
		std::set<ComponentIndex> get_component_ids() const;
		/// Returns a vector of element pointers that belong to the component
		/// Note: only includes elements in this mesh object
		std::vector<Element_ *> get_component(ComponentIndex component_id) const;

		//////////////////////////////
		// Distributed mesh methods //
		//////////////////////////////

		/// Returns the number of partitions (assuming they are contiguous 0-based when assigned)
		template<HAS_LOCAL_NUMBERING_DECLARE(Element_)>
		unsigned int get_num_partitions() const;

		/// This method will go through the nodes and elements and update
		/// the global numbering appropriately.  This should be called after
		/// an operation of adding nodes or elements to the mesh.  It is not
		/// called automatically for performance reasons.
		/// This will also set the local numbering to be the same as the global,
		/// which is neccessary for non-distributed I/O.
		template<HAS_LOCAL_NUMBERING_DECLARE(Element_)>
		void update_global_numbering();
		template<MISSING_LOCAL_NUMBERING_DECLARE(Element_)>
		void update_global_numbering();

		/// Renumbers the nodes and element local numbers for the current partition.
		/// for node not : the current parition, the local number is unchanged.
		/// Therefore, never rely on the local number to see if it belons to this partition;
		/// instead, use get_parent_partition() for the elements and nodes to test.
		/// Returns the number of nodes and elements respectively for the current partition
		/// Template note: template siables this method if the Element_ type does not contain a
		/// method specifically called get_local_number()
		template<HAS_LOCAL_NUMBERING_DECLARE(Element_)>
		std::pair<std::set<Node_ *>, std::set<Element_ *>> update_local_numbering(const PartitionIndex current_partition = 0);

		/// Sort the nodes and elements by the parent partition id (ascending)
		template<HAS_LOCAL_NUMBERING_DECLARE(Element_)>
		void sort_by_partition();

		///////////
		// Misc. //
		///////////

		/// Manually marks the nodes container as changed, which will force updates to
		/// subcomponent lists or KD trees when they are requested.
		/// The nodes container tries to track when nodes could be modified, but you
		/// get a raw pointer from a method and then change it, this should be called.
		void mark_nodes_as_changed();
		/// Manually marks the elements container as changed, which will force updates to
		/// subcomponent lists or KD trees when they are requested.
		/// The elements container tries to track when elements could be modified, but you
		/// get a raw pointer from a method and then change it, this should be called.
		void mark_elements_as_changed();

	private:
		// Containers
		MeshEntityContainer<Node_> nodes_;
		MeshEntityContainer<Element_> elements_;
	};

	template<typename Node_, typename Element_>
	inline Mesh<Node_, Element_>::Mesh(const Mesh & other)
	{
		add_to_self(other);
	}

	template<typename Node_, typename Element_>
	inline Mesh<Node_, Element_>::Mesh(Mesh && other) noexcept
	{
		*this = std::move(other);
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::add_to_self(const Mesh & other)
	{
		mark_nodes_as_changed();
		mark_elements_as_changed();

		std::map<Node_ *, Node_ *> other_to_new_node_map;
		for (auto& other_n_ptr : other.nodes_) {
			auto new_node = std::make_unique<Node_>(*other_n_ptr);
			this->nodes_.add(new_node);
			other_to_new_node_map[other_n_ptr->get()] = nodes_.back();
		}
		for (auto& other_e_ptr : other.elements_) {
			auto new_element = this->element_.add(std::make_unique<Element_>(*other_e_ptr));
			new_element->remap_element_nodes(other_to_new_node_map);
		}
		this->update_global_numbering();
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::move_to_self(Mesh & other)
	{
		mark_nodes_as_changed();
		mark_elements_as_changed();

		std::move(other.nodes_.objects_.begin(), other.nodes_.objects_.end(), std::back_inserter(this->nodes_.objects_));
		std::move(other.elements_.objects_.begin(), other.elements_.objects_.end(), std::back_inserter(this->elements_.objects_));
		update_global_numbering();
		other.clear();
	}

	template<typename Node_, typename Element_>
	inline Mesh<Node_, Element_>& Mesh<Node_, Element_>::operator+=(const Mesh & other) noexcept
	{
		this->add_to_self(other);
		return *this;
	}

	template<typename Node_, typename Element_>
	Mesh<Node_, Element_>& Mesh<Node_, Element_>::operator=(const Mesh & other) noexcept
	{
		this->clear();
		this->add_to_self(other);
		this->mark_nodes_as_changed();
		this->mark_elements_as_changed();
		return *this;
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::clear()
	{
		nodes_.objects_.clear();
		elements_.objects_.clear();
		mark_nodes_as_changed();
		mark_elements_as_changed();
	}

	template<typename Node_, typename Element_>
	inline const MeshEntityContainer<Node_>& Mesh<Node_, Element_>::nodes() const {
		return nodes_;
	}

	template<typename Node_, typename Element_>
	inline MeshEntityContainer<Node_>& Mesh<Node_, Element_>::nodes() {
		return nodes_;
	}

	template<typename Node_, typename Element_>
	inline Node_ * Mesh<Node_, Element_>::add_node_by_coord(const Point & coord) {
		auto n = std::unique_ptr<Node_>(new Node_(coord));
		return nodes_.add(std::move(n));
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::set_compare_tolerance(double tol) {
		Node_::set_compare_tolerance(tol);
	}

	template<typename Node_, typename Element_>
	inline double Mesh<Node_, Element_>::get_compare_tolerance() const {
		return Node_::get_compare_tolerance();
	}

	template<typename Node_, typename Element_>
	inline const MeshEntityContainer<Element_>& Mesh<Node_, Element_>::elements() const {
		return elements_;
	}

	template<typename Node_, typename Element_>
	inline MeshEntityContainer<Element_>& Mesh<Node_, Element_>::elements() {
		return elements_;
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::reset_element_faces_and_edges() {
		for (auto& e_ptr : elements_)
			e_ptr->reset_faces_and_edges();
	}

	template<typename Node_, typename Element_>
	inline std::set<ComponentIndex> Mesh<Node_, Element_>::get_component_ids() const
	{
		std::set<ComponentIndex> component_ids;
		for (auto& e_ptr : elements_)
			component_ids.insert(e_ptr->get_component_id());
		return component_ids;
	}

	template<typename Node_, typename Element_>
	inline std::vector<Element_*> Mesh<Node_, Element_>::get_component(ComponentIndex component_id) const
	{
		std::vector<Element_*> elements_in_component;
		for (auto& e_ptr : elements_) {
			if (e_ptr->get_component_id() == component_id) {
				elements_in_component.push_back(e_ptr.get());
			}
		}
		return elements_in_component;
	}

	template<typename Node_, typename Element_>
	template<HAS_LOCAL_NUMBERING_DEFINE(Element_)>
	inline unsigned int Mesh<Node_, Element_>::get_num_partitions() const
	{
		// Get the number of partitions first
		std::set<PartitionIndex> partitions;
		for (auto& n_ptr : nodes())
			partitions.insert(n_ptr->get_parent_partition());
		// Should we spend the time to ensure that the partitions are contiguous (no gaps)
		return static_cast<unsigned int>(partitions.size());
	}

	template<typename Node_, typename Element_>
	template<HAS_LOCAL_NUMBERING_DEFINE(Element_)>
	inline void Mesh<Node_, Element_>::update_global_numbering()
	{
		GlobalIndex current_num = 0;
		for (auto& n_ptr : nodes_) {
			n_ptr->set_global_number(current_num);
			n_ptr->set_local_number(current_num);
			current_num++;
		}
		current_num = 0;
		for (auto& e_ptr : elements_) {
			e_ptr->set_global_number(current_num);
			e_ptr->set_local_number(current_num);
			current_num++;
		}
	}

	template<typename Node_, typename Element_>
	template<MISSING_LOCAL_NUMBERING_DEFINE(Element_)>
		inline void Mesh<Node_, Element_>::update_global_numbering()
	{
		GlobalIndex current_num = 0;
		for (auto n_ptr : nodes_) {
			n_ptr->set_global_number(current_num);
			current_num++;
		}
		current_num = 0;
		for (auto e_ptr : elements_) {
			e_ptr->set_global_number(current_num);
			current_num++;
		}
	}

	template<typename Node_, typename Element_>
	template<HAS_LOCAL_NUMBERING_DEFINE(Element_)>
	inline std::pair<std::set<Node_ *>, std::set<Element_ *>>
		Mesh<Node_, Element_>::update_local_numbering(const PartitionIndex current_partition)
	{
		std::set<Node_ *> nodes_in_part;
		std::set<Element_ *> elements_in_part;
		// Loop over the elements
		for (auto& e_ptr : elements_) {
			// If one of the nodes is owned by this partition,
			// then all the nodes in the element to the partition need to be in
			// the partition's mesh file, and the element needs to be in the
			// partition's mesh file as well (regardless if the element is owned or not).
			for (auto& n_ptr : *e_ptr) {
				if (n_ptr->get_parent_partition() == current_partition) {
					elements_in_part.insert(e_ptr.get());
					for (auto& n_ptr2 : *e_ptr)
						nodes_in_part.insert(n_ptr2);
					break;
				}
			}
		} // End creating sets

		// Now renumber the nodes and elements in the sets
		unsigned int current_num = 0;
		for (auto node : nodes_in_part)
			node->set_local_number(current_num++);

		current_num = 0;
		for (auto element : elements_in_part)
			element->set_local_number(current_num++);

		return{ nodes_in_part, elements_in_part };
	}

	template<typename Node_, typename Element_>
	template<HAS_LOCAL_NUMBERING_DEFINE(Element_)>
	inline void Mesh<Node_, Element_>::sort_by_partition() {
		using NodeUniquePtr = std::unique_ptr<Node_>;
		using ElementUniquePtr = std::unique_ptr<Element_>;

		std::sort(nodes_.objects_.begin(), nodes_.objects_.end(),
			[](const NodeUniquePtr& n1, const NodeUniquePtr& n2) -> bool {
			return n1->get_parent_partition() < n2->get_parent_partition();
		});

		std::sort(elements_.objects_.begin(), elements_.objects_.end(),
			[](const ElementUniquePtr& e1, const ElementUniquePtr& e2) -> bool {
			return e1->get_parent_partition() < e2->get_parent_partition();
		});
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::mark_nodes_as_changed()
	{
		nodes_.mark_as_changed();
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::mark_elements_as_changed()
	{
		elements_.mark_as_changed();
	}
}

#undef HAS_LOCAL_NUMBERING_DECLARE
#undef HAS_LOCAL_NUMBERING_DEFINE
#undef MISSING_LOCAL_NUMBERING_DECLARE
#undef MISSING_LOCAL_NUMBERING_DEFINE
