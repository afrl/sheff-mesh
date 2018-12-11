#pragma once

#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/Exceptions.h"

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace flex_mesh_detail
{
// This is a bit of complicated macro that tests if a class has a particular method
// at compile time, allowing specialized methods based on whether the method exists.
// A great idea borrowed from https://en.wikibooks.org/wiki/More_C++_Idioms/Member_Detector
// Consider moving this to Common.h?
#define GENERATE_HAS_MEMBER(member)                                               \
                                                                                  \
template < class T >                                                              \
class HasMember_##member                                                          \
{                                                                                 \
private:                                                                          \
    using Yes = char[2];                                                          \
    using  No = char[1];                                                          \
                                                                                  \
    struct Fallback { int member; };                                              \
    struct Derived : T, Fallback { };                                             \
                                                                                  \
    template < class U >                                                          \
    static No& test ( decltype(U::member)* );                                     \
    template < typename U >                                                       \
    static Yes& test ( U* );                                                      \
                                                                                  \
public:                                                                           \
    static constexpr bool RESULT = sizeof(test<Derived>(nullptr)) == sizeof(Yes); \
};                                                                                \
                                                                                  \
template < class T >                                                              \
struct has_##member                                                        \
: public std::integral_constant<bool, HasMember_##member<T>::RESULT>              \
{                                                                                 \
};

GENERATE_HAS_MEMBER(get_local_number)  // Creates 'has_get_local_number'.

}

// Tests to see if elements contain distributed information
#define HAS_LOCAL_NUMBERING_DECLARE(T) typename T_ = T, typename std::enable_if<flex_mesh_detail::has_get_local_number<T_>::value>::type* = nullptr
#define HAS_LOCAL_NUMBERING_DEFINE(T) typename T_, typename std::enable_if<flex_mesh_detail::has_get_local_number<T_>::value>::type*
#define MISSING_LOCAL_NUMBERING_DECLARE(T) typename T_ = T, typename std::enable_if<!flex_mesh_detail::has_get_local_number<T_>::value>::type* = nullptr
#define MISSING_LOCAL_NUMBERING_DEFINE(T) typename T_, typename std::enable_if<!flex_mesh_detail::has_get_local_number<T_>::value>::type*

namespace flex_mesh
{
	template<typename T>
	class UniquePtrContainer
	{
	public:
		using ObjectUniquePtr = std::unique_ptr<T>;
		typedef typename std::vector<ObjectUniquePtr>::iterator Iterator;
		typedef typename std::vector<ObjectUniquePtr>::const_iterator ConstIterator;

		UniquePtrContainer() {};
		virtual ~UniquePtrContainer() {};

		/// Returns the size of the container
		virtual size_t size() const { return objects_.size(); }
		/// Returns the begin iterator of the container
		virtual auto begin() { return objects_.begin(); }
		/// Returns the begin iterator of the container
		virtual auto begin() const { return objects_.begin(); }
		/// Returns the end iterator of the container
		virtual auto end() const { return objects_.end(); }
		/// Returns a raw pointer to the first object added
		T * const front() const { return  !objects_.empty() ? objects_.front().get() : nullptr; }
		/// Returns a raw pointer to the last object added
		T * const back() const { return !objects_.empty() ? objects_.back().get() : nullptr; }
		/// Add an object to the container
		/// Returns a raw pointer to that object
#ifndef SWIG
		T * add(ObjectUniquePtr obj) {
			objects_.push_back(std::move(obj));
			return back();
		};
#endif
		/// Add an object to the container and have it call its default constructor
		/// Returns a raw pointer to that object
		T * add() {
			objects_.push_back(std::make_unique<T>());
			return back();
		};
		/// Remove an object from the container given an iterator
		/// Returns the next iterator
		Iterator remove(Iterator it) {
			*it = std::move(objects_.back());
			objects_.resize(objects_.size() - 1);
			return it++;
		};
		/// Remove a list of objects from the container given iterators
		void remove(std::list<Iterator> iter_list) {
			for (auto it : iter_list)
				this->remove(it);
		}
		/// Given a set of object pointers, this will remove all of them from the container.
		/// Note: if the container stores elements, it does not update element connectivities,
		/// that must be done in the algorithm that call this.
		void remove_nodes(std::set<T *> objects_to_remove) {
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

	protected:
		std::vector<ObjectUniquePtr> objects_;
	};

	template<typename Element_, typename PolyhedralElement_>
	class ElementContainer : public UniquePtrContainer<Element_>
	{
	public:
		using Parent = UniquePtrContainer<Element_>;
		using PolyhedronUniquePtr = std::unique_ptr<PolyhedralElement_>;

		ElementContainer() : Parent() {};
		virtual ~ElementContainer() override {};



	protected:
		std::vector<PolyhedronUniquePtr> polyhedral_elements_;
	};

	// Iterator adapter for elements to iterate through elements then polyhedrals?



	/// A base container for NodeTopology and ElementTopology objects
	/// that compose a mesh. This is meant to be very simple and only
	/// contain methods related to accessing and manipulating a mesh.
	/// All algorithms will be found in other classes that derive
	template<typename Node_, typename Element_>
	class Mesh
	{
	public:
		using NodeUniquePtr = std::unique_ptr<Node_>;
		using ElementUniquePtr = std::unique_ptr<Element_>;

		using NodeVector = std::vector<NodeUniquePtr>;
		typedef typename NodeVector::iterator NodeVectorIterator;
		typedef typename NodeVector::const_iterator NodeVectorConstIterator;
		using ElementVector = std::vector<ElementUniquePtr>;
		typedef typename ElementVector::iterator ElementVectorIterator;
		typedef typename ElementVector::const_iterator ElementVectorConstIterator;

		//////////////////
		// Constructors //
		//////////////////

		Mesh() = default;
		/// A deepcopy constructor
		Mesh(const Mesh& other);
		/// Move constructor
		Mesh(Mesh&& other) noexcept;
		/// Move assignment operator.
		Mesh& operator=(Mesh&&) noexcept = default;

		/// Clears the mesh object of its contents
		void clear();

		/////////////////////
		// Node operations //
		/////////////////////

		/// Returns the size of the node container
		size_t nodes_size() const { return nodes_.size(); }
		/// Returns the begin iterator of the node container
		auto nodes_begin() { return nodes_.begin(); }
		/// Returns the begin iterator of the node container
		auto nodes_begin() const { return nodes_.begin(); }
		/// Returns the end iterator of the node container
		auto nodes_end() const { return nodes_.end(); }
		/// Returns a pointer to the first node added
		Node_ * nodes_front() const { return !nodes_.empty() ? nodes_.front().get() : nullptr; }
		/// Returns a pointer to the last node added
		Node_ * nodes_back() const { return !nodes_.empty() ? nodes_.back().get() : nullptr; }
		/// Returns a const reference to the nodes in the mesh (useful for range based loops).  Example:
		/// \code{.cpp} for (auto& n : mesh.nodes()) n->do_something(); \endcode
		const NodeVector& nodes() const { return nodes_; }
		/// Add a node to the container
		/// Returns a pointer to that node
#ifndef SWIG
		Node_ * add_node(NodeUniquePtr& n) {
			nodes_.push_back(std::move(n));
			return nodes_back();
		};
#endif
        /// Add a node via the coordinate
        Node_ * add_node_by_coord(const Point& coord) {
            auto n = std::unique_ptr<Node_>(new Node_(coord));
            return add_node(n);
        };
		/// Remove a node from the container given an iterator
        /// Returns the next iterator
        typename NodeVectorIterator remove_node(NodeVectorIterator n_it) {
			*n_it = std::move(nodes_.back());
			nodes_.resize(nodes_.size() - 1);
            return n_it++;
        };
		/// Given a set of node pointers, this will remove all of them from the mesh.
		/// Note: it does not update element connectivities, that must be done in the
		/// algorithm that call this.
		void remove_nodes(std::set<Node_ *> nodes_to_remove) {
			for (auto nIt = nodes_begin(); nIt != nodes_end(); ) {
				auto setIt = nodes_to_remove.find(nIt->get());
				if (setIt != nodes_to_remove.end()) {
					nodes_to_remove.erase(setIt);
					nIt = remove_node(nIt);
				}
				else
					nIt++;
			}
		};
		void set_compare_tolerance(double tol) {
			Node_::set_compare_tolerance(tol);
		}
        double get_compare_tolerance() const {
            return Node_::get_compare_tolerance();
        }

		////////////////////////
		// Element operations //
		////////////////////////

		/// Returns the size of the element container
		size_t elements_size() const { return elements_.size(); }
		/// Returns the begin iterator of the element container
		ElementVectorIterator elements_begin() { return elements_.begin(); }
		/// Returns the begin iterator of the element container
		ElementVectorConstIterator elements_begin() const { return elements_.begin(); }
		/// Returns the end iterator of the element container
		ElementVectorConstIterator elements_end() const { return elements_.end(); }
		/// Returns a pointer to the first element added
		Element_ * elements_front() const { return  !elements_.empty() ? elements_.front().get() : nullptr; }
		/// Returns a pointer to the last element added
		Element_ * elements_back() const { return !elements_.empty() ? elements_.back().get() : nullptr; }
		/// Returns a const reference to the nodes in the mesh (useful for range based loops).  Example:
		/// \code{.cpp} for (auto& e : mesh.elements()) e->do_something(); \endcode
		const ElementVector& elements() const { return elements_; }
		/// Add an element to the container
		/// Returns a shared_ptr to that element
#ifndef SWIG
		Element_ * add_element(ElementUniquePtr e) {
			elements_.push_back(std::move(e));
			return elements_back();
		};
#endif
        /// Add a default element to the container
        /// Returns a shared_ptr to that element
		Element_ * add_element() {
            elements_.push_back(std::make_unique<Element_>());
            return elements_back();
        };
		/// Remove an element from the container given an iterator
        /// Returns the next iterator
		ElementVectorIterator remove_element(ElementVectorIterator e_it) {
			*e_it = std::move(elements_.back());
			elements_.resize(elements_.size() - 1);
			return e_it++;
        };
		/// Remove a list of elements from the container given iterators
		void remove_elements(std::list<ElementVectorIterator> iter_list) {
			for (auto e_it : iter_list)
				this->remove_element(e_it);
		}
        /// Reset element face/edges
        void reset_element_faces_and_edges() {
            for (auto& e_ptr : elements_)
				e_ptr->reset_faces_and_edges();
        }

		//////////////////////////
		// Subcomponent methods //
		//////////////////////////

		/// Updates teh subcomponents of the mesh.
		/// Since the component ID of elements can be modified outside the methods in the mesh,
		/// this method must be called to update the element groups based on component ID's.
		void update_subcomponents();
		/// Returns a vector of element const pointers to the subcomponent.
		/// Be aware that calling this may result in unexpected elements if you modify the
		/// component id of elements via the Element::set_component_id(...) method
		/// However, the removal and addition of elements are internally tracked to keep the
		/// subcomponents up to date as elements are inserted and removed.
		const std::vector<Element_ * const>& get_subcomponent(const ComponentIndex component_id);

		//////////////////////////
		// Mesh-Mesh operations //
		//////////////////////////

		/// Add another to this mesh, making a deepcopy of nodes and elements
		void add_to_self(const Mesh & other);
        /// Move another mesh into this mesh
        /// \important Note: the other mesh will be empty after this
        void move_to_self(Mesh & other);
		/// Destroys the current mesh and copies the other mesh.
		Mesh<Node_, Element_>& operator=(const Mesh& other) {
			this->clear();
			this->add_to_self(other);
			return *this;
		}

		//////////////////////////////
		// Distributed mesh methods //
		//////////////////////////////

		// Note to devs: should this be cleaned up by creating a derived type?  This is clever
		// in that methods "disappear" if the elements do not support the distributed functionality,
		// but it does introduce some complicated templates and type traits.  C++20 should clean this
		// up significantly with constexpr (if the standard passes as expected as of writing this).
		// One could create a derived mesh type and add a static_assert to check if the nodes/elements
		// contain a get_local_number method instead of using SPINAE to get rid of the methods.

		/// Returns the number of partitions (assuming they are contiguous 0-based when assigned)
		// TODO add documentation explaining types of valid elements
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
		
	protected:
		// Containers
		std::vector<NodeUniquePtr> nodes_;
		std::vector<ElementUniquePtr> elements_;

		bool subcomponts_up_to_date;
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
	inline void Mesh<Node_, Element_>::clear()
	{
		nodes_.clear();
		elements_.clear();
	}

	template<typename Node_, typename Element_>
	inline void Mesh<Node_, Element_>::add_to_self(const Mesh & other)
	{
		std::map<Node_ *, Node_ *> other_to_new_node_map;
		for (auto& other_n_ptr : other.nodes_) {
			auto new_node = std::make_unique<Node_>(*other_n_ptr);
			this->add_node(new_node);
			other_to_new_node_map[other_n_ptr->get()] = nodes_back();
		}
		for (auto& other_e_ptr : other.elements_) {
			auto new_element = this->add_element(std::make_unique<Element_>(*other_e_ptr));
			new_element->remap_element_nodes(other_to_new_node_map);
		}
		this->update_global_numbering();
	}

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::move_to_self(Mesh & other)
    {
		std::move(other.nodes_.begin(), other.nodes_.end(), std::back_inserter(this->nodes_));
		std::move(other.elements_.begin(), other.elements_.end(), std::back_inserter(this->elements_));
        update_global_numbering();
        other.clear();
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
		return partitions.size();
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
		std::sort(nodes_.begin(), nodes_.end(),
			[](const NodeUniquePtr& n1, const NodeUniquePtr& n2) -> bool {
			return n1->get_parent_partition() < n2->get_parent_partition();
		});

		std::sort(elements_.begin(), elements_.end(),
			[](const ElementUniquePtr& e1, const ElementUniquePtr& e2) -> bool {
			return e1->get_parent_partition() < e2->get_parent_partition();
		});
	}
	
#undef HAS_LOCAL_NUMBERING
#undef MISSING_LOCAL_NUMBERING

};