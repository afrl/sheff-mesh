#pragma once
/*
#include "PrimativeMesh.h"
#include "nanoflann/nanoflannHelpers.h"

namespace BetaMesh
{
	/// An extension of the mesh class that provides:
	/// - Vector storage of nodes and elments that can be turned
	///   on and off, allowing fast access by local index. Off by default.
	/// - A map of nodes to their parent element, which can be turned
	///   on and off.  Off by default.
	/// - Groups
	template<typename NodeType_, typename ElementType_>
	class Mesh : public PrimativeMesh<NodeType_, ElementType_>
	{
	public:
		typedef PrimativeMesh<NodeType_, ElementType_> SuperClass;
        typedef NodeType_ NodeType;
        typedef ElementType_ ElementType;
        typedef typename ElementType_::FaceType FaceType;
        typedef typename ElementType_::EdgeType EdgeType;
        typedef NodeType* NodePtr;
        typedef EdgeType* EdgePtr;
        typedef FaceType* FacePtr;
        typedef ElementType* ElementPtr;

		typedef std::unique_ptr<NodeType> NodeUniquePtr;
		typedef std::unique_ptr<ElementType> ElementUniquePtr;

        typedef std::list<std::unique_ptr<NodeType>> UniquePtrNodeList;
        typedef typename std::list<std::unique_ptr<NodeType>>::iterator UniquePtrNodeListIt;
        typedef std::list<std::unique_ptr<ElementType>> UniquePtrElementList;
        typedef typename std::list<std::unique_ptr<ElementType>>::iterator UniquePtrElementListIt;

		typedef std::vector<NodePtr> NodeVector;
		typedef typename std::vector<NodePtr>::iterator NodeVectorIt;
		typedef std::vector<ElementPtr> ElementVector;
		typedef typename std::vector<ElementPtr>::iterator ElementVectorIt;

		typedef std::forward_list<ElementPtr> ElementForwardList;

		typedef std::set<NodePtr, ObjectComparisonWithNumber<NodeType>> NodeGroup;
		typedef std::map<std::string, NodeGroup> NodeGroupMap;
		typedef typename NodeGroupMap::iterator NodeGroupMapIt;
		typedef std::set<ElementPtr> ElementGroup;
		typedef std::map<std::string, ElementGroup> ElementGroupMap;
		typedef typename std::map<std::string, ElementGroup>::iterator ElementGroupMapIt;

		// For KD trees
		// Note: the template requires the dimensionality a priori, so we will tell it that
		// it is 3D, and the adapter will fake the node coordinates being 3D by giving a 0
		// for z when requested.
		typedef nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<double, KDTreeAdapter<NodePtr, std::vector<double>> >,
			KDTreeAdapter<NodePtr, std::vector<double>>,
			3 > NodeKDTreeType; // dimension
        typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<double, KDTreeAdapter<ElementPtr, std::vector<double>> >,
            KDTreeAdapter<ElementPtr, std::vector<double>>,
            3  > ElementKDTreeType; // dimension

	protected:
		// These vectors may or may not contain anything depending on storing_vectors_
		NodeVector nodes_vec_;
		ElementVector elements_vec_;
		// Indicates if pointers will be stored in vectors as well as lists.
		bool storing_vectors_;

		// Note: KD structures will use unique pointers to maintain ownership
		// and allow return value optimization to work.  The rvalue move does not
		// know what to do with the KD tree since they did not implement an rvalue
		// move themselves.  So, we will use a pointer to make sure there is one
		// object but bypass the needed to implement our own move.
		// A KD adaptor that builds on top of the vector
		std::unique_ptr<KDTreeAdapter<NodePtr, std::vector<double>>> node_vec_kd_adapter_;
		// A KD tree structure for fast access to nodes nearest a point, nearest
		// neighbor searches, and radius searches around a point.
		std::unique_ptr<NodeKDTreeType> node_kd_tree_;
		// Storing this makes it far faster for subsequent calls to get the nearest neighbor to a point
		std::unique_ptr<KDSingleNeighbor> node_kd_single_nn_;
		// A bool flag inidicating if the KD index structure is up to date with the vector
		bool node_kd_up_to_date_;
        // A KD adaptor that builds on top of the vector
        std::unique_ptr<KDTreeAdapter<ElementPtr, std::vector<double>>> element_vec_kd_adapter_;
        // A KD tree structure for fast access to elements nearest a point, nearest
        // neighbor searches, and radius searches around a point.
        std::unique_ptr<ElementKDTreeType> element_kd_tree_;
        // Storing this makes it far faster for subsequent calls to get the nearest neighbor to a point
        std::unique_ptr<KDSingleNeighbor> element_kd_single_nn_;
        // A bool flag inidicating if the KD index structure is up to date with the vector
        bool element_kd_up_to_date_;

		// Containers for parent relations
		std::map<NodePtr, ElementForwardList> node_parent_elements_;
		bool parents_up_to_date_;

		// Groups
		std::map<std::string, NodeGroup> node_groups_;
		std::map<std::string, ElementGroup> element_groups_;

	public:
		// Constructors
		/// Default constructor
		Mesh() : PrimativeMesh<NodeType_, ElementType_>(), storing_vectors_(false),
			parents_up_to_date_(false),
            node_vec_kd_adapter_(new KDTreeAdapter<NodePtr, std::vector<double>>(nodes_vec_, KDFunctors::KDNodeFunctor<NodePtr>())),
			node_kd_tree_(new NodeKDTreeType(3 /*dim*/, *node_vec_kd_adapter_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */))),
			node_kd_single_nn_(new KDSingleNeighbor()),
            node_kd_up_to_date_(false),
            element_vec_kd_adapter_(new KDTreeAdapter<ElementPtr, std::vector<double>>(elements_vec_, KDFunctors::KDElementFunctor<ElementPtr>())),
            element_kd_tree_(new ElementKDTreeType(3 /*dim*/, *element_vec_kd_adapter_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */))),
            element_kd_single_nn_(new KDSingleNeighbor()),
            element_kd_up_to_date_(false) {};
		/// Copy constructor
		Mesh(const Mesh<NodeType_, ElementType_>& other) : Mesh() {
			this->add_to_self(other);
		}
		/// Move constructor
		Mesh(Mesh&& other) noexcept {
			*this = std::move(other);
		}
		/// Move assignment operator.
		Mesh& operator=(Mesh&&) noexcept = default;

        /// Clears the mesh object of the vectors, parent lookups,
        /// and both KD trees.
        void clear_vectors_parents_and_kd() {
            do_not_store_vectors();
            node_parent_elements_.clear();
            parents_up_to_date_ = false;
            free_kd();
        }

		/// Clears the mesh object of its contents
		void clear() {
			nodes_vec_.clear();
			elements_vec_.clear();
			node_parent_elements_.clear();
			node_groups_.clear();
			element_groups_.clear();
			parents_up_to_date_ = false;
			free_kd();
			node_kd_up_to_date_ = false;
            element_kd_up_to_date_ = false;

			SuperClass::clear();
		}

		///////////////////
		// KD Operations //
		///////////////////

	protected:
		/// A method for rebuilding the node KD tree, which will be called by the query methods below.
		/// I decided to make this method protected becuase I do not see a reason to expose it outside
		/// this class, since it will be called as needed internally.
		void rebuild_node_kd() {
			if (!storing_vectors_) // If vectors are not available, then we can't do a KD search.
				throw _UTIL_MSG_EXCEPT("Cannot build KD tree if vectors are not being stored, since it cannot operator on a linked list.  Call build_vectors first.");
			if (nodes_vec_.size() == 0)
				throw _UTIL_MSG_EXCEPT("Cannot build KD tree if the vector of nodes is empty.  Does the mesh have nodes?");
			node_vec_kd_adapter_->initialize(SuperClass::nodes_front()->get_dimensionality());
			node_kd_tree_->buildIndex();
			node_kd_up_to_date_ = true;
		}

        /// A method for rebuilding the KD tree, which will be called by the query methods below.
        /// I decided to make this method protected becuase I do not see a reason to expose it outside
        /// this class, since it will be called as needed internally.
        void rebuild_element_kd() {
            if (!storing_vectors_) // If vectors are not available, then we can't do a KD search.
                throw _UTIL_MSG_EXCEPT("Cannot build KD tree if vectors are not being stored, since it cannot operator on a linked list.  Call build_vectors first.");
            if (elements_vec_.size() == 0)
                throw _UTIL_MSG_EXCEPT("Cannot build KD tree if the vector of elements is empty.  Does the mesh have elements?");
            element_vec_kd_adapter_->initialize(SuperClass::nodes_front()->get_dimensionality());
            element_kd_tree_->buildIndex();
            element_kd_up_to_date_ = true;
        }

	public:
		/// Forces the KD tree to free its memory, but it will be rebuilt with the next KD call.
		void free_kd() {
			node_kd_tree_->freeIndex();
			node_kd_up_to_date_ = false;
            element_kd_tree_->freeIndex();
            element_kd_up_to_date_ = false;
		}

		/// Get the nearest node to a point
		/// IMPORTANT: queryPoint MUST be an array of length 3!!!
		NodePtr get_nearest_node(const double* queryPoint) {
			if (!node_kd_up_to_date_) rebuild_node_kd();
			node_kd_tree_->findNeighbors(node_kd_single_nn_->result_set_, &queryPoint[0], nanoflann::SearchParams(10));
			return nodes_vec_[node_kd_single_nn_->index()];
		}

		/// Get the nearest node to a point, but this returns the node pointer and the distance from the query point
		/// IMPORTANT: queryPoint MUST be an array of length 3!!!
		std::pair<NodePtr, double> get_nearest_node_with_dist(const double* queryPoint) {
			if (!node_kd_up_to_date_) rebuild_node_kd();
			node_kd_tree_->findNeighbors(node_kd_single_nn_->result_set_, &queryPoint[0], nanoflann::SearchParams(10));
			return std::pair<NodePtr, double>(nodes_vec_[node_kd_single_nn_->index()], node_kd_single_nn_->distance());
		}

		/// Get an arbitrary number of nodes nearest a point.  Returns the node pointers.
		/// Note: it is not guarenteed to find as many as you request
		/// IMPORTANT: queryPoint MUST be an array of length 3!!!
		NodeVector get_nearest_nodes(const double* queryPoint, const size_t numNeighbors) {
			if (!node_kd_up_to_date_) rebuild_node_kd();
			// Perform the search
			std::vector<size_t> nodeIndices(numNeighbors);
			std::vector<double> nodeDistances(numNeighbors);
			const size_t numResults = node_kd_tree_->knnSearch(&queryPoint[0], numNeighbors, &nodeIndices[0], &nodeDistances[0]);
			// Create the resulting vector of node pointers
			NodeVector nearestNodes(numResults);
			for (size_t i = 0; i < numResults; i++)
				nearestNodes[i] = nodes_vec_[nodeIndices[i]];
			return nearestNodes;
		}

		/// Get the nodes within a given radius of a point.  Returns the node pointers.
		/// Note: does not include nodes that are right on the search radius.
		/// IMPORTANT: queryPoint MUST be an array of length 3!!!
		NodeVector get_nodes_within_radius(const double* queryPoint, const double radius, const bool sortedByDistance = false) {
			if (!node_kd_up_to_date_) rebuild_node_kd();
			// Perform the search
			std::vector<std::pair<size_t, double> > matches;
			nanoflann::SearchParams params;
			params.sorted = sortedByDistance;
			const size_t numResults = node_kd_tree_->radiusSearch(&queryPoint[0], radius, matches, params);
			// Create the resulting vector of node pointers
			NodeVector nearestNodes(numResults);
			for (size_t i = 0; i < numResults; i++)
				nearestNodes[i] = nodes_vec_[matches[i].first];
			return nearestNodes;
		}

        /// Get the nearest element to a point
        /// IMPORTANT: queryPoint MUST be an array of length 3!!!
        ElementPtr get_nearest_element(const double* queryPoint) {
            if (!element_kd_up_to_date_) rebuild_element_kd();
            element_kd_tree_->findNeighbors(element_kd_single_nn_->result_set_, &queryPoint[0], nanoflann::SearchParams(10));
            return elements_vec_[element_kd_single_nn_->index()];
        }

        /// Get the nearest element to a point, but this returns the element pointer and the distance from the query point
        /// IMPORTANT: queryPoint MUST be an array of length 3!!!
        std::pair<ElementPtr, double> get_nearest_element_with_dist(const double* queryPoint) {
            if (!element_kd_up_to_date_) rebuild_element_kd();
            element_kd_tree_->findNeighbors(element_kd_single_nn_->result_set_, &queryPoint[0], nanoflann::SearchParams(10));
            return std::pair<ElementPtr, double>(elements_vec_[element_kd_single_nn_->index()], element_kd_single_nn_->distance());
        }

        /// Get an arbitrary number of elements nearest a point.  Returns the element pointers.
        /// Note: it is not guarenteed to find as many as you request
        /// IMPORTANT: queryPoint MUST be an array of length 3!!!
        ElementVector get_nearest_elements(const double* queryPoint, const size_t numNeighbors) {
            if (!element_kd_up_to_date_) rebuild_element_kd();
            // Perform the search
            std::vector<size_t> elementIndices(numNeighbors);
            std::vector<double> elementDistances(numNeighbors);
            const size_t numResults = element_kd_tree_->knnSearch(&queryPoint[0], numNeighbors, &elementIndices[0], &elementDistances[0]);
            // Create the resulting vector of element pointers
            ElementVector nearestElements(numResults);
            for (size_t i = 0; i < numResults; i++)
                nearestElements[i] = elements_vec_[elementIndices[i]];
            return nearestElements;
        }

        /// Get the elements within a given radius of a point.  Returns the element pointers.
        /// Note: does not include elements that are right on the search radius.
        /// IMPORTANT: queryPoint MUST be an array of length 3!!!
        ElementVector get_elements_within_radius(const double* queryPoint, const double radius, const bool sortedByDistance = false) {
            if (!element_kd_up_to_date_) rebuild_element_kd();
            // Perform the search
            std::vector<std::pair<size_t, double> > matches;
            nanoflann::SearchParams params;
            params.sorted = sortedByDistance;
            const size_t numResults = element_kd_tree_->radiusSearch(&queryPoint[0], radius, matches, params);
            // Create the resulting vector of element pointers
            ElementVector nearestElements(numResults);
            for (size_t i = 0; i < numResults; i++)
                nearestElements[i] = elements_vec_[matches[i].first];
            return nearestElements;
        }

		////////////////////////////
		// Vector storage methods //
		////////////////////////////

		/// Normally, the lists will be used, but in anticipation of applications where
		/// nodes and elements will need to be reference by their local numbering,
		/// this will build a vector of pointers rather than a list.  This effectively
		/// provides a container with very fast random access (which the list class does not
		/// not provide), but removing nodes and elements are disabled.
		/// Note: the vectors will never be built unless something calls this method first,
		/// it is not a performance hit if you do not need to access things by local number.
		void build_vectors() {
			if (storing_vectors_) return;
			storing_vectors_ = true;
			nodes_vec_.resize(SuperClass::nodes_.size());
			for (auto& n : SuperClass::nodes_)
				nodes_vec_[n->get_local_number()] = n.get();
			elements_vec_.resize(SuperClass::elements_.size());
			for (auto& e : SuperClass::elements_)
				elements_vec_[e->get_local_number()] = e.get();

			// Check consistency
			for (auto n : nodes_vec_)
				if (n == nullptr) throw _BETA_MSG_EXCEPT_C2("Something was wrong with the local numbering the resulted in null nodes in the vector.");
			for (auto e : elements_vec_)
				if (e == nullptr) throw _BETA_MSG_EXCEPT_C2("Something was wrong with the local numbering the resulted in null nodes in the vector.");
		}

		/// Once build_vectors is called, they disable removing elements, so this method
		/// can be used to throw away the vectors and only use the lists.
		void do_not_store_vectors() {
			storing_vectors_ = false;
			nodes_vec_.clear();
			elements_vec_.clear();
		}

		bool is_storing_vectors() { return storing_vectors_; };

		// Node operation overrides for nodes

		/// Override to ensure vectors do not get messed up
#ifndef SWIG
		NodePtr add_node(NodeUniquePtr n) {
			if (parents_up_to_date_) parents_up_to_date_ = false;
			auto raw_ptr = PrimativeMesh<NodeType, ElementType>::add_node(n);
			if (storing_vectors_) nodes_vec_.push_back(raw_ptr);
			return raw_ptr;
		};
#endif
		/// Override to ensure vectors do not get messed up
        typename UniquePtrNodeList::iterator remove_node(typename UniquePtrNodeList::iterator n_it) {
			if (storing_vectors_) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "Cannot remove an object when storing vectors.");
			if (parents_up_to_date_) parents_up_to_date_ = false;
			return PrimativeMesh<NodeType, ElementType>::remove_node(n_it);
		};
		/// Returns a node by its local index
		NodePtr get_node(GlobalIndex index) {
			if (!storing_vectors_) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "The vectors were never built, call build_vectors first.");
			return nodes_vec_[index];
		}

		// Element operation overrides for vectors

		/// Override to ensure vectors do not get messed up
#ifndef SWIG
		ElementPtr add_element(ElementUniquePtr e) {
			if (parents_up_to_date_) parents_up_to_date_ = false;
			auto e_ptr = PrimativeMesh<NodeType, ElementType>::add_element(std::move(e));
			if (storing_vectors_) elements_vec_.push_back(e_ptr);
			return e_ptr;
		};
#endif
		/// Override to ensure vectors do not get messed up
		ElementPtr add_element() {
			if (parents_up_to_date_) parents_up_to_date_ = false;
			auto e = PrimativeMesh<NodeType, ElementType>::add_element();
			if (storing_vectors_) elements_vec_.push_back(e);
			return e;
		};
		/// Override to ensure vectors do not get messed up
		UniquePtrElementListIt remove_element(UniquePtrElementListIt e_it) {
			if (storing_vectors_) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "Cannot remove an object when storing vectors.");
			if (parents_up_to_date_) {
				// Remove this element from the parents of each of it's nodes
				for (auto& n : e_it->get()->get_conn())
					node_parent_elements_[n].remove(e_it->get());
			}
			return PrimativeMesh<NodeType, ElementType>::remove_element(e_it);
		};

        /// Removes all elements that satisfy a criteria, which is tested by the functor provided
        void remove_element_if(std::function<bool(const ElementUniquePtr&)> criteria) {
            if (storing_vectors_) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "Cannot remove an object when storing vectors.");
            if (parents_up_to_date_) parents_up_to_date_ = false;
            SuperClass::elements_.remove_if(criteria);
        }

		/// Returns a node by its local index
		ElementPtr get_element(GlobalIndex index) {
			if (!storing_vectors_) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "The vectors were never built, call build_vectors first.");
			return elements_vec_[index];
		}

		/////////////////////////
		// Node-Parent methods //
		/////////////////////////

		// Design note: for CZ's we will need the ability to update one node at a time (TODO)

        /// Updates the parent elements for each node (clears the old list completely)
        /// \attention { This will not be maintained during add or removal of nodes, so the algorithm using this
        /// functionality is responsible for making sure updates are called appropriately. }
        void update_parents()
        {
            node_parent_elements_.clear();
            for (auto& n : SuperClass::nodes_)
                node_parent_elements_[n.get()] = std::forward_list<ElementPtr>();
            for (auto& e : SuperClass::elements_) {
                for (auto n : e->get_conn())
                    node_parent_elements_[n].push_front(e.get());
            }
            parents_up_to_date_ = true;
        }

        /// Returns the nodes parent elements by const reference
        const ElementForwardList& get_parent_elements_of_node(const NodePtr& n) {
            if (node_parent_elements_.size() == 0) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "To get the parents of a node, update parents must be called first.");
            return node_parent_elements_[n];
        };

        /// Just indicates if the current node parents map is up to date
        bool is_parents_up_to_date() { return parents_up_to_date_; };

        /// When doing duplicate node removal or combining nodes within an algorithms, call this
        /// will allow you to preserve the parent list
        void adjust_parents_for_node_removal(const NodePtr& oldNode, const NodePtr& newNode) {
            if (node_parent_elements_.size() == 0) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "To adjust the parents of a node, update parents must be called first.");
            std::set<ElementPtr> combinedParents;
            for (auto& p : node_parent_elements_[oldNode])
                combinedParents.insert(p);
            for (auto& p : node_parent_elements_[newNode])
                combinedParents.insert(p);
            node_parent_elements_.erase(oldNode);
            node_parent_elements_[newNode].clear();
            for (auto& p : combinedParents)
                node_parent_elements_[newNode].push_front(p);
        };

		/// Updates the parent elements for each node (clears the old list completely)
		/// \attention { This will not be maintained during add or removal of nodes, so the algorithm using this
		/// functionality is responsible for making sure updates are called appropriately. }
		void update_parents()
		{
            for (auto n : SuperClass::nodes_)
                n->reset_parents();
			for (auto& e : SuperClass::elements_) {
				for (auto& n : e->get_conn())
                    n->add_parent_element<ElementType_>(e.get());
			}
			parents_up_to_date_ = true;
		}

		/// Returns the nodes parent elements by const reference
		const ElementForwardList& get_parent_elements_of_node(NodePtr n) {
			if (node_parent_elements_.size() == 0) throw _BETA_MSG_EXCEPT(BetaErrorCode::ERR_LOGIC_GENERAL, "To get the parents of a node, update parents must be called first.");
			return node_parent_elements_[n];
		};

		/// Just indicates if the current node parents map is up to date
		bool is_parents_up_to_date() { return parents_up_to_date_; };

        //////////////////////
		// Grouping methods //
		//////////////////////

		// Node Groups

		/// Create a new node group.  If it already exists, then it will do nothing.  Returns if one is created.
		bool create_node_group(const std::string& groupName)
		{
			auto group = node_groups_.find(groupName);
			if (group != node_groups_.end())
				return false;
			node_groups_[groupName] = NodeGroup();
			return true;
		}
		/// Gets rid of an node group.  Returns true if there was a group that was deleted
		bool remove_node_group(const std::string& groupName)
		{
			auto group = node_groups_.find(groupName);
			if (group == node_groups_.end())
				return false;
			node_groups_.erase(group);
			return true;
		}
		/// Adds a single node to a group.  Will throw if the group doesn't exist.
		void add_node_to_group(const std::string& groupName, NodePtr n)
		{
			auto group = node_groups_.find(groupName);
			if (group == node_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			group->second.insert(n);
		}
		/// Adds a set of nodes to a group.  Will throw if the group doesn't exist.
		void add_nodes_to_group(const std::string& groupName, std::set<NodePtr> nSet)
		{
			auto group = node_groups_.find(groupName);
			if (group == node_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			for (auto n : nSet)
				group->second.insert(n);
		}
		/// Adds all the elements in the mesh to a group.  Will throw if the group doesn't exist.
		void add_all_nodes_to_group(const std::string& groupName)
		{
			auto group = node_groups_.find(groupName);
			if (group == node_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			for (auto& n : SuperClass::nodes_)
				group->second.insert(n.get());
		}
		/// Gets rid of an node group.  Will throw if the group doesn't exist.
		void remove_node_from_group(const std::string& groupName, NodePtr e)
		{
			auto group = node_groups_.find(groupName);
			if (group == node_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			group->second.erase(e);
		}
		/// Returns the node by const reference.  It can be looped over and the nodes can be modified
		/// but the set itself cannot be modified.
		/// Because in distributed runs, some groups may or may not exist, this method will not throw
		/// if the group does not exist.  But it will output a warning and create an empty group.
		const NodeGroup& get_node_group(const std::string& groupName)
		{
			auto group = node_groups_.find(groupName);
			if (group == node_groups_.end()) {
				//std::cout << "Warning! Group " + groupName + " does not exist, but was requested.  Creating an empty group." << std::endl;
				node_groups_[groupName] = NodeGroup();
				return node_groups_[groupName];
			}
			return group->second;
		}
		/// Returns the node group names
		const std::vector<std::string> get_node_group_names()
		{
			std::vector<std::string> names;
			for (auto pair : node_groups_)
				names.push_back(pair.first);
			return names;
		}
		/// Returns the iterator of the full map of groups, which is usually useful for I/O
		NodeGroupMapIt node_groups_begin() {
			return node_groups_.begin();
		};
		NodeGroupMapIt node_groups_end() {
			return node_groups_.end();
		};
		/// Clears all node groups
		void clear_node_groups() {
			node_groups_.clear();
		}

		// Element Groups

		/// Create a new element group.  If it already exists, it does nothing.  Returns if one was created.
		bool create_element_group(const std::string& groupName)
		{
			auto group = element_groups_.find(groupName);
			if (group != element_groups_.end())
				return false;
			element_groups_[groupName] = std::set<ElementPtr>();
			return true;
		}
		/// Gets rid of an element group.  Returns true if there was a group that was deleted
		bool remove_element_group(const std::string& groupName)
		{
			auto group = element_groups_.find(groupName);
			if (group == element_groups_.end())
				return false;
			element_groups_.erase(group);
			return true;
		}
		/// Adds a single element to a group.  Will throw if the group doesn't exist.
		void add_element_to_group(const std::string& groupName, ElementPtr e)
		{
			auto group = element_groups_.find(groupName);
			if (group == element_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			group->second.insert(e);
		}
		/// Adds a set of elements to a group.  Will throw if the group doesn't exist.
		void add_elements_to_group(const std::string& groupName, const std::set<ElementPtr>& eSet)
		{
			auto group = element_groups_.find(groupName);
			if (group == element_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			for (auto e : eSet)
				group->second.insert(e);
		}
		/// Adds all the elements in the mesh to a group.  Will throw if the group doesn't exist.
		void add_all_elements_to_group(const std::string& groupName)
		{
			auto group = element_groups_.find(groupName);
			if (group == element_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			for (auto& e : SuperClass::elements_)
				group->second.insert(e.get());
		}
		/// Gets rid of an element group.  Will throw if the group doesn't exist.
		void remove_element_from_group(const std::string& groupName, ElementPtr e)
		{
			auto group = element_groups_.find(groupName);
			if (group == element_groups_.end())
				throw _UTIL_MSG_EXCEPT("Group " + groupName + " does not exist.");
			group->second.erase(e);
		}
		/// Returns the element by const reference.  It can be looped over and the elements can be modified
		/// but the set itself cannot be modified.
		const ElementGroup& get_element_group(const std::string& groupName)
		{
			auto group = element_groups_.find(groupName);
			if (group == element_groups_.end()) {
				//std::cout << "Warning! Group " + groupName + " does not exist, but was requested.  Creating an empty group." << std::endl;
				element_groups_[groupName] = std::set<ElementPtr>();
				return element_groups_[groupName];
			}
			return group->second;
		}
		/// Returns the node group names
		const std::vector<std::string> get_element_group_names()
		{
			std::vector<std::string> names;
			for (auto pair : element_groups_)
				names.push_back(pair.first);
			return names;
		}
		/// Returns the iterator of the full map of groups, which is usually useful for I/O
		ElementGroupMapIt element_groups_begin() {
			return element_groups_.begin();
		};
		ElementGroupMapIt element_groups_end() {
			return element_groups_.end();
		};

		void clear_element_groups() {
			element_groups_.clear();
		}

		//////////////////
		// Misc methods //
		//////////////////

		/// Overload the PrimitiveMesh implementation to also copy the group information
		void add_to_self(const Mesh &Other);
        /// Overload the PrimitiveMesh implementation to also copy the group information
        void move_to_self(Mesh & OtherMesh);
		// Destroys the current mesh and copies the other mesh.
		Mesh<NodeType_, ElementType_>& operator=(const Mesh<NodeType_, ElementType_>& other) {
			this->clear();
			this->add_to_self(other);
			return *this;
		};

		/// Add another to this mesh, making a deepcopy of nodes and elements
		/// Also add any new nodes and elements to groups
		void add_to_self_with_groups(
			const Mesh & OtherMesh,
			std::string nodeGroupName,
			std::string elementGroupName);

        /// Move another mesh to this mesh
        /// \important Note: the other mesh will be empty after this
        void move_to_self_with_groups(
            Mesh & OtherMesh,
            std::string nodeGroupName,
            std::string elementGroupName);
	};

	template<typename NodeType_, typename ElementType_>
	inline void Mesh<NodeType_, ElementType_>::add_to_self(const Mesh & OtherMesh)
	{
        if (parents_up_to_date_) parents_up_to_date_ = false;

		// Call base class implementation
		Mesh<NodeType_, ElementType_>::PrimativeMesh::add_to_self(OtherMesh);
		// TODO transfer over group information
	}

    template<typename NodeType_, typename ElementType_>
    inline void Mesh<NodeType_, ElementType_>::move_to_self(Mesh & OtherMesh)
    {
        if (parents_up_to_date_) parents_up_to_date_ = false;

        // Call base class implementation
        Mesh<NodeType_, ElementType_>::PrimativeMesh::move_to_self(OtherMesh);
        // TODO transfer over group information
    }

	template<typename NodeType_, typename ElementType_>
	inline void Mesh<NodeType_, ElementType_>::add_to_self_with_groups(const Mesh & OtherMesh, std::string nodeGroupName, std::string elementGroupName)
	{
        if (parents_up_to_date_) parents_up_to_date_ = false;

		create_node_group(nodeGroupName);
		create_element_group(elementGroupName);

		std::map<NodePtr, NodePtr> OtherToThisNodeMap;
		for (auto OtherMeshIter = OtherMesh.nodes_begin(); OtherMeshIter != OtherMesh.nodes_end(); ++OtherMeshIter) {
			auto newNode = add_node(std::make_unique<NodeType>(**OtherMeshIter));
			add_node_to_group(nodeGroupName, newNode);
			OtherToThisNodeMap[OtherMeshIter->get()] = newNode;
		}
		for (auto OtherMeshIter = OtherMesh.elements_begin(); OtherMeshIter != OtherMesh.elements_end(); ++OtherMeshIter) {
			auto newElement = add_element(std::make_unique<ElementType>(**OtherMeshIter));
			add_element_to_group(elementGroupName, newElement);
			newElement->remap_element_nodes(OtherToThisNodeMap);
		}
		Mesh<NodeType_, ElementType_>::PrimativeMesh::update_global_numbering();
	}
    template<typename NodeType_, typename ElementType_>
    inline void Mesh<NodeType_, ElementType_>::move_to_self_with_groups(Mesh & OtherMesh, std::string nodeGroupName, std::string elementGroupName)
    {
        if (parents_up_to_date_) parents_up_to_date_ = false;

        create_node_group(nodeGroupName);
        create_element_group(elementGroupName);

        for (auto OtherMeshIter = OtherMesh.nodes_begin(); OtherMeshIter != OtherMesh.nodes_end(); ++OtherMeshIter) {
            auto newNode = add_node(std::move(*OtherMeshIter));
            add_node_to_group(nodeGroupName, newNode);
        }
        for (auto OtherMeshIter = OtherMesh.elements_begin(); OtherMeshIter != OtherMesh.elements_end(); ++OtherMeshIter) {
            auto newElement = add_element(std::move(*OtherMeshIter));
            add_element_to_group(elementGroupName, newElement);
        }
        Mesh<NodeType_, ElementType_>::PrimativeMesh::update_global_numbering();
        OtherMesh.clear();
    }
};
*/