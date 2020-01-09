#pragma once

#include "flex_mesh/common/Defines.h"
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
        using Iterator = typename std::vector<std::unique_ptr<T>>::iterator;
        using ConstIterator = typename std::vector<std::unique_ptr<T>>::const_iterator;

        MeshEntityContainer() {};
        MeshEntityContainer(MeshEntityContainer&& other) noexcept = default;
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
        auto end() { return objects_.end(); }
        /// Returns the end iterator of the container
        auto end() const { return objects_.end(); }
        /// Returns a raw pointer to the first object added
        T * front() const { return  !objects_.empty() ? objects_.front().get() : nullptr; }
        /// Returns a raw pointer to the last object added
        T * back() const { return !objects_.empty() ? objects_.back().get() : nullptr; }
        /// Access a raw pointer to an object by index
        T * operator[] (std::size_t index) const
        {
            return objects_[index].get();
        }
        /// Access a raw pointer to an object by index
        T * operator[] (std::size_t index)
        {
            return objects_[index].get();
        }
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
            for (auto it = begin(); it != end(); ) {
                auto setIt = objects_to_remove.find(it->get());
                if (setIt != objects_to_remove.end()) {
                    objects_to_remove.erase(setIt);
                    it = remove_node(it);
                }
                else {
                    it++;
                }
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
        using ElementContainer = MeshEntityContainer<Element_>;
        using NodeGroup = std::set<Node_ *>;
        using NodeGroupMap = std::map<std::string, NodeGroup>;
        using NodeGroupMapIt = typename NodeGroupMap::iterator;
        using ElementGroup = std::set<Element_>;
        using ElementGroupMap = std::map<std::string, ElementGroup>;
        using ElementGroupMapIt = typename std::map<std::string, ElementGroup>::iterator;

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
        Mesh& operator=(Mesh&&) noexcept;

        /// Clears the mesh object of its contents
        void clear();

        /////////////////////
        // Node operations //
        /////////////////////

        /// Access element container by non-const reference
        MeshEntityContainer<Node_>& nodes();
        /// Access element container by const reference
        const MeshEntityContainer<Node_>& nodes() const;
        /// Add a node via the coordinate and return a non-owning pointer to it
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
        MeshEntityContainer<Element_>& elements();
        /// Access element container by const reference
        const MeshEntityContainer<Element_>& elements() const;
        /// Returns the set of component IDs for all the elements
        std::set<ComponentIndex> get_component_ids() const;
        /// Returns a vector of element pointers that belong to the component
        /// Note: only includes elements in this mesh object
        std::vector<Element_ *> get_component(ComponentIndex component_id) const;

        //////////////////////
        // Grouping methods //
        //////////////////////

        // Node Groups

        /// Create a new node group.  If it already exists, then it will do nothing.  Returns if one is created.
        bool create_node_group(const std::string& group_name);
        /// Gets rid of an node group.  Returns true if there was a group that was deleted
        bool remove_node_group(const std::string& group_name);
        /// Adds a single node to a group.  Will throw if the group doesn't exist.
        void add_node_to_group(const std::string& group_name, Node_ * node_ptr);
        /// Adds a set of nodes to a group.  Will throw if the group doesn't exist.
        void add_nodes_to_group(const std::string& group_name, std::set<Node_ *> node_set);
        /// Adds all the nodes in the mesh to a group.  Will throw if the group doesn't exist.
        void add_all_nodes_to_group(const std::string& group_name);
        /// Gets rid of an node group.  Will throw if the group doesn't exist.
        void remove_node_from_group(const std::string& group_name, Node_ * node_ptr);
        /// Returns the node by const reference.  It can be looped over and the nodes can be modified
        /// but the set itself cannot be modified.
        /// Because in distributed runs, some groups may or may not exist, this method will not throw
        /// if the group does not exist.  But it will output a warning and create an empty group.
        const NodeGroup& get_node_group(const std::string& group_name);
        /// Returns the node group names
        const std::vector<std::string> get_node_group_names();
        /// Returns the iterator of the full map of groups, which is usually useful for I/O
        NodeGroupMapIt node_groups_begin();
        NodeGroupMapIt node_groups_end();
        /// Clears all node groups
        void clear_node_groups();

        // Element Groups

        /// Create a new element group.  If it already exists, it does nothing.  Returns if one was created.
        bool create_element_group(const std::string& group_name);
        /// Gets rid of an element group.  Returns true if there was a group that was deleted
        bool remove_element_group(const std::string& group_name);
        /// Adds a single element to a group.  Will throw if the group doesn't exist.
        void add_element_to_group(const std::string& group_name, Element_ * e_ptr);
        /// Adds a set of elements to a group.  Will throw if the group doesn't exist.
        void add_elements_to_group(const std::string& group_name, const std::set<Element_ *>& e_ptr_set);
        /// Adds all the elements in the mesh to a group.  Will throw if the group doesn't exist.
        void add_all_elements_to_group(const std::string& group_name);
        /// Gets rid of an element group.  Will throw if the group doesn't exist.
        void remove_element_from_group(const std::string& group_name, Element_ * e_ptr);
        /// Returns the element by const reference.  It can be looped over and the elements can be modified
        /// but the set itself cannot be modified.
        const ElementGroup& get_element_group(const std::string& group_name);
        /// Returns the node group names
        const std::vector<std::string> get_element_group_names();
        /// Returns the iterator of the full map of groups, which is usually useful for I/O
        ElementGroupMapIt element_groups_begin();
        ElementGroupMapIt element_groups_end();

        void clear_element_groups();

        ///////////
        // Misc. //
        ///////////

        /// This method will go through the nodes and elements and update
        /// the global numbering appropriately.  This should be called after
        /// an operation of adding nodes or elements to the mesh.  It is not
        /// called automatically for performance reasons.
        /// \param set_local_number_to_global - indicates if the local node/element numbers
        ///     should be set to the same value as the global numbers.
        void update_global_numbering(const bool set_local_number_to_global = true);

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

    protected:
        // Containers
        MeshEntityContainer<Node_> nodes_;
        MeshEntityContainer<Element_> elements_;
        // TODO the type should be PolyhedronElement<Element_> but PolyhedronElement needs
        // to be changed to take the base class as the template.
        //MeshEntityContainer<Element_> polyhedral_elements_;

        // Groups
        std::map<std::string, NodeGroup> node_groups_;
        std::map<std::string, ElementGroup> element_groups_;
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
            this->nodes_.add(std::move(new_node));
            other_to_new_node_map[other_n_ptr.get()] = nodes_.back();
        }
        for (auto& other_e_ptr : other.elements_) {
            auto new_element_unique = std::make_unique<Element_>(*other_e_ptr);
            auto new_element = this->elements_.add(std::move(new_element_unique));
            new_element->remap_element_nodes(other_to_new_node_map);
        }
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::move_to_self(Mesh & other)
    {
        mark_nodes_as_changed();
        mark_elements_as_changed();

        std::move(other.nodes_.objects_.begin(), other.nodes_.objects_.end(), std::back_inserter(this->nodes_.objects_));
        std::move(other.elements_.objects_.begin(), other.elements_.objects_.end(), std::back_inserter(this->elements_.objects_));
        other.clear();
    }

    template<typename Node_, typename Element_>
    inline Mesh<Node_, Element_>& Mesh<Node_, Element_>::operator+=(const Mesh & other) noexcept
    {
        this->add_to_self(other);
        return *this;
    }

    template<typename Node_, typename Element_>
    inline Mesh<Node_, Element_>& Mesh<Node_, Element_>::operator=(const Mesh & other) noexcept
    {
        this->clear();
        this->add_to_self(other);
        this->mark_nodes_as_changed();
        this->mark_elements_as_changed();
        return *this;
    }

    template<typename Node_, typename Element_>
    inline Mesh<Node_, Element_>& Mesh<Node_, Element_>::operator=(Mesh&& other) noexcept
    {
        if (this != &other)
        {
            this->move_to_self(other);
        }
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
    inline bool Mesh<Node_, Element_>::create_node_group(const std::string & group_name)
    {
        auto group = node_groups_.find(group_name);
        if (group != node_groups_.end()) {
            return false;
        }
        node_groups_[group_name] = NodeGroup();
        return true;
    }

    template<typename Node_, typename Element_>
    inline bool Mesh<Node_, Element_>::remove_node_group(const std::string & group_name)
    {
        auto group = node_groups_.find(group_name);
        if (group == node_groups_.end()) {
            return false;
        }
        node_groups_.erase(group);
        return true;
    }


    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::add_node_to_group(const std::string & group_name, Node_ * node_ptr)
    {
        auto group = node_groups_.find(group_name);
        if (group == node_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        group->second.insert(node_ptr);
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::add_nodes_to_group(const std::string & group_name, std::set<Node_*> node_set)
    {
        auto group = node_groups_.find(group_name);
        if (group == node_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        for (auto n : node_set) {
            group->second.insert(n);
        }
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::add_all_nodes_to_group(const std::string & group_name)
    {
        auto group = node_groups_.find(group_name);
        if (group == node_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        for (auto& n : nodes_) {
            group->second.insert(n.get());
        }
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::remove_node_from_group(const std::string & group_name, Node_ * node_ptr)
    {
        auto group = node_groups_.find(group_name);
        if (group == node_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        group->second.erase(node_ptr);
    }

    template<typename Node_, typename Element_>
    inline const typename Mesh<Node_, Element_>::NodeGroup & Mesh<Node_, Element_>::get_node_group(const std::string & group_name)
    {
        auto group = node_groups_.find(group_name);
        if (group == node_groups_.end()) {
            //std::cout << "Warning! Group " + group_name + " does not exist, but was requested.  Creating an empty group." << std::endl;
            node_groups_[group_name] = NodeGroup();
            return node_groups_[group_name];
        }
        return group->second;
    }

    template<typename Node_, typename Element_>
    inline const std::vector<std::string> Mesh<Node_, Element_>::get_node_group_names()
    {
        std::vector<std::string> names;
        for (auto pair : node_groups_)
            names.push_back(pair.first);
        return names;
    }

    template<typename Node_, typename Element_>
    inline typename Mesh<Node_, Element_>::NodeGroupMapIt Mesh<Node_, Element_>::node_groups_begin()
    {
        return node_groups_.begin();
    }

    template<typename Node_, typename Element_>
    inline typename Mesh<Node_, Element_>::NodeGroupMapIt Mesh<Node_, Element_>::node_groups_end()
    {
        return node_groups_.end();
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::clear_node_groups()
    {
        node_groups_.clear();
    }

    template<typename Node_, typename Element_>
    inline bool Mesh<Node_, Element_>::create_element_group(const std::string & group_name)
    {
        auto group = element_groups_.find(group_name);
        if (group != element_groups_.end()) {
            return false;
        }
        element_groups_[group_name] = ElementGroup();
        return true;
    }

    template<typename Node_, typename Element_>
    inline bool Mesh<Node_, Element_>::remove_element_group(const std::string & group_name)
    {
        auto group = element_groups_.find(group_name);
        if (group == element_groups_.end()) {
            return false;
        }
        element_groups_.erase(group);
        return true;
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::add_element_to_group(const std::string & group_name, Element_ * e_ptr)
    {
        auto group = element_groups_.find(group_name);
        if (group == element_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        group->second.insert(e_ptr);
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::add_elements_to_group(const std::string & group_name, const std::set<Element_*>& e_ptr_set)
    {
        auto group = element_groups_.find(group_name);
        if (group == element_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        for (auto e_ptr : e_ptr_set) {
            group->second.insert(e_ptr);
        }
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::add_all_elements_to_group(const std::string & group_name)
    {
        auto group = element_groups_.find(group_name);
        if (group == element_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        for (auto& e : elements_) {
            group->second.insert(e.get());
        }
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::remove_element_from_group(const std::string & group_name, Element_ * e_ptr)
    {
        auto group = element_groups_.find(group_name);
        if (group == element_groups_.end()) {
            throw _UTIL_MSG_EXCEPT("Group " + group_name + " does not exist.");
        }
        group->second.erase(e_ptr);
    }

    template<typename Node_, typename Element_>
    inline const typename Mesh<Node_, Element_>::ElementGroup & Mesh<Node_, Element_>::get_element_group(const std::string & group_name)
    {
        auto group = element_groups_.find(group_name);
        if (group == element_groups_.end()) {
            //std::cout << "Warning! Group " + group_name + " does not exist, but was requested.  Creating an empty group." << std::endl;
            element_groups_[group_name] = ElementGroup();
            return element_groups_[group_name];
        }
        return group->second;
    }

    template<typename Node_, typename Element_>
    inline const std::vector<std::string> Mesh<Node_, Element_>::get_element_group_names()
    {
        std::vector<std::string> names;
        for (auto pair : element_groups_)
            names.push_back(pair.first);
        return names;
    }

    template<typename Node_, typename Element_>
    inline typename Mesh<Node_, Element_>::ElementGroupMapIt Mesh<Node_, Element_>::element_groups_begin()
    {
        return element_groups_.begin();
    }

    template<typename Node_, typename Element_>
    inline typename Mesh<Node_, Element_>::ElementGroupMapIt Mesh<Node_, Element_>::element_groups_end()
    {
        return element_groups_.end();
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::clear_element_groups()
    {
        element_groups_.clear();
    }

    template<typename Node_, typename Element_>
    inline void Mesh<Node_, Element_>::update_global_numbering(const bool set_local_number_to_global)
    {
        GlobalIndex current_num = 0;
        for (auto& n_ptr : nodes_) {
            n_ptr->set_global_number(current_num);
            if (set_local_number_to_global) {
                n_ptr->set_local_number(current_num);
            }
            current_num++;
        }
        current_num = 0;
        for (auto& e_ptr : elements_) {
            e_ptr->set_global_number(current_num);
            if (set_local_number_to_global) {
                e_ptr->set_local_number(current_num);
            }
            current_num++;
        }
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
