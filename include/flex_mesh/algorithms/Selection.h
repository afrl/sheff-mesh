#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include <utility>
#include <vector>

namespace flex_mesh
{
    /// \defgroup Selecting Node/Elements in a Mesh
    /// @{

    /// \brief A class defining a selection of nodes/elements from the mesh.
    /// 
    /// A subset of nodes or elements are needed in a finite element analysis for:
    /// 1. Assigning boundary conditions to a subset of nodes
    /// 2. Assigning constitutive relations to a subset of elements
    /// 3. Assigning internal state variables for a subset of elements
    /// However, there is the issue that when the user wants to make these selections, the
    /// mesh may not be initialized, the degrees of freedom may not be enumerated, or the element
    /// formulations may not have been created yet.  Consequently, this class stores lists of
    /// criteria that must be satisfied, and the selection is not actually done until all the
    /// relevant data is known.  Several methods in the fea_analysis::FeaAnalysis class that
    /// takes a node or element selection as an argument. \see fea_analysis::FeaAnalysis
    ///
    /// A Selection can have multiple of two types of criteria: a LocalCriteria or NonLocalCriteria.
    /// A LocalCriteria is a criteria that is based on only the data available to the node or element.
    /// A NonLocalCriteria is a criteria that is based on data available to the whole mesh.  Both of
    /// these types of criteria should be lambda functions (or std::function objects).  To see how
    /// these criteria should be created, see the documentation for the Selection::add_criteria methods.
    template<typename MeshType, typename ObjectType>
    struct Selection
    {
    public:
        using LocalCriteria = std::function<bool(const ObjectType& obj)>;
        using NonLocalCriteria = std::function<void(const MeshType& mesh, std::set<ObjectType *>& object_set)>;

        //////////////////
        // Constructors //
        //////////////////

        Selection() {}
        /// Copy constructor takes the criteria from another Selection object, but the selected objects
        /// will be empty, requiring a call to Selection::run_selection.
        Selection(const Selection& other);
        /// Move constructor takes the criteria AND selected objects from another Selection object.
        Selection(Selection&& other) noexcept;
        virtual ~Selection() {}

        /// Add a local criteria to the Selection object.
        /// For an ElementSelection, a criteria selecting elements with component ID of 0 would be:
        /// \code{cpp}
        /// ElementSelection selection;
        /// selection.add_criteria([](const flex_mesh::Element<flex_mesh::Node>& element) -> bool {
        ///		return element.get_component_id() == 0;
        /// });
        /// \endcode
        /// For a NodeSelection, a criteria selecting nodes with an x value <= 0 would be:
        /// \code{cpp}
        /// NodeSelection selection;
        /// selection.add_criteria([](const flex_mesh::Node& node) -> bool {
        ///		return node.x() <= 0;
        /// });
        /// \endcode
        void add_criteria(const LocalCriteria& criteria);
        /// Add a non-local criteria to the Selection object.  The lambda for a NonLocalCriteria
        /// should take const-reference to a mesh and a reference to a set of objects.  Given a
        /// mesh object, it should add any nodes or elements that fit the criteria to the set and
        /// return void.
        /// For example, selecting nodes in a particular node group would look like:
        /// \code{cpp}
        /// ElementSelection selection;
        /// selection.add_criteria(
        ///     [](const MeshType& mesh, std::set<flex_mesh::Node *>& object_set) -> void {
        ///     for (auto n_ptr : mesh.get_node_group("SomeNodeGroup"))
        ///         object_set.insert(n_ptr);
        /// });
        /// \endcode
        void add_criteria(const NonLocalCriteria& criteria);

        /// Run the selection criteria over a mesh, which will fill the selected_objects_ property.
        /// The selected objects can be accessed const-ref via Selection::get_selection().
        virtual void run_selection(const MeshType& mesh) = 0;
        /// Accesses the selected objects const-ref.  Note: this will return an empty set until
        /// Selection::run_selection() is called.
        const std::set<ObjectType *>& get_selection() const;

    protected:
        std::set<ObjectType *> selected_objects_;
        std::vector<LocalCriteria> local_criteria_;
        std::vector<NonLocalCriteria> non_local_criteria_;
    };

    /// A class defining a selection of elements. \see Selection
    template<typename MeshTypename, typename ElementTypename>
    struct ElementSelection : public Selection<MeshTypename, ElementTypename>
    {
    public:
        virtual ~ElementSelection() {}

        /// A specialization of Selection::run_selection for elements.
        virtual void run_selection(const MeshTypename& mesh) override;
    };

    /// A class defining a selection of nodes. \see Selection
    template<typename MeshTypename, typename NodeTypename>
    struct NodeSelection : public Selection<MeshTypename, NodeTypename>
    {
    public:
        virtual ~NodeSelection() {}

        /// A specialization of Selection::run_selection for nodes.
        virtual void run_selection(const MeshTypename& mesh) override;
    };

    /////////////////
    // Definitions //
    /////////////////

    template<typename MeshType, typename ObjectType>
    inline Selection<MeshType, ObjectType>::Selection(const Selection & other)
        : local_criteria_(other.local_criteria_), non_local_criteria_(other.non_local_criteria_)
    {
    }

    template<typename MeshType, typename ObjectType>
    inline Selection<MeshType, ObjectType>::Selection(Selection && other) noexcept
    {
        selected_objects_.swap(other.selected_objects_);
        local_criteria_.swap(other.local_criteria_);
        non_local_criteria_.swap(other.non_local_criteria_);
    }

    template<typename MeshType, typename ObjectType>
    inline void Selection<MeshType, ObjectType>::add_criteria(const LocalCriteria & criteria)
    {
        local_criteria_.push_back(criteria);
    }

    template<typename MeshType, typename ObjectType>
    inline void Selection<MeshType, ObjectType>::add_criteria(const NonLocalCriteria & criteria)
    {
        non_local_criteria_.push_back(criteria);
    }

    template<typename MeshType, typename ObjectType>
    inline const std::set<ObjectType*>& Selection<MeshType, ObjectType>::get_selection() const
    {
        return selected_objects_;
    }

    template<typename MeshTypename, typename ElementTypename>
    inline void ElementSelection<MeshTypename, ElementTypename>::run_selection(const MeshTypename & mesh)
    {
        this->selected_objects_.clear();

        if (this->non_local_criteria_.empty())
        {
            for (auto& e_ptr : mesh.elements())
            {
                bool meets_criteria = true;
                for (auto& criterion : this->local_criteria_)
                {
                    if (!criterion(*e_ptr))
                    {
                        meets_criteria = false;
                        break;
                    }
                }
                if (meets_criteria)
                {
                    this->selected_objects_.insert(e_ptr.get());
                }
            }
        }
        else
        {
            // Get the intersection of all non-local criteria
            std::set<ElementTypename *> non_local_criteria_intersection;
            for (auto& criterion : this->non_local_criteria_)
            {
                std::set<ElementTypename *> criteria_set;
                criterion(mesh, criteria_set);

                // Take the intersection of criteria_set and the current intersection of previous
                // criteria (non_local_criteria_intersection) and put the result in tmp_intersection.
                std::set<ElementTypename *> tmp_intersection;
                std::set_intersection(
                    non_local_criteria_intersection.begin(), non_local_criteria_intersection.end(),
                    criteria_set.begin(), criteria_set.end(),
                    std::inserter(tmp_intersection, tmp_intersection.begin()));

                // We don't need tmp_intersection so just swap with non_local_criteria_intersection
                // rather than copy
                non_local_criteria_intersection.swap(tmp_intersection);
            }

            // Now run the local criteria over the intersection of the non-local criteria to get
            // the final set.
            for (auto& e_ptr : non_local_criteria_intersection)
            {
                bool meets_criteria = true;
                for (auto& criterion : this->local_criteria_)
                {
                    if (!criterion(*e_ptr))
                    {
                        meets_criteria = false;
                        break;
                    }
                }
                if (meets_criteria)
                {
                    this->selected_objects_.insert(e_ptr);
                }
            }
        }
    }

    template<typename MeshTypename, typename NodeTypename>
    inline void NodeSelection<MeshTypename, NodeTypename>::run_selection(const MeshTypename & mesh)
    {
        this->selected_objects_.clear();

        if (this->non_local_criteria_.empty())
        {
            for (auto& n_ptr : mesh.nodes())
            {
                bool meets_criteria = true;
                for (auto& criterion : this->local_criteria_)
                {
                    if (!criterion(*n_ptr))
                    {
                        meets_criteria = false;
                        break;
                    }
                }
                if (meets_criteria)
                {
                    this->selected_objects_.insert(n_ptr.get());
                }
            }
        }
        else
        {
            // Get the intersection of all non-local criteria
            std::set<NodeTypename *> non_local_criteria_intersection;
            for (auto& criterion : this->non_local_criteria_)
            {
                std::set<NodeTypename *> criteria_set;
                criterion(mesh, criteria_set);

                // Take the intersection of criteria_set and the current intersection of previous
                // criteria (non_local_criteria_intersection) and put the result in tmp_intersection.
                std::set<NodeTypename *> tmp_intersection;
                std::set_intersection(
                    non_local_criteria_intersection.begin(), non_local_criteria_intersection.end(),
                    criteria_set.begin(), criteria_set.end(),
                    std::inserter(tmp_intersection, tmp_intersection.begin()));

                // We don't need tmp_intersection so just swap with non_local_criteria_intersection
                // rather than copy
                non_local_criteria_intersection.swap(tmp_intersection);
            }

            // Now run the local criteria over the intersection of the non-local criteria to get
            // the final set.
            for (auto& n_ptr : non_local_criteria_intersection)
            {
                bool meets_criteria = true;
                for (auto& criterion : this->local_criteria_)
                {
                    if (!criterion(*n_ptr))
                    {
                        meets_criteria = false;
                        break;
                    }
                }
                if (meets_criteria)
                {
                    this->selected_objects_.insert(n_ptr);
                }
            }
        }
    }

    /// @}
}