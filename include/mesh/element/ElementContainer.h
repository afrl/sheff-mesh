#pragma once

#include "mesh/Defines.h"
#include "mesh/common/ElementTypes.h"
#include "mesh/common/QuaternionContainer.h"
#include "mesh/element/ElementGroup.h"
#include "mesh/node/NodeContainer.h"
#include "utility/Slice.h"
#include "utility/serialization/NamedValue.h"

#include <memory>
#include <vector>

namespace mesh
{
template<typename NC, typename EC>
class ElementAccessorBase;

template<typename ElementAccessorType>
class ElementIteratorBase;

/// @brief Stores information pertaining to a set of elements, including indices, geometric types,
/// connectivities, local coordinate systems, and data neccessary for FEM/XFEM. This class along
/// with an NodeContainer object would compose information about a mesh.
class ElementContainer
{
  public:
    using ElementAccessor      = ElementAccessorBase<NodeContainer, ElementContainer>;
    using ElementConstAccessor = ElementAccessorBase<const NodeContainer, const ElementContainer>;
    using ElementIterator      = ElementIteratorBase<ElementAccessor>;
    using ElementConstIterator = ElementIteratorBase<ElementConstAccessor>;

    /// A reference to the node container held to allow the construction the element
    /// accessors/iterators, but it is not used within this class.
    ElementContainer(NodeContainer& nodes);

    bool operator==(const ElementContainer& other) const;

    /// @name Accessors
    ///@{
    /// Enables ranged-based loop over elements.
    ElementIterator begin();
    /// Enables ranged-based loop over elements (const).
    ElementConstIterator begin() const;
    /// Enables ranged-based loop over elements.
    ElementConstIterator end() const;

    /// Access a particular element by local index
    ElementAccessor operator[](LocalIndex i);
    /// Const access a particular element by local index
    ElementConstAccessor operator[](LocalIndex i) const;

    /// Access the last element
    ElementAccessor back();
    /// Const access the last element
    ElementConstAccessor back() const;

    /// Returns the size of the container.
    LocalIndex size() const;

    /// Access a particular group by its name
    ElementGroup& group(const std::string group_name);
    /// Const access a particular group by its name
    const ElementGroup& group(const std::string group_name) const;
    /// Access the collection of node groups
    std::vector<ElementGroup>& groups();
    /// Const access the collection of node groups
    const std::vector<ElementGroup>& groups() const;
    ///@}

    /// @name Mutators
    ///@{
    /// Adds an element to the end of the container.  Assigns the global index of the new node to
    /// the next highest number (on this rank).
    void add(ElementType type, const std::vector<LocalIndex>& connectivity);
    /// Adds an element to the end of the container.  Assigns parent partition to 0.
    void add(ElementType type, const std::vector<LocalIndex>& connectivity, GlobalIndex i);
    /// Adds an element to the end of the container.
    void add(ElementType type,
             const std::vector<LocalIndex>& connectivity,
             GlobalIndex i,
             PartitionIndex parent_partition);
    /// Adds an element to the end of the container.  Assigns the global index of the new node to
    /// the next highest number (on this partition).
    void add(ElementType type, const utility::Slice<const LocalIndex>& connectivity);
    /// Adds an element to the end of the container.  Assigns parent partition to 0.
    void add(ElementType type, const utility::Slice<const LocalIndex>& connectivity, GlobalIndex i);
    /// Adds an element to the end of the container.
    void add(ElementType type,
             const utility::Slice<const LocalIndex>& connectivity,
             GlobalIndex i,
             PartitionIndex parent_partition);

    /// Adds an element that is a copy of another element.
    ElementAccessor add_copy(const ElementConstAccessor& other_element);

    /// Swaps all data with the other container, but does not change the reference to the
    /// NodeContainer, so be sure to swap the nodes as well to keep the mesh in a valid state.
    void swap(ElementContainer& other);

    /// Permutes the elements in the container to rearrange them with the given order.
    /// @param backward_permutation_indices - a map that goes from the new index to the old index of
    ///        the reordering.
    void permute(const std::vector<LocalIndex>& backward_permutation_indices);
    /// Permutes the connectivities of the elements to reflect a permutation of the NodeContainer
    /// that this object references.  Be sure to call this if NodeContainer::permute is called.
    /// @param node_forward_permutation_indices - a map that goes from the old index to the new
    ///        index of the reordering of nodes.  Note: this is the inverse of the map passed to
    ///        NodeContainer::permute.
    void permute_connectivity(const std::vector<LocalIndex>& node_forward_permutation_indices);

    /// @brief Deletes all elements from the container.
    void clear();
    ///@}


    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("n_elements_in_partition", global_indices_.size()),
          utility::named_value("global_indices", utility::slice(global_indices_)),
          utility::named_value("types", utility::slice(types_)),
          utility::named_value("n_connectivity_entries", connectivities_.size()),
          utility::named_value("connectivities", utility::slice(connectivities_)),
          utility::named_value("connectivity_offsets", utility::slice(connectivity_offsets_)),
          utility::named_value("parent_partitions", utility::slice(parent_partitions_)));

        s(utility::named_value("has_lcs", lcs_.size() > 0));
        if (lcs_.size() > 0)
        {
            s(utility::named_value("lcs", utility::slice(lcs_)));
        }
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        std::size_t n_elements;
        s(n_elements);
        global_indices_.resize(n_elements);
        types_.resize(n_elements);
        s(utility::slice(global_indices_));
        s(utility::slice(types_));

        std::size_t n_connectivity_entries;
        s(n_connectivity_entries);
        connectivities_.resize(n_connectivity_entries);
        connectivity_offsets_.resize(n_elements + 1);
        s(utility::slice(connectivities_), utility::slice(connectivity_offsets_));

        parent_partitions_.resize(n_elements);
        s(utility::slice(parent_partitions_));
        
        bool has_lcs;
        s(has_lcs);
        if (has_lcs)
        {
            lcs_.resize(n_elements);
            s(utility::slice(lcs_));
        }
    }

    friend class ElementAccessor;
    friend class ElementConstAccessor;

  private:
    /// Reference held for construction of element accessors/iterators
    NodeContainer& nodes_;

    /// @name Required members
    ///@{
    /// Stores a unique identifier for each element.  The global set of indices is a 0-based,
    /// continuous range of integers.  In a distributed analysis, two different partition may have
    /// the same element, though only one of the partitions owns the element.
    /// Expected size: # of elements on partition.
    std::vector<GlobalIndex> global_indices_;
    /// Stores the geometric type of each element.  @see ElementType
    /// Expected size: # of elements on partition.
    std::vector<ElementType> types_;
    /// Stores the connectivities of all elements on the partition as a 1D array.  The connectivity
    /// is described as an ordered list of local indices of nodes.  The vector connectivity_offsets_
    /// indicates which subrange of values corresponds to a given element.  For example, element i
    /// has the connectivity of connectivities_[connectivity_offsets_[i]] ...
    /// connectivities_[connectivity_offsets_[i+1]]. Expected size: sum of # of nodes in all
    /// elements on partition.
    std::vector<LocalIndex> connectivities_;
    /// Expected size: # of elements on partition + 1.  @see ElementContainer::connectivities_ for
    /// more details
    std::vector<std::size_t> connectivity_offsets_;
    /// Stores the parent partition of each element.  Each element is owned by exactly one
    /// partition, while all elements on a partition that are not owned are considered ghost
    /// elements.  During meshing, the parent partitions are not determined, so ths vector may be
    /// empty. Expected size: # of elements on partition.
    std::vector<PartitionIndex> parent_partitions_;
    ///@}

    /// @name Group members
    ///@{
    ///
    std::vector<ElementGroup> groups_;
    ///@}

    /// @name Optional members
    /// Descr.
    ///@{
    /// Stores the local coordinate systems for the element.  Each QuaternionContainer object
    /// may contain multiple sets of quaternions depending on the LCSType.  This is typically
    /// only used during mesh generation as it does not allow different fields within a given
    /// element to have different local coordinate systems.  During an analysis, the FieldLayout
    /// objects will contain BasalField objects (one for each basal field) that contain their
    /// own local coordinate system definition, and this vector will be empty.
    /// Expected size: # of elements on partition.
    std::vector<mesh::QuaternionContainer> lcs_;
    /// Contains information pertaining to fields and enrichments during a finite element analysis.
    /// During meshing or visualizing the geometry, this vector will be empty.  Note that each
    /// object will be a polymorphic type deriving from FieldLayout. Expected size: # of elements on
    /// partition.
    // std::vector<std::unique_ptr<FieldLayoutBase>> field_layouts_;
    ///@}
};

}