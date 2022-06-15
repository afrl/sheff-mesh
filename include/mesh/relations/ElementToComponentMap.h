#pragma once

#include "mesh/Mesh.h"
#include "mesh/element/ElementSubrangeIterator.h"

#include <utility>
#include <vector>

namespace mesh
{
struct ElementComponentInfo
{
    ComponentIndex component;
    ComponentIndex subcomponent;
    LocalIndex index_in_subcomponent;
};

class ElementToComponentMap
{
  public:
    ElementToComponentMap();

    /// @brief Initializes the maps
    void update(const Mesh& mesh);

    /// @brief Returns the list of indices of elements are part of the given component.
    /// The elements essentially compose the subcomponent corresponding to the given component in
    /// this partition.  The returned list is ordered by the elements' indices in the subcomponent.
    std::vector<LocalIndex> get_element_indices(const ComponentIndex& component) const;

    /// @brief Returns a range that iterates over the elements that are part of the given componet.
    /// An empty range will be returned if no elements exist in the given component.
    IteratorRange<ElementSubrangeConstIterator> get_elements(const ComponentIndex& component) const;

    /// @brief Returns a bool indicating if the component exists in the partition and the
    /// subcomponent index of the given component.
    ///
    /// The component ID is a global identifier for a component, but some partitions may not contain
    /// a component.  A subcomponent is a piece of a component that appears within a partition, and
    /// the subcomponent index is its index in the list of components that appear within the
    /// partition.
    std::pair<bool, ComponentIndex> get_subcomponent_index(const ComponentIndex& component) const;

    /// @brief Returns the component identifier corresponding to the given subcomponent
    ComponentIndex get_component(const ComponentIndex& subcomponent_index) const;

    /// @brief Returns the number of subcomponents in this partition
    ComponentIndex get_n_subcomponents() const;

    /// @brief Given an element, returns a struct containing the component ID, subcomponent index,
    /// and the index of the element in the subcomponent.
    ElementComponentInfo get_element_component_info(const LocalIndex element_index) const;

    /// @brief Given an element index, this returns the connectivity of for the element in
    /// terms of how the subcomponent numbers the nodes.
    utility::Slice<const LocalIndex, ConnIndex> get_element_subcomponent_connectivity(
      const LocalIndex element_index) const;

    /// Returns the number of nodes in the subcomponent, given the component id
    LocalIndex get_n_nodes_in_component(const ComponentIndex& component) const;

    /// Returns the number of nodes in the subcomponent (fast access)
    LocalIndex get_n_nodes_in_subcomponent(const ComponentIndex& subcomponent_index) const;

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("n_subcomponents", subcomponent_to_component_map_.size()),
          utility::named_value("subcomponent_to_component_map",
                               utility::slice_c(subcomponent_to_component_map_)),
          utility::named_value("element_to_subcomponent_map",
                               utility::slice_c(element_to_subcomponent_map_)));
    }

    /// @brief Sets the mesh without building anything maps (only useful before a deserialize)
    void set_mesh(const Mesh& mesh);

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        _UTIL_ASSERT(mesh_ptr_ != nullptr,
                     "Deserialize of ComponentMap requires setting the mesh (via set_mesh()) "
                     "before data is deserialized.");

        std::size_t n_subcomponents;
        s(n_subcomponents);
        subcomponent_to_component_map_.resize(n_subcomponents);
        element_to_subcomponent_map_.resize(mesh_ptr_->elements().size());
        s(utility::Slice<ComponentIndex>(subcomponent_to_component_map_.data(),
                                         subcomponent_to_component_map_.size()),
          utility::Slice<ComponentIndex>(element_to_subcomponent_map_.data(),
                                         element_to_subcomponent_map_.size()));
    }

  private:
    /// Once subcomponent_to_component_map_ and element_to_subcomponent_map_,
    /// this method builds the rest of the data.
    void build_maps_from_component_assignments();

    /// Returns the list of elements in the given component ordered by their local index.
    /// The order in which they appear in the returned list should match their index in the
    /// subcomponent as given by element_to_index_in_subcomponent_.
    std::vector<LocalIndex> extract_elements_in_component(const ComponentIndex& component) const;

    /// Returns the set of nodes in the given subcomponent.
    std::set<LocalIndex> extract_nodes_in_subcomponent(
      const ComponentIndex& subcomponent_index) const;

    // Important note: essential piece of data
    // We need to know the subcomponent index for each component and the subcomponent numbering for
    // each element.  Size: # of subcomponents
    std::vector<ComponentIndex> subcomponent_to_component_map_;
    // A map from the element local index to the subcomponent index that contains the element.
    // Important note: essential piece of data
    // Size: # of elements in mesh
    std::vector<ComponentIndex> element_to_subcomponent_map_;
    // The rest of the properties can be recovered from the above two

    // Size: # of elements in mesh
    std::vector<LocalIndex> element_to_index_in_subcomponent_map_;

    // Caches the number of nodes in each subcomponent to avoid reconstructing the set in
    // extract_nodes_in_subcomponent() whenever only the size is needed.
    std::vector<LocalIndex> n_nodes_per_subcomponent_;

    // Important note: only filled if the connectivity of the subcomponents differs from the mesh.
    // The relationship between a node and a subcomponent is not unique.  A single node can appear
    // in more than one component if it lies on the boundary of several components.  Consequently,
    // this map stores the connectivity of each element in terms of the subcomponent numbering.
    // It is ordered by subcomponent then by element in the subcomponent (IT IS NOT ORDERED BY
    // ELEMENT #).
    std::vector<LocalIndex> element_subcomponent_connectivities_;
    // Important note: only filled if the connectivity of the subcomponents differs from the mesh.
    // Offsets in element_subcomponent_connectivities_ for each element.  Similar to how
    // connectivities are stored in the mesh::Mesh class.  The end of a connectivity can deduced
    // from the connectivity size of the element in the mesh (since the size will match between the
    // mesh and subcomponent numbering).  Size: # elements + 1
    std::vector<std::size_t> connectivity_element_offsets_;

    /// A non-owning pointer that is null until update is called
    /// Needed for get_elements() to return the element accessors
    const Mesh* mesh_ptr_;
};
}