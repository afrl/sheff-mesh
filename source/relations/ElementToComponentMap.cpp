#include "..\..\include\mesh\relations\ElementToComponentMap.h"

#include "mesh/relations/ElementToComponentMap.h"

namespace mesh
{
ElementToComponentMap::ElementToComponentMap()
  : mesh_ptr_(nullptr)
{
}

void ElementToComponentMap::update(const Mesh& mesh)
{
    mesh_ptr_ = &mesh;

    std::set<ComponentIndex> component_ids_set;
    for (auto& element : mesh.elements())
    {
        component_ids_set.insert(element.component_id());
    }

    subcomponent_to_component_map_.clear();
    subcomponent_to_component_map_.resize(component_ids_set.size());
    ComponentIndex subcomponent_index = 0;
    std::map<ComponentIndex, ComponentIndex> component_to_subcomponent_map;
    for (auto& component_id : component_ids_set)
    {
        subcomponent_to_component_map_[subcomponent_index] = component_id;
        component_to_subcomponent_map[component_id]        = subcomponent_index;
        subcomponent_index++;
    }

    element_to_subcomponent_map_.resize(mesh.elements().size());
    for (auto& element : mesh.elements())
    {
        element_to_subcomponent_map_[element.local_index()] =
          component_to_subcomponent_map[element.component_id()];
    }

    build_maps_from_component_assignments();
}

std::vector<LocalIndex> ElementToComponentMap::get_element_indices(
  const ComponentIndex& component) const
{
    _UTIL_ASSERT_C(mesh_ptr_ != nullptr, "initialize() must be called first.");

    return extract_elements_in_component(component);
}

IteratorRange<ElementSubrangeConstIterator> ElementToComponentMap::get_elements(
  const ComponentIndex& component) const
{
    _UTIL_ASSERT_C(mesh_ptr_ != nullptr, "initialize() must be called first.");

    auto element_indices = extract_elements_in_component(component);
    return {
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     element_indices), // begin
        ElementSubrangeConstIterator(mesh_ptr_->nodes(),
                                     mesh_ptr_->elements(),
                                     element_indices,
                                     static_cast<LocalIndex>(element_indices.size())) // end
    };
}


std::pair<bool, ComponentIndex> ElementToComponentMap::get_subcomponent_index(
  const ComponentIndex& component) const
{
    _UTIL_ASSERT_C(mesh_ptr_ == nullptr, "initialize() must be called first.");

    for (auto subcomponent_index :
         utility::range_as<ComponentIndex>(subcomponent_to_component_map_.size()))
    {
        if (subcomponent_to_component_map_[subcomponent_index] == component)
        {
            return { true, subcomponent_index };
        }
    }
    return { false, -1 };
}

ComponentIndex ElementToComponentMap::get_component(const ComponentIndex& subcomponent_index) const
{
    return subcomponent_to_component_map_[subcomponent_index];
}

ComponentIndex ElementToComponentMap::get_n_subcomponents() const
{
    return static_cast<ComponentIndex>(subcomponent_to_component_map_.size());
}

ElementComponentInfo ElementToComponentMap::get_element_component_info(
  const LocalIndex element_index) const
{
    _UTIL_ASSERT_C(mesh_ptr_ != nullptr, "initialize() must be called first.");
    _UTIL_ASSERT_C(element_index < element_to_subcomponent_map_.size() && element_index >= 0,
                   "Element index is out of range.");

    return { subcomponent_to_component_map_[element_to_subcomponent_map_[element_index]],
             element_to_subcomponent_map_[element_index],
             element_to_index_in_subcomponent_map_[element_index] };
}

utility::Slice<const LocalIndex, ConnIndex>
ElementToComponentMap::get_element_subcomponent_connectivity(const LocalIndex element_index) const
{
    // See initialize(), but element_subcomponent_connectivities_ is only filled if the connectivity
    // differs from that in the mesh.
    if (element_subcomponent_connectivities_.empty())
    {
        return mesh_ptr_->elements()[element_index].connectivity();
    }
    else
    {
        auto component_info      = get_element_component_info(element_index);
        auto& element_conn_start = connectivity_element_offsets_[element_index];
        auto element_conn_end =
          element_conn_start + mesh_ptr_->elements()[element_index].connectivity().size();
        return utility::Slice<const LocalIndex, ConnIndex>(
          &element_subcomponent_connectivities_[element_conn_start],
          static_cast<ConnIndex>(element_conn_end - element_conn_start));
    }
}

LocalIndex ElementToComponentMap::get_n_nodes_in_component(const ComponentIndex& component) const
{
    _UTIL_ASSERT_C(mesh_ptr_ != nullptr, "initialize() must be called first.");

    for (auto subcomponent_index : utility::range(subcomponent_to_component_map_.size()))
    {
        if (subcomponent_to_component_map_[subcomponent_index] == component)
        {
            return n_nodes_per_subcomponent_[subcomponent_index];
        }
    }
    return 0;
}

LocalIndex ElementToComponentMap::get_n_nodes_in_subcomponent(
  const ComponentIndex& subcomponent_index) const
{
    _UTIL_ASSERT_C(mesh_ptr_ != nullptr, "initialize() must be called first.");
    return n_nodes_per_subcomponent_[subcomponent_index];
}

void ElementToComponentMap::set_mesh(const Mesh& mesh)
{
    mesh_ptr_ = &mesh;
}

void ElementToComponentMap::build_maps_from_component_assignments()
{

    // Update element_to_index_in_subcomponent_map_
    std::vector<LocalIndex> element_count_by_subcomponent(subcomponent_to_component_map_.size(), 0);
    element_to_index_in_subcomponent_map_.resize(mesh_ptr_->elements().size());
    for (auto element_index : utility::range(element_to_subcomponent_map_.size()))
    {
        element_to_index_in_subcomponent_map_[element_index] =
          element_count_by_subcomponent[element_to_subcomponent_map_[element_index]]++;
    }

    if (subcomponent_to_component_map_.size() == 1)
    {
        n_nodes_per_subcomponent_ = { mesh_ptr_->nodes().size() };

        // Leave element_subcomponent_connectivities_ and connectivity_element_offsets_
        // empty since the connectivity is the same as that in the mesh.
        element_subcomponent_connectivities_.clear();
        connectivity_element_offsets_.clear();
    }
    else
    {
        // Update element_subcomponent_connectivities_, connectivity_element_offsets_, and
        // n_nodes_per_subcomponent_
        element_subcomponent_connectivities_.clear();
        connectivity_element_offsets_.resize(element_to_subcomponent_map_.size() + 1, 0);
        n_nodes_per_subcomponent_.resize(subcomponent_to_component_map_.size());
        for (auto subcomponent_index :
             utility::range_as<ComponentIndex>(subcomponent_to_component_map_.size()))
        {
            // Create a map from the local index to the subcomponent index
            LocalIndex count = 0;
            std::map<LocalIndex, LocalIndex> local_index_to_subcomponent_index;
            for (auto& node_index : extract_nodes_in_subcomponent(subcomponent_index))
            {
                local_index_to_subcomponent_index[node_index] = count++;
            }
            n_nodes_per_subcomponent_[subcomponent_index] = count;

            // Add the remapped connectivites for the elements in the subcomponent
            for (auto element_index :
                 utility::range_as<LocalIndex>(element_to_subcomponent_map_.size()))
            {
                if (element_to_subcomponent_map_[element_index] == subcomponent_index)
                {

                    connectivity_element_offsets_[element_index] =
                      element_subcomponent_connectivities_.size();

                    for (auto& node_index : mesh_ptr_->elements()[element_index].connectivity())
                    {
                        element_subcomponent_connectivities_.push_back(
                          local_index_to_subcomponent_index[node_index]);
                    }
                }
            }
        }
    }
}

std::vector<LocalIndex> ElementToComponentMap::extract_elements_in_component(
  const ComponentIndex& component) const
{
    for (auto subcomponent_index : utility::range(subcomponent_to_component_map_.size()))
    {
        if (subcomponent_to_component_map_[subcomponent_index] == component)
        {
            std::vector<LocalIndex> element_indices;
            for (auto i : utility::range(element_to_subcomponent_map_.size()))
            {
                if (element_to_subcomponent_map_[i] == subcomponent_index)
                {
                    element_indices.push_back(static_cast<LocalIndex>(i));
                }
            }
            return element_indices;
        }
    }
    return std::vector<LocalIndex>();
}

std::set<LocalIndex> ElementToComponentMap::extract_nodes_in_subcomponent(
  const ComponentIndex& subcomponent_index) const
{
    std::set<LocalIndex> nodes_in_subcomponent;
    for (auto element_index : utility::range_as<LocalIndex>(element_to_subcomponent_map_.size()))
    {
        if (element_to_subcomponent_map_[element_index] == subcomponent_index)
        {
            for (auto& node_index : mesh_ptr_->elements()[element_index].connectivity())
            {
                nodes_in_subcomponent.insert(node_index);
            }
        }
    }
    return nodes_in_subcomponent;
}
}