#include "mesh/relations/ElementToComponentMap.h"

namespace mesh
{
ElementToComponentMap::ElementToComponentMap()
  : mesh_ptr_(nullptr)
{
}

void ElementToComponentMap::initialize(const Mesh& mesh)
{
    mesh_ptr_ = &mesh;

    subcomponent_to_component_map_ = { 0 };
    element_to_subcomponent_map_.resize(mesh.elements().size());
    for (auto& i : element_to_subcomponent_map_)
    {
        i = 0;
    }
    element_to_index_in_subcomponent_map_.resize(mesh.elements().size());
    std::iota(element_to_index_in_subcomponent_map_.begin(),
              element_to_index_in_subcomponent_map_.end(),
              0);

    build_maps_from_component_assignments();
}

void ElementToComponentMap::initialize(
  const Mesh& mesh,
  const std::vector<std::pair<std::string, ComponentIndex>>& group_component_assignments)
{
    mesh_ptr_ = &mesh;

    element_to_subcomponent_map_.resize(mesh.elements().size());
    for (auto& i : element_to_subcomponent_map_)
    {
        i = 0;
    }

    // Update element_to_subcomponent_map_ and subcomponent_to_component_map_
    std::map<ComponentIndex, ComponentIndex> component_to_subcomponent_map;
    subcomponent_to_component_map_.clear();
    for (auto& p : group_component_assignments)
    {
        ComponentIndex subcomponent_index;
        auto find_iter = component_to_subcomponent_map.find(p.second);
        if (find_iter == component_to_subcomponent_map.end())
        {
            subcomponent_index = static_cast<ComponentIndex>(component_to_subcomponent_map.size());
            component_to_subcomponent_map[p.second] = subcomponent_index;
            subcomponent_to_component_map_.push_back(p.second);
        }
        else
        {
            subcomponent_index = find_iter->second;
        }

        auto& group = mesh.elements().group(p.first);
        for (auto& element_index : group.element_indices)
        {
            element_to_subcomponent_map_[element_index] = subcomponent_index;
        }
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

utility::Slice<const LocalIndex> ElementToComponentMap::get_element_subcomponent_connectivity(
  const LocalIndex element_index) const
{
    // See initialize(), but element_subcomponent_connectivities_ is only filled if the connectivity
    // differs from that in the mesh.
    if (element_subcomponent_connectivities_.size())
    {
        return mesh_ptr_->elements()[element_index].connectivity();
    }
    else
    {
        auto component_info      = get_element_component_info(element_index);
        auto& element_conn_start = connectivity_element_offsets_[element_index];
        auto element_conn_end =
          element_conn_start + mesh_ptr_->elements()[element_index].connectivity().size();
        return utility::Slice<const LocalIndex>(
          &element_subcomponent_connectivities_[element_conn_start],
          element_conn_end - element_conn_start);
    }
}

LocalIndex ElementToComponentMap::get_n_nodes_in_subcomponent(const ComponentIndex& component) const
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
                    for (auto& node_index : mesh_ptr_->elements()[element_index].connectivity())
                    {
                        element_subcomponent_connectivities_.push_back(
                          local_index_to_subcomponent_index[node_index]);
                    }

                    connectivity_element_offsets_[element_index] =
                      element_subcomponent_connectivities_.size();
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