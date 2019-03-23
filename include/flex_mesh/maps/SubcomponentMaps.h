#pragma once

#include "flex_mesh/common/Defines.h"
#include "flex_mesh/algorithms/Search.h"
#include "utility/Range.h"

#include <map>

namespace flex_mesh
{
    class SubcomponentMap
    {
    public:

        /// Note: update_maps to create the maps from partition to subcomponent numbering
        SubcomponentMap(const ComponentIndex component_id,
            const ComponentIndex subcomponent_index);

        const ComponentIndex get_component_id() const { return component_id_; }

        /// The subcomponent index of a node given the node's local index (partition-based index)
        const LocalIndex get_subcomponent_index(const LocalIndex node_local_index) const;
        const std::map<LocalIndex, LocalIndex> & get_map() const { return local_to_subcomponent_node_map_; }

        /// Update the index map given a the nodes and elements of the partition
        template<typename Node_, typename ElementContainer>
        void update_maps(const ElementContainer & elements);

    private:
        /// The component id of this subcomponent
        ComponentIndex component_id_;
        /// On the partition, this is the index of the subcomponent
        ComponentIndex subcomponent_index_;
        /// A map from the local number of a node to the index within this subcomponent.
        /// Note: only nodes that appear in an element within the component appear in this map.
        std::map<LocalIndex, LocalIndex> local_to_subcomponent_node_map_;
    };

    class SubcomponentMaps
    {
    public:
        SubcomponentMaps() {}

        const ComponentIndex get_n_subcomponents() const { return static_cast<flex_mesh::ComponentIndex>(component_id_to_subcomponent_index_.size()); }
        const ComponentIndex get_subcomponent_index(const ComponentIndex component_id) const;
        const ComponentIndex get_component_id(const LocalIndex subcomponent_index) const;
        const LocalIndex get_n_nodes_in_subcomponent(const ComponentIndex subcomponent_index) const;
        const std::map<LocalIndex, LocalIndex> & get_map(const ComponentIndex component_id) const;

        /// Update the index map given a the nodes and elements of the partition
        template<typename MeshTypename>
        void update_maps(const MeshTypename & mesh);
    private:
        std::vector<SubcomponentMap> subcomponent_maps_;
        std::map<ComponentIndex, ComponentIndex> component_id_to_subcomponent_index_;
    };

    /////////////////
    // Definitions //
    /////////////////

    inline SubcomponentMap::SubcomponentMap(const ComponentIndex component_id, const ComponentIndex subcomponent_index)
        : component_id_(component_id), subcomponent_index_(subcomponent_index)
    {
    }

    inline const LocalIndex SubcomponentMap::get_subcomponent_index(const LocalIndex node_local_index) const
    {
        return local_to_subcomponent_node_map_.find(node_local_index)->second;
    }

    template<typename Node_, typename ElementContainer>
    inline void SubcomponentMap::update_maps(const ElementContainer & elements)
    {
        local_to_subcomponent_node_map_.clear();
        auto nodes_in_subcomponent_pair = find_nodes_in_component<Node_>(elements, component_id_);
        for (auto i : utility::indices(nodes_in_subcomponent_pair.first))
        {
            local_to_subcomponent_node_map_[nodes_in_subcomponent_pair.first[i]] =
                nodes_in_subcomponent_pair.second[i]->get_local_number();
        }
    }

    inline const ComponentIndex SubcomponentMaps::get_subcomponent_index(const ComponentIndex component_id) const
    {
        return component_id_to_subcomponent_index_.find(component_id)->second;
    }

    inline const ComponentIndex SubcomponentMaps::get_component_id(const LocalIndex subcomponent_index) const
    {
        return subcomponent_maps_[subcomponent_index].get_component_id();
    }

    inline const LocalIndex SubcomponentMaps::get_n_nodes_in_subcomponent(const ComponentIndex subcomponent_index) const
    {
        return static_cast<LocalIndex>(subcomponent_maps_[subcomponent_index].get_map().size());
    }

    inline const std::map<LocalIndex, LocalIndex>& SubcomponentMaps::get_map(const ComponentIndex component_id) const
    {
        auto subcomponent_index = get_subcomponent_index(component_id);
        return subcomponent_maps_[subcomponent_index].get_map();
    }

    template<typename MeshTypename>
    inline void SubcomponentMaps::update_maps(const MeshTypename & mesh)
    {
        subcomponent_maps_.clear();
        component_id_to_subcomponent_index_.clear();

        // Find what components are in this partition
        std::set<ComponentIndex> component_ids_in_partition;
        for (auto& e_ptr : mesh.elements())
        {
            component_ids_in_partition.insert(e_ptr->get_component_id());
        }

        // Figure out what each subcomponent contains in this partition
        ComponentIndex subcomponent_index = 0;
        for (auto& component_id : component_ids_in_partition)
        {
            component_id_to_subcomponent_index_[component_id] = subcomponent_index;
            subcomponent_maps_.push_back(SubcomponentMap(component_id, subcomponent_index));
            auto elements_in_subcomponent_pair = find_elements_in_component<typename MeshTypename::ElementTypename>(mesh.elements(), component_id);
            subcomponent_maps_.back().update_maps<typename MeshTypename::NodeTypename>(elements_in_subcomponent_pair.second);
        }
    }
}