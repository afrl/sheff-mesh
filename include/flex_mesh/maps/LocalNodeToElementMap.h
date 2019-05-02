#pragma once

#include "flex_mesh/common/Defines.h"
#include "utility/DynamicArray.h"
#include "utility/VectorSet.h"

namespace flex_mesh
{
    /// Stores a map from the local index number to the local indices of the elements that
    /// contain the respective node.
    template<typename MeshTypename>
    class LocalNodeToElementMap
    {
    public:
        /// Update the map given a mesh
        void update(const MeshTypename & mesh);
        /// Returns the list of indices of parent elements given a local index for a node
        const std::vector<LocalIndex>& get_element_indices(const LocalIndex node_index) const;

    private:
        std::vector<std::vector<LocalIndex>> parent_element_indices_;
    };

    /////////////////
    // Definitions //
    /////////////////

    template<typename MeshTypename>
    inline void LocalNodeToElementMap<MeshTypename>::update(const MeshTypename & mesh)
    {
        parent_element_indices_.resize(mesh.nodes().size());

        for (auto& v : parent_element_indices_)
        {
            v.clear();
        }

        for (auto& e_ptr : mesh.elements())
        {
            for (auto n : *e_ptr)
            {
                parent_element_indices_[n->get_local_number()].push_back(e_ptr->get_local_number());
            }
        }
    }

    template<typename MeshTypename>
    inline const std::vector<LocalIndex>& LocalNodeToElementMap<MeshTypename>::get_element_indices(const LocalIndex node_index) const
    {
        return parent_element_indices_[node_index];
    }
}