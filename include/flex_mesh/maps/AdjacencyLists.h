#pragma once

#include "LocalNodeToElementMap.h"

namespace flex_mesh
{
    template<typename VertexIndexType>
    struct AdjacencyVertex
    {
        AdjacencyVertex() : skip_in_neighbor_flag(false) {}

        /// A flag that indicates if the this vertex should be skipped in
        /// an accumulation of neighbors search
        bool skip_neighbor_flag;
        /// The local index of the neighbors that are adjacent to this vertex
        /// Should I use vector set?
        utility::VectorSet<VertexIndexType> neighbors;
    };

    template<typename VertexType, typename VertexPtrContainer,
        typename VertexIndexType, typename MeshType, typename AdjacencySizeType>
        class AdjacencyList
    {
    public:
        using AdjacencyVertexType = AdjacencyVertex<VertexIndexType>;

        AdjacencyList() : underlying_vertex_pointers_(nullptr) {}

        virtual void initialize(const MeshType & mesh,
            const VertexPtrContainer & vertex_pointers) = 0;

        template<typename IndexContainer>
        void skip_in_neighbor_search(const IndexContainer & vertex_indices_to_skip);

        std::vector<VertexType *> get_neighbors(const VertexIndexType & index) const;

    private:

        std::vector<AdjacencyVertexType> adjacency_list_;
        const VertexPtrContainer * underlying_vertex_pointers_;
    };

    template<typename MeshType>
    class ElementAdjacencyList : public AdjacencyList<typename MeshType::ElementTypename,
        typename MeshType::ElementContainer, LocalIndex, MeshType,
        typename MeshType::ElementTypename::ConnSize>
    {
    public:
        virtual void initialize(const MeshType & mesh) override;
    };

    /////////////////
    // Definitions //
    /////////////////

    template<typename VertexType, typename VertexPtrContainer, typename VertexIndexType, typename MeshType, typename AdjacencySizeType>
    template<typename IndexContainer>
    inline void AdjacencyList<VertexType, VertexPtrContainer, VertexIndexType, MeshType, AdjacencySizeType>::skip_in_neighbor_search(const IndexContainer & vertex_indices_to_skip)
    {
        for (auto v : adjacency_list_)
        {
            v.skip_neighbor_flag = false;
        }

        for (auto i : vertex_indices_to_skip)
        {
            adjacency_list_[i].skip_neighbor_flag = true;
        }
    }

    template<typename VertexType, typename VertexPtrContainer, typename VertexIndexType, typename MeshType, typename AdjacencySizeType>
    inline std::vector<VertexType*> AdjacencyList<VertexType, VertexPtrContainer, VertexIndexType, MeshType, AdjacencySizeType>::get_neighbors(const VertexIndexType & index) const
    {
        std::vector<VertexType> neighbors;
        for (auto neighbor_index : adjacency_list_[index].neighbors)
        {
            if (!adjacency_list_[neighbor_index].skip_neighbor_flag)
            {
                // Extra indirection to avoid issue with container of unique_ptr's
                neighbors.push_back(&(*(*underlying_vertex_objects_)[neighbor_index]));
            }
        }
        return neighbors;
    }

    template<typename MeshType>
    inline void ElementAdjacencyList<MeshType>::initialize(const MeshType & mesh)
    {
        this->adjacency_list_.resize(mesh.elements().size());

        LocalNodeToElementMap<MeshType> node_parent_map;
        node_parent_map.update(mesh);

        for (auto& n_ptr : mesh.nodes())
        {
            auto& parent_indices = this->node_parent_map.get_element_indices(
                n_ptr->get_local_number());

            for (auto i : parent_indices)
            {
                for (auto j : parent_indices)
                {
                    if (i != j)
                    {
                        this->adjacency_list_[i].neighbors.insert(j);
                    }
                }
            }
        }
    }
}