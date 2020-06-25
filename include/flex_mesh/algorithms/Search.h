#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/common/Utility.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include <functional>
#include <utility>
#include <unordered_set>
#include <vector>

namespace flex_mesh
{
    /// Finds the list of nodes and elements that are part of a given component in the given
    /// mesh.  Returns a pair of vectors.  The first is the local numbers of the nodes, and
    /// the second is the local numbers of the elements.
    template<typename MeshType>
    static std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> find_in_component(
        const MeshType& mesh,
        const ComponentIndex component_id);

    /// Similar to the definition above but searches a subset of elements.
    template<typename MeshType>
    static std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> find_in_component(
        const MeshType& mesh,
        const std::vector<LocalIndex>& subset_nodes_to_search,
        const std::vector<LocalIndex>& subset_elements_to_search,
        const ComponentIndex component_id);

    /// Returns the set of component ids that are referenced in the elements provided.
    template<typename ElementContainer>
    static std::set<ComponentIndex> get_component_ids(const ElementContainer& elements);

    /// Returns a pair with the first being the number of nodes in the component and
    /// the second being the number of elements in the component.
    template<typename Node_, typename ElementContainer>
    static std::pair<LocalIndex, LocalIndex> get_n_nodes_and_elements_in_component(
        const ElementContainer& elements,
        const ComponentIndex component_id);

    /// Returns the minimum coordinate (minimum x, minimum, y, minimum z) based
    /// on a collectino of node pointers.
    template<typename NodeContainer>
    static std::array<FloatType, 3> min_coord(const NodeContainer& nodes);

    /// Returns the maximum coordinate (maximum x, maximum, y, maximum z) based
    /// on a collectino of node pointers.
    template<typename NodeContainer>
    static std::array<FloatType, 3> max_coord(const NodeContainer& nodes);

    /// Returns the number of partitions (assuming they are contiguous 0-based when assigned)
    template<typename MeshType>
    PartitionIndex get_num_partitions(const MeshType& mesh);

    /// Sort the nodes and elements by the parent partition id (ascending)
    template<typename MeshType>
    void sort_by_partition(const MeshType& mesh);

    template<typename MeshType>
    static typename MeshType::ElementTypename* find_element_with_point(
        const MeshType& mesh,
        const simple_math::Point x);

    /// A helper struct for find_faces_from_node_set
    struct ElementFacePair
    {
        LocalIndex local_element_index;
        ConnIndex face_index;
    };

    /// Find all faces such that all nodes in the face appears in the given node set.
    template<typename MeshType, typename NodeTypename>
    static std::vector<ElementFacePair> find_faces_from_node_set(
        const MeshType& mesh,
        const std::set<NodeTypename*>& node_set);

    /////////////////
    // Definitions //
    /////////////////

    template<typename MeshType>
    std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> find_in_component(
        const MeshType& mesh,
        const ComponentIndex component_id)
    {
        std::unordered_set<LocalIndex> nodes_in_component_set;
        std::vector<LocalIndex> elements_in_component;
        for (auto& e_ptr : mesh.elements())
        {
            if (e_ptr->get_component_id() == component_id)
            {
                for (auto& n_ptr : *e_ptr)
                {
                    nodes_in_component_set.insert(n_ptr->get_local_number());
                }
                elements_in_component.push_back(e_ptr->get_local_number());
            }
        }

        return {
            std::vector<LocalIndex>(nodes_in_component_set.begin(), nodes_in_component_set.end()),
            elements_in_component
        };
    }

    template<typename MeshType>
    std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> find_in_component(
        const MeshType& mesh,
        const std::vector<LocalIndex>& subset_nodes_to_search,
        const std::vector<LocalIndex>& subset_elements_to_search,
        const ComponentIndex component_id)
    {
        std::map<LocalIndex, LocalIndex> mesh_to_local_numbering_map;
        for (auto i : utility::indices(subset_nodes_to_search))
        {
            mesh_to_local_numbering_map[subset_nodes_to_search[i]] = i;
        }   

        std::unordered_set<LocalIndex> nodes_in_component_set;
        std::vector<LocalIndex> elements_in_component;
        for (auto i : utility::indices(subset_elements_to_search))
        {
            auto& e_local_id = subset_elements_to_search[i];
            auto e_ptr = mesh.elements()[e_local_id];

            if (e_ptr->get_component_id() == component_id)
            {
                for (auto& n_ptr : *e_ptr)
                {
                    nodes_in_component_set.insert(mesh_to_local_numbering_map[n_ptr->get_local_number()]);
                }
                elements_in_component.push_back(i);
            }
        }

        return {
            std::vector<LocalIndex>(nodes_in_component_set.begin(), nodes_in_component_set.end()),
            elements_in_component
        };
    }

    template<typename ElementContainer>
    std::set<ComponentIndex> get_component_ids(const ElementContainer & elements)
    {
        std::set<ComponentIndex> components_set;
        for (auto& e_ptr : elements)
        {
            components_set.insert(e_ptr->get_component_id());
        }
        return components_set;
    }

    template<typename Node_, typename ElementContainer>
    std::pair<LocalIndex, LocalIndex> get_n_nodes_and_elements_in_component(const ElementContainer & elements, const ComponentIndex component_id)
    {
        std::set<Node_ *> nodes_in_component_set;
        LocalIndex n_elements_in_component = 0;
        for (auto& e_ptr : elements) {
            if (e_ptr->get_component_id() == component_id) {
                for (auto& n_ptr : *e_ptr) {
                    nodes_in_component_set.insert(n_ptr);
                }
                n_elements_in_component++;
            }
        }
        return { static_cast<LocalIndex>(nodes_in_component_set.size()), n_elements_in_component };
    }

    template<typename NodeContainer>
    std::array<FloatType, 3> min_coord(const NodeContainer & nodes)
    {
        if (nodes.size() > 0)
        {
            Point min = nodes.front()->get_coord();
            for (auto& n_ptr : nodes) {
                for (auto i = 0; i < 3; i++) {
                    if ((*n_ptr)[i] < min[i]) {
                        min[i] = (*n_ptr)[i];
                    }
                }
            }
            return min;
        }
        else
        {
            return { 0, 0, 0 };
        }
    }
    template<typename NodeContainer>
    std::array<FloatType, 3> max_coord(const NodeContainer & nodes)
    {
        if (nodes.size() > 0)
        {
            Point max = nodes.front()->get_coord();
            for (auto& n_ptr : nodes) {
                for (auto i = 0; i < 3; i++) {
                    if ((*n_ptr)[i] > max[i]) {
                        max[i] = (*n_ptr)[i];
                    }
                }
            }
            return max;
        }
        else
        {
            return { 0, 0, 0 };
        }
    }

    template<typename MeshType>
    PartitionIndex get_num_partitions(const MeshType & mesh)
    {
        // Get the number of partitions first
        std::set<PartitionIndex> partitions;
        for (auto& n_ptr : mesh.nodes())
        {
            partitions.insert(n_ptr->get_parent_partition());
        }

        // Ensure the partitions are contiguous
        PartitionIndex curr_id = 0;
        for (auto id : partitions)
        {
            if (id != curr_id)
            {
                throw _UTIL_MSG_EXCEPT("Partition " + std::to_string(curr_id) + " is missing from the mesh.");
            }
            curr_id++;
        }

        return static_cast<PartitionIndex>(partitions.size());
    }

    template<typename MeshType>
    void sort_by_partition(const MeshType & mesh)
    {
        using NodeUniquePtr = std::unique_ptr<typename MeshType::NodeTypename>;
        using ElementUniquePtr = std::unique_ptr<typename MeshType::ElementTypename>;

        std::sort(mesh.nodes().begin(), mesh.nodes().end(),
            [](const NodeUniquePtr& n1, const NodeUniquePtr& n2) -> bool {
            return n1->get_parent_partition() < n2->get_parent_partition();
        });

        std::sort(mesh.nodes().begin(), mesh.nodes().end(),
            [](const ElementUniquePtr& e1, const ElementUniquePtr& e2) -> bool {
            return e1->get_parent_partition() < e2->get_parent_partition();
        });
    }

    template<typename MeshType>
    typename MeshType::ElementTypename* find_element_with_point(
        const MeshType& mesh,
        const simple_math::Point x)
    {
        for (auto& e_ptr : mesh.elements())
        {
            if (e_ptr->contains_coordinate(x))
            {
                return e_ptr.get();
            }
        }
        return nullptr;
    }

    template<typename MeshType, typename NodeTypename>
    std::vector<ElementFacePair> find_faces_from_node_set(
        const MeshType& mesh,
        const std::set<NodeTypename*>& node_set)
    {
        std::vector<ElementFacePair> faces_in_node_set;
        for (auto& e_ptr : mesh.elements())
        {
            auto& element_conn = e_ptr->get_conn();
            const std::vector<FaceInfo> face_info = get_face_info(e_ptr->get_type());
            for (auto face_index : utility::indices(face_info))
            {
                auto& face_node_indices = face_info[face_index].connectivity;
                if (std::all_of(
                    face_node_indices.begin(), face_node_indices.end(),
                    [&element_conn, &node_set](const ConnIndex& i) -> bool {
                        return node_set.count(element_conn[i]) != 0;
                    }))
                {
                    faces_in_node_set.push_back({ e_ptr->get_local_number(), static_cast<ConnIndex>(face_index) });
                }
            }
        }
        return faces_in_node_set;
    }
}