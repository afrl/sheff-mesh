#pragma once

#include "mesh/Defines.h"

#include <string>
#include <vector>

namespace mesh
{
/// @brief A simple struct to describe nodes in a group
struct NodeGroup
{
    std::string name;
    std::vector<LocalIndex> node_indices;

    bool operator==(const NodeGroup& other) const
    {
        return name == other.name && node_indices == other.node_indices;
    }

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("name", name),
          utility::named_value("n_nodes", node_indices.size()),
          utility::named_value("node_indices", utility::slice_c(node_indices)));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(name);
        std::size_t n_nodes;
        s(n_nodes);
        node_indices.resize(n_nodes);
        s(utility::Slice<LocalIndex>(node_indices.data(), node_indices.size()));
    }
};
}
