#pragma once

#include "mesh/Defines.h"

#include <string>
#include <vector>

namespace mesh
{
/// @brief A simple struct to describe elements in a group
struct ElementGroup
{
    std::string name;
    std::vector<LocalIndex> element_indices;

    bool operator==(const ElementGroup& other) const
    {
        return name == other.name && element_indices == other.element_indices;
    }

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("name", name),
          utility::named_value("n_elements", element_indices.size()),
          utility::named_value("element_indices", utility::slice(element_indices)));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(name);
        std::size_t n_elements;
        s(n_elements);
        element_indices.resize(n_elements);
        s(utility::Slice<LocalIndex>(element_indices.data(), element_indices.size()));
    }
};
}