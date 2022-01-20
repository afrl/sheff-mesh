#pragma once

namespace mesh
{
/// @brief The types of data fields supported by MeshAssociatedData
enum class DataFieldType
{
    NODE,
    ELEMENT,
    NODE_BY_ELEMENT,
    NODE_BY_COMPONENT,
    QUAD
};

}