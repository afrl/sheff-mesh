#pragma once

#include "simple_math/Defines.h"
#include "simple_math/EigenDefines.h"

#include <array>
#include <cinttypes>
#include <limits>

namespace mesh
{
using ConnIndex      = uint8_t;
using LocalIndex     = uint32_t;
using GlobalIndex    = uint32_t;
using PartitionIndex = uint16_t;
using ComponentIndex = uint16_t;
using FloatType      = simple_math::FloatType;
using Point          = simple_math::Point;
using Vector3D       = simple_math::Vector3D;

/// A value used to indicate an index is not set
constexpr LocalIndex NOT_SET_LocalIndex = std::numeric_limits<LocalIndex>::max();
constexpr GlobalIndex NOT_SET_GlobalIndex = std::numeric_limits<GlobalIndex>::max();
constexpr ComponentIndex NOT_SET_ComponentIndex = std::numeric_limits<ComponentIndex>::max();

enum class TopologyType
{
    NODE,
    EDGE,
    FACE,
    ELEMENT,
    NONE ///< Exists to allow an empty/default ElementTopologyComponent
};

/// @brief Describes a component of an element's topology (nodes, edges, and faces)
struct ElementTopologyComponent
{
    /// Will either be Node, Edge, or Face since this is a topology part of an element
    const TopologyType type;
    /// The index of the respective node, face or edge relative to the element
    const ConnIndex index;

    ElementTopologyComponent()
      : type(TopologyType::NONE)
      , index(0)
    {
    }
    ElementTopologyComponent(const TopologyType t, const ConnIndex i)
      : type(t)
      , index(i)
    {
    }
};
}