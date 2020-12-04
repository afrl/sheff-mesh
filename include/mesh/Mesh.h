#pragma once

#include "mesh/element/Element.h"
#include "mesh/node/NodeAccessor.h"
#include "mesh/node/NodeContainer.h"
#include "mesh/node/NodeIterator.h"

#include <set>

namespace mesh
{
/// @brief A simple class providing an interface to a mesh composed a NodeContainer and
/// ElementContainer.
class Mesh
{
  public:
    Mesh();
    /// @brief Move constructor (swaps with other)
    Mesh(Mesh&& other) noexcept;

    bool operator==(const Mesh& other) const;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Accessor methods

    inline NodeContainer& nodes()
    {
        return nodes_;
    }
    inline const NodeContainer& nodes() const
    {
        return nodes_;
    }
    inline ElementContainer& elements()
    {
        return elements_;
    }
    inline const ElementContainer& elements() const
    {
        return elements_;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Mutator methods

    /// Copies a subset of elements from another mesh and adds them to the end of this container and
    /// also copies any nodes that appear in the element connectivies being copied.
    /// @important No checks are performed to see if a node already exists, so duplicate nodes could
    /// be created.
    void copy(const Mesh& other, const std::set<LocalIndex> element_indices_to_copy);

    /// Copies all elements and nodes from another mesh and adds them to the end of this container
    /// and also copies any nodes that appear in the element connectivies being copied.
    /// @important No checks are performed to see if a node already exists, so duplicate nodes could
    /// be created.
    void copy(const Mesh& other);

    /// Sorts nodes by global number and sorts elements to arrange the owned elements first,
    /// followed by ghost elements sorted by partition.
    /// @param this_partition - determines which partitions elements gets put at the front.
    ///        If you want to simply sort elements by partition, this should be 0.
    void sort(const PartitionIndex this_partition);

    void clear();

    Mesh& operator=(const Mesh& other);
    Mesh& operator=(Mesh&& other) noexcept;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Serailization methods

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("nodes", nodes_), utility::named_value("elements", elements_));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(nodes_, elements_);
    }

  private:
    NodeContainer nodes_;
    ElementContainer elements_;
};
}