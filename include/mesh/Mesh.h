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
    Mesh(const Mesh& other);
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

    /// Indicates if this instance contains the nodes/elements for the entire global mesh or just a
    /// partition.
    bool contains_all_partitions() const;

    /// Const access to the number of partitions
    const PartitionIndex& n_partitions() const;

    /// Const access to the partition that this instance represents.  Note: if this instance
    /// represents all partitions, then the returned value will be
    /// numeric_limits<PartitionIndex>::max(), but contains_all_partitions provides a more
    /// convenient test for this case.
    const PartitionIndex& partition_id() const;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Mutator methods

    /// Copies a subset of elements from another mesh and adds them to the end of this container and
    /// also copies any nodes that appear in the element connectivies being copied.
    /// @important No checks are performed to see if a node already exists, so duplicate nodes could
    /// be created.
    /// @important Calling this will invalidate n_global_nodes_ and n_global_elements_.  Be sure to
    /// call update_global_counts after modifications to the mesh to restore the correct # of
    /// nodes/elements.
    void copy(const Mesh& other, const std::vector<LocalIndex> element_indices_to_copy);

    /// Copies a subset of elements from another mesh and adds them to the end of this container and
    /// also copies any nodes that appear in the element connectivies being copied.
    /// @important No checks are performed to see if a node already exists, so duplicate nodes could
    /// be created.
    /// @important Calling this will invalidate n_global_nodes_ and n_global_elements_.  Be sure to
    /// call update_global_counts after modifications to the mesh to restore the correct # of
    /// nodes/elements.
    void copy(const Mesh& other, const std::set<LocalIndex> element_indices_to_copy);

    /// @brief See other definitions of copy for documentation (this only copies one element).
    void copy(const Mesh& other, const LocalIndex& element_index_to_copy);

    /// Copies all elements and nodes from another mesh and adds them to the end of this container
    /// and also copies any nodes that appear in the element connectivies being copied.
    /// @important No checks are performed to see if a node already exists, so duplicate nodes could
    /// be created.
    /// @important Calling this will invalidate n_global_nodes_ and n_global_elements_.  Be sure to
    /// call update_global_counts after modifications to the mesh to restore the correct # of
    /// nodes/elements.
    void copy(const Mesh& other);

    /// Sorts nodes by global number and sorts elements to arrange the owned elements first,
    /// followed by ghost elements sorted by partition.
    /// Note: this_partition_id_ determines which partitions elements gets put at the front.
    void sort();

    /// @important If distributed, calling this will only clear this rank's mesh and n_global_nodes_
    /// and n_global_elements_ will be incorrect.  Be sure to call update_global_counts after
    /// modifications to the mesh to restore the correct # of nodes/elements.
    void clear();

    /// @important If distributed, calling this will only affect this rank's mesh and
    /// n_global_nodes_ and n_global_elements_ will be incorrect.  Be sure to call
    /// update_global_counts after modifications to the mesh to restore the correct # of
    /// nodes/elements.
    Mesh& operator=(const Mesh& other);
    /// @important If distributed, calling this will only affect this rank's mesh and
    /// n_global_nodes_ and n_global_elements_ will be incorrect.  Be sure to call
    /// update_global_counts after modifications to the mesh to restore the correct # of
    /// nodes/elements.
    Mesh& operator=(Mesh&& other) noexcept;

    /// Mutable access to the number of partitions (can be used to update this if the mesh is be
    /// partitioned)
    PartitionIndex& n_partitions();

    /// Mutable access to the partition that this instance represents (can be used to update this if
    /// the mesh is be partitioned).  Note: if this instance represents all partitions, then the
    /// returned value will be numeric_limits<PartitionIndex>::max(), but contains_all_partitions
    /// provides a more convenient test for this case.
    PartitionIndex& partition_id();

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Serailization methods

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("n_partitions", n_partitions_),
          utility::named_ref("this_partition_id", this_partition_id_),
          utility::named_value("nodes", nodes_),
          utility::named_value("elements", elements_));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(n_partitions_,
          this_partition_id_,
          nodes_,
          elements_);
    }

  private:
    /// Indicates how many partitions the global mesh is divided into (this object may or may not
    /// contain all the partitions)
    PartitionIndex n_partitions_;
    /// Indicates which partition this mesh object represents.  If this mesh object contains the
    /// entire mesh, then the value will be numeric_limits<PartitionIndex>::max()
    PartitionIndex this_partition_id_;


    NodeContainer nodes_;
    ElementContainer elements_;
};
}