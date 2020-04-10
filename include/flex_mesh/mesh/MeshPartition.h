#pragma once

#include "Mesh.h"

#include "utility/MPIUtils.h"
#include "utility/GlobalVariables.h"

#include "flex_mesh/algorithms/Search.h"

#ifdef USE_FEAIO
#    include "fea_io/mesh/MeshEnums.h"
#    include "fea_io/mesh/MeshInfoStructs.h"
#endif

namespace flex_mesh
{
    /// An extension of the Mesh class for distributed meshes, where each
    /// instance of MeshPartition only stores a partition of the total mesh.
    template<typename Node_, typename Element_>
    class MeshPartition : public Mesh<Node_, Element_>
    {
    public:
        using NodeTypename = Node_;
        using ElementTypename = Element_;
        using NodeGroup = std::set<Node_ *>;
        using NodeGroupMap = std::map<std::string, NodeGroup>;
        using NodeGroupMapIt = typename NodeGroupMap::iterator;
        using ElementGroup = std::set<Element_>;
        using ElementGroupMap = std::map<std::string, ElementGroup>;
        using ElementGroupMapIt = typename std::map<std::string, ElementGroup>::iterator;

        MeshPartition();

        PartitionIndex get_partition_id() const { return this_partition_; }
        /// Returns the total number of elements across all ranks, since last call of sync()
        GlobalIndex get_n_total_elements() const { return n_total_elements_; }
        /// Returns the total number of nodes across all ranks, since last call of sync()
        GlobalIndex get_n_total_nodes() const { return n_total_nodes_; }
        /// Retursn the number of nodes owned by this rank
        GlobalIndex get_n_owned_elements() const;
        GlobalIndex get_n_owned_nodes() const;

        void set_partition_id(const PartitionIndex id) { this_partition_ = id; }
        void set_n_total_elements(const GlobalIndex n_elements) { n_total_elements_ = n_elements; }
        void set_n_total_nodes(const GlobalIndex n_nodes) { n_total_nodes_ = n_nodes; }

        #ifdef USE_FEAIO
        /// Returns the information need about the mesh to either write the mesh in a distribtued
        /// manner or write data for the mesh.
        /// Note: Right now, the MeshInfo::components is left blank, which is fine for
        ///       for writing data, but is not appropriate for writing the mesh.
        fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_mesh_info();

        fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_partition_info();
        #endif

        /// Sync the ghost elements/nodes across all of the ranks.  Algorithm is:
        /// 1. Check if any rank added or removed nodes/elements (if not return)
        /// 2. Update the total number of nodes/elements
        ///
        /// Note: currently, the algorithm ends here, but in the future it should probably:
        /// 3. Communicate the global node/element numbers that this rank requires to all other ranks
        /// 4. Loop over each other rank and serialize any new nodes/elements that rank created but
        ///    that will be needed by that other rank and send them.
        /// 5. Loop over each other rank and receive any node/elements that are owned by the other
        ///    rank but needed as a ghost for this rank.
        ///
        /// Returns a bool indicating if the total number of nodes or elements have
        /// changed since the last call to sync()
        /// \todo finish implementation
        bool sync();

    protected:
        PartitionIndex this_partition_;
        /// The number of elements in the total mesh since the last call to sync()
        GlobalIndex n_total_elements_;
        /// The number of nodes in the total mesh since the last call to sync()
        GlobalIndex n_total_nodes_;
    };

    template<typename Node_, typename Element_>
    inline MeshPartition<Node_, Element_>::MeshPartition()
    {
        this_partition_ = utility::MPIUtils::get_local_rank();
    }

    template<typename Node_, typename Element_>
    inline GlobalIndex MeshPartition<Node_, Element_>::get_n_owned_elements() const
    {
        GlobalIndex n_owned_elements = 0;
        for (auto& e : this->elements())
        {
            if (e->get_parent_partition() == this_partition_)
            {
                n_owned_elements++;
            }
        }
        return n_owned_elements;
    }

    template<typename Node_, typename Element_>
    inline GlobalIndex MeshPartition<Node_, Element_>::get_n_owned_nodes() const
    {
        GlobalIndex n_owned_nodes = 0;
        for (auto& n : this->nodes())
        {
            if (n->get_parent_partition() == this_partition_)
            {
                n_owned_nodes++;
            }
        }
        return n_owned_nodes;
    }

    #ifdef USE_FEAIO
    template<typename Node_, typename Element_>
    inline fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> MeshPartition<Node_, Element_>::get_fea_io_mesh_info()
    {
        sync();

        fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> mesh_info;
        mesh_info.total_num_nodes = n_total_nodes_;
        mesh_info.total_num_elements = n_total_elements_;
        mesh_info.dimensionality = utility::GlobalVariables::analysis_dimensionality__;
        mesh_info.extent = fea_io::Extent<FloatType>();
        mesh_info.connectivity_scheme = fea_io::ConnectivityScheme::FlexMesh;
        mesh_info.components = {};
        mesh_info.total_num_partitions = utility::MPIUtils::get_num_local_ranks();
        return mesh_info;
    }

    template<typename Node_, typename Element_>
    inline fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> MeshPartition<Node_, Element_>::get_fea_io_partition_info()
    {
        fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> partition_info;

        partition_info.partition_id = this_partition_;
        partition_info.num_nodes = static_cast<LocalIndex>(this->nodes().size());
        partition_info.num_elements = static_cast<LocalIndex>(this->elements().size());

        partition_info.total_num_connectivity_entries = 0;
        for (auto& e_ptr : this->elements())
        {
            partition_info.total_num_connectivity_entries += e_ptr->get_num_nodes();
        }

        auto min = min_coord(this->nodes());
        auto max = max_coord(this->nodes());
        partition_info.extent = fea_io::Extent<FloatType>(min, max);

        auto component_ids = flex_mesh::get_component_ids(this->elements());
        partition_info.num_subcomponents = static_cast<LocalIndex>(component_ids.size());
        for (auto component_id : component_ids)
        {
            auto n_nodes_elements = get_n_nodes_and_elements_in_component<Node_>(this->elements(), component_id);
            partition_info.subcomponents.push_back(fea_io::SubcomponentInfo<LocalIndex>(
                component_id, n_nodes_elements.first, n_nodes_elements.second));
        }

        return partition_info;
    }
    #endif

    template<typename Node_, typename Element_>
    inline bool MeshPartition<Node_, Element_>::sync()
    {
        if (utility::MPIUtils::is_distributed())
        {
            bool n_total_changed = false;

            auto curr_total_nodes = utility::MPIUtils::sum_across_ranks(get_n_owned_nodes());
            n_total_changed |= curr_total_nodes != n_total_nodes_;
            n_total_nodes_ = curr_total_nodes;

            auto curr_total_elements = utility::MPIUtils::sum_across_ranks(get_n_owned_elements());
            n_total_changed |= curr_total_elements != n_total_elements_;
            n_total_elements_ = curr_total_elements;

            //throw _UTIL_MSG_EXCEPT_C("TODO finish implementation");

            return n_total_changed;
        }
        else
        {
            bool n_total_changed = false;
            n_total_changed |= n_total_nodes_ != get_n_owned_nodes();
            n_total_nodes_ = get_n_owned_nodes();
            n_total_changed |= n_total_elements_ != get_n_owned_elements();
            n_total_elements_ = get_n_owned_elements();
            return n_total_changed;
        }

        return false;
    }
}
