#pragma once

#include "utility/Exceptions.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/Range.h"
#include "utility/MPIUtils.h"

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/mesh/MeshPartition.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/maps/AdjacencyLists.h"

#include "flex_mesh/algorithms/Search.h"

#ifdef SWIG
#define USE_FEAIO
#endif

#ifdef USE_FEAIO
#include "fea_io/meta/MetaWriter.h"
#include "fea_io/mesh/MeshWriter.h"
#include "fea_io/data/DataWriter.h"
#include "fea_io/mesh/MeshReader.h"
#include "fea_io/meta/MetaReader.h"

#ifdef HAVE_FILESYSTEM
#include<experimental/filesystem>
#else
#include "fea_io/common/Filesystem.h"
#endif

#endif

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef USE_FEAIO

namespace flex_mesh
{
    namespace fs = std::experimental::filesystem;

    /// This collection of mesh algorithms pertains to reading/writing with the fea_io format
    template<typename MeshType>
    class FEAIOInterface
    {
    public:
        using NodePtr = typename MeshType::NodeTypename*;
        using ElementPtr = typename MeshType::ElementTypename*;
        using NodeTypename = typename MeshType::NodeTypename;
        using ElementTypename = typename MeshType::ElementTypename;
        using NodeDataCallback = std::function<std::vector<double>(std::string, NodePtr)>;
        using ElementDataCallback = std::function<std::vector<double>(std::string, ElementPtr)>;
        using NodeByElemDataCallback = std::function<std::vector<std::vector<double>>(std::string, ElementPtr)>;

        using Partition = fea_io::Partition<GlobalIndex, LocalIndex, FloatType>;
        using MeshReader = fea_io::MeshReader<GlobalIndex, LocalIndex, FloatType>;

        /// Writes the mesh, material information as a element group file, and
        /// the orientation information as the new compact, binary Beta IO Format.
        /// The files are naturally partitioned and ready for HPC.
        /// Note: the mesh should be partitioned before this method is called,
        ///  otherwise, only 1 partition will be created.
        /// Pass the base file name, and the correct extension will
        /// be appended (mesh: *.pfecm, group information for materials: *.pfecg,
        /// orientation: *.pfeco).  It will create a meta file (*.pfecx) with a single
        /// time step and references null data to allow processing, reading
        /// into beta, or visualization.
        static void write(MeshType& mesh,
            const std::string & base_directory,
            const std::string & meta_file_name,
            const std::string & mesh_name,
            const std::map<ComponentIndex, std::string> & component_names = {});

        static void write(MeshType & mesh,
            const std::string mesh_filepath,
            const std::map<ComponentIndex, std::string> & component_names = {});

        static void write_multipart(MeshType & mesh,
            const std::string & base_directory,
            const std::string & meta_file_name,
            const std::map<ComponentIndex, std::string> & component_names,
            const LocalIndex n_ghost_levels = 1);

        /*
        static void write(
            MeshType& mesh,
            const std::string& base_file_path,
            const std::string& meta_name,
            const std::map<ComponentIndex, std::string>& component_names,
            std::vector<fea_io::ArrayInfo> array_info,
            NodeDataCallback node_data_callback,
            ElementDataCallback element_data_callback,
            NodeByElemDataCallback node_by_elem_data_callback);
            */

        static void read_partition(MeshType& mesh,
            const std::string & meta_filepath,
            const PartitionIndex partition_id = 0,
            const int time_step = 0);

        /// Computes the information for the mesh needed for writing with fea_io
        static fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_mesh_info(
            const MeshType& mesh, const PartitionIndex n_partitions,
            const std::map<ComponentIndex, std::string>& component_names);
        /// Computes the partition information needed for writing in parallel with fea_io
        /// Note: currently, this interface just writes in serial.
        static std::vector<fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType>> get_fea_io_partition_infos(
            MeshType& mesh, const std::vector<PartitionIndex>& partition_ids);

    private:
        static flex_mesh::LCSAssignment convert_lcs_assignment(fea_io::LCSAssignment lcs_assignment);
        static fea_io::LCSAssignment convert_lcs_assignment(flex_mesh::LCSAssignment lcs_assignment);
        static flex_mesh::ElementType convert_element_type(fea_io::ElementType type);
        static fea_io::ElementType convert_element_type(flex_mesh::ElementType type);

        /// Gets the fea_io::Partition object that contains all of the data needed to write the
        /// mesh for this parititon.
        /// \attention nodes_in_partition is non-const because their local numbering will be updated.
        /// Note: The compiler should be able to deduce NodeContainer and ElementContainer.
        template<typename NodeContainer, typename ElementContainer>
        static Partition get_fea_io_partition(
            NodeContainer & nodes_in_partition,
            const ElementContainer & elements_in_partition,
            const PartitionIndex partition_id);

        /// Given a set of nodes and elements, this method will add ghost elements to the set
        /// and add any new nodes that appear in the ghosts to set of nodes.
        static void add_ghost_layers(const LocalIndex n_ghost_levels,
            const ElementAdjacencyList<ElementTypename> adjacency_list,
            std::set<NodeTypename> & nodes, std::set<ElementTypename> & elements);

        static void write_kernel(MeshType& partition_mesh,
            const PartitionIndex & partition_id,
            fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer);
    };

    /////////////////
    // Definitions //
    /////////////////

    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::write(MeshType & mesh, const std::string & base_directory,
        const std::string & meta_file_name, const std::string & mesh_name,
        const std::map<ComponentIndex, std::string> & component_names)
    {
        // Start a new meta file
        fea_io::MetaWriter meta_writer;
        std::string meta_filepath = (fs::path(base_directory) / meta_file_name).string();
        meta_writer.open(meta_filepath, meta_file_name);
        std::string mesh_filepath = meta_writer.add_mesh_path(0, 0, mesh_name);
        write(mesh, mesh_filepath);
    }

    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::write(MeshType & mesh, const std::string mesh_filepath,
        const std::map<ComponentIndex, std::string> & component_names)
    {
        // If the program is not distributed, then we can assume this process will write all of the partitions
        // and is given the entire mesh.
        if (!utility::MPIUtils::is_distributed())
        {
            // This method treats the mesh as a single partition.
            const PartitionIndex n_partitions = 1;
            const PartitionIndex partition_id = 0;

            fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer;
            auto mesh_info = get_fea_io_mesh_info(mesh, n_partitions, component_names);
            mesh_writer.initialize(mesh_filepath, mesh_info);

            // For a serial write, we do not need to queue up all of the partitions for writing, since
            // the file does not have to be preallocated for serial writing.
            auto partition = get_fea_io_partition(mesh.nodes(), mesh.elements(), partition_id);

            // Overwrite the parent partitions
            for (auto& n_part : partition.node_parent_partitions)
            {
                n_part = 0;
            }
            for (auto& e_part : partition.element_parent_partitions)
            {
                e_part = 0;
            }

            mesh_writer.write_partition(partition);
        }
        else
        {
            throw _UTIL_MSG_EXCEPT("FEAIOInterface::write cannot be called while running in distributed, consider calling FEAIOInterface::write_multipart instead.");
        }
    }

    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::write_multipart(MeshType & mesh, const std::string & base_directory,
        const std::string & meta_file_name, const std::map<ComponentIndex, std::string> & component_names,
        const LocalIndex n_ghost_levels)
    {
        // If the program is not distributed, then we can assume this process will write all of the partitions
        // and is given the entire mesh.
        if (!utility::MPIUtils::is_distributed())
        {
            // Get the number of partitions first
            PartitionIndex num_partitions = get_num_partitions(mesh);

            std::vector<PartitionIndex> partition_ids(num_partitions);
            // Fills vector with 0..num_partitions-1
            std::iota(partition_ids.begin(), partition_ids.end(), 0);

            // Start a new meta file
            fea_io::MetaWriter meta_writer;
            meta_writer.open(base_directory, meta_file_name);
            std::string mesh_filepath = meta_writer.add_mesh_path(0, 0, meta_file_name);

            fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer;
            auto mesh_info = get_fea_io_mesh_info(mesh, component_names);
            mesh_writer.initialize(mesh_filepath, mesh_info);

            // Create an adjacency list to quickly get the neighbors of elements while finding the
            // layer of ghost elements.
            ElementAdjacencyList<MeshType> adjacency_list;
            adjacency_list.initialize(mesh);

            for (auto& partition_id : partition_ids)
            {
                // Get the nodes and elements that should appear in the partition corresponding to
                // partition_id, including ghost elements.
                auto nodes_and_elements = get_partition_nodes_and_elements(mesh, partition_id);
                add_ghost_layers(1, adjacency_list, nodes_and_elements.first, nodes_and_elements.second);

                // For a serial write, we do not need to queue up all of the partitions for writing, since
                // the file does not have to be preallocated for serial writing.
                auto partition = get_fea_io_partition(nodes_and_elements.first, nodes_and_elements.second, partition_id);
                mesh_writer.write_partition(partition);
            }
        }
        
    }

    /*
    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::write(MeshType & mesh, const std::string & base_file_path, const std::string& meta_name, const std::map<ComponentIndex, std::string>& component_names, std::vector<fea_io::ArrayInfo> array_info, NodeDataCallback node_data_callback, ElementDataCallback element_data_callback, NodeByElemDataCallback node_by_elem_data_callback)
    {
        std::sort(nodes_in_component.begin(), nodes_in_component.end(),
            [](const NodePtr& n1, const NodePtr& n2) -> bool {
            return n1->get_local_number() < n2->get_local_number();
        });
        std::sort(elements_in_component.begin(), elements_in_component.end(),
            [](const ElementPtr& e1, const ElementPtr& e2) -> bool {
            return e1->get_global_number() < e2->get_global_number();
        });

        SubcomponentInfo<int16_t> sc0 = { 0, 4, 1 };
        SubcomponentInfo<int16_t> sc1 = { 0, 4, 1 };
        SubcomponentInfo<int16_t> sc2 = { 1, 4, 1 };
        SubcomponentInfo<int16_t> sc3 = { 1, 4, 1 };

        PartitionInfo<int16_t, int16_t, float> part0 = {
            0,
            6,
            2,
            8,
            Extent<float>(),
            2,
            { sc0, sc2 }
        };

        PartitionInfo<int16_t, int16_t, float> part1 = {
            1,
            6,
            2,
            8,
            Extent<float>(),
            2,
            { sc1, sc3 }
        };

        partitions = { part0, part1 };

        get_test_info(mesh_info, all_partitions);

        MeshWriter<int16_t, int16_t, float> writer;
        writer.initialize("test_serial", mesh_info);

        writer.write_partition(get_part0());
        writer.write_partition(get_part1());

        // Set up the group file if it is needed
        std::unique_ptr<PFECNativeInterface::GroupWriter> groupWriterPtr;
        std::vector<std::string> nodeGroupNames = mesh.get_node_group_names();
        std::vector<std::string> elementGroupNames = mesh.get_element_group_names();
        if (nodeGroupNames.size() != 0 || elementGroupNames.size() != 0) {
            groupWriterPtr = metaWriter.get_group_writer(BaseFileName + ".pfecg");
            groupWriterPtr->write_header_info(partNums, nodeGroupNames, elementGroupNames);
        }

        // Set up the data file if it is needed
        std::unique_ptr<PFECNativeInterface::DataWriter> dataWriterPtr;
        if (outputMeta.size() > 0) {
            dataWriterPtr = metaWriter.get_data_writer(BaseFileName + ".pfecd");
            dataWriterPtr->initialize_streams(partNums, outputMeta);
        }

        // Do a special case if there is only one partition, to avoid creating a lot of over head in the partitioning
        if (numParts == 1) {
            // Update all nodes to have the same local numbering, since the elements will use these numbers
            mesh.update_global_numbering();

            // Write mesh file
            typedef typename MeshType::UniquePtrNodeListItConst MeshNodeIter;
            typedef typename MeshType::UniquePtrElementListItConst MeshElemIter;
            meshWriterPtr->write_part_nodal2(0, PFECNodeIterator<MeshNodeIter>(mesh.nodes_begin()), PFECNodeIterator<MeshNodeIter>(mesh.nodes_end()));
            meshWriterPtr->write_part_element2(0, PFECElementIterator<MeshElemIter>(mesh.elements_begin()), PFECElementIterator<MeshElemIter>(mesh.elements_end()));

            // Write orientation data
            if (needOrientations) {
                std::vector<std::vector<std::array<double, 4>>> localQuaterions;
                localQuaterions.resize(mesh.elements_size());
                for (auto elemIt = mesh.elements_begin(); elemIt != mesh.elements_end(); elemIt++) {
                    auto localNum = (*elemIt)->get_local_number();
                    auto& orientations = (*elemIt)->get_orientations();
                    localQuaterions[localNum].resize(orientations.size());
                    for (auto j = 0; j < orientations.size(); j++)
                        localQuaterions[localNum][j] = orientations[j].get_quaterion();
                }
                orientWriterPtr->write_part_orientation_data(0, 0, localQuaterions);
            }

            if (groupWriterPtr != 0) {
                std::vector<PFECNativeInterface::PFEC_IND> indices;
                indices.reserve(mesh.nodes_size());
                for (std::string groupName : nodeGroupNames) {
                    for (auto node : mesh.get_node_group(groupName))
                        indices.push_back(node->get_global_number());
                    groupWriterPtr->write_part_node_group(0, groupName, indices, true);
                    indices.clear();
                }
                indices.reserve(mesh.elements_size());
                for (std::string groupName : elementGroupNames) {
                    for (auto element : mesh.get_element_group(groupName))
                        indices.push_back(element->get_global_number());
                    groupWriterPtr->write_part_elem_group(0, groupName, indices, true);
                    indices.clear();
                }
            }

            // Write data using callback if there is any
            for (auto meta : outputMeta) {
                if (meta.type == PFECNativeInterface::DataType::Nodal) {
                    UniqueDataVec data;
                    for (auto& n : mesh.nodes())
                        data.push_back(nodeDataCallback(meta.id, n.get()));
                    dataWriterPtr->write_part_unique_data(0, meta.id, data);
                }
                else if (meta.type == PFECNativeInterface::DataType::Element) {
                    UniqueDataVec data;
                    for (auto& e : mesh.elements())
                        data.push_back(elementDataCallback(meta.id, e.get()));
                    dataWriterPtr->write_part_unique_data(0, meta.id, data);
                }
                else if (meta.type == PFECNativeInterface::DataType::NodalByElement) {
                    NodeByElemDataVec data;
                    for (auto& e : mesh.elements())
                        data.push_back(nodalByElemDataCallback(meta.id, e.get()));
                    dataWriterPtr->write_part_by_elem_data(0, meta.id, data);
                }
                else {
                    auto dataNameMap = PFECNativeInterface::get_data_type_name_map();
                    throw _BETA_UTIL_MSG_EXCEPT2("Writing this type of data is not supported yet in the write_pfec_files method: " + dataNameMap[meta.type]);
                }
            }
        }
        else { // There is more than one partition
                // Loop over partitions and write the mesh/group information
            for (int i = 0; i < numParts; i++) {
                // Update the local number and get the nodes/elements in the partition
                auto partObjects = mesh.update_local_numbering(i);
                NodeSet& localNodes = partObjects.first;
                ElementSet& localElements = partObjects.second;

                // Write the mesh
                typedef typename NodeSet::iterator NodeSetIter;
                typedef typename ElementSet::iterator ElemSetIter;
                meshWriterPtr->write_part_nodal2(i, PFECNodeIterator<NodeSetIter>(localNodes.begin()), PFECNodeIterator<NodeSetIter>(localNodes.end()));
                meshWriterPtr->write_part_element2(i, PFECElementIterator<ElemSetIter>(localElements.begin()), PFECElementIterator<ElemSetIter>(localElements.end()));

                // Write orientation data
                if (needOrientations) {
                    std::vector<std::vector<std::array<double, 4>>> localQuaterions;
                    localQuaterions.resize(localElements.size());
                    for (auto elemIt = localElements.begin(); elemIt != localElements.end(); elemIt++) {
                        auto localNum = (*elemIt)->get_local_number();
                        auto& orientations = (*elemIt)->get_orientations();
                        localQuaterions[localNum].resize(orientations.size());
                        for (auto j = 0; j < orientations.size(); j++)
                            localQuaterions[localNum][j] = orientations[j].get_quaterion();
                    }
                    orientWriterPtr->write_part_orientation_data(i, 0, localQuaterions);
                }

                // Write data using callback if there is any
                for (auto meta : outputMeta) {
                    if (meta.type == PFECNativeInterface::DataType::Nodal) {
                        UniqueDataVec data;
                        for (auto& n : mesh.nodes())
                            data.push_back(nodeDataCallback(meta.id, n.get()));
                        dataWriterPtr->write_part_unique_data(i, meta.id, data);
                    }
                    else if (meta.type == PFECNativeInterface::DataType::Element) {
                        UniqueDataVec data;
                        for (auto& e : localElements)
                            data.push_back(elementDataCallback(meta.id, e));
                        dataWriterPtr->write_part_unique_data(i, meta.id, data);
                    }
                    else if (meta.type == PFECNativeInterface::DataType::NodalByElement) {
                        NodeByElemDataVec data;
                        for (auto& e : localElements)
                            data.push_back(nodalByElemDataCallback(meta.id, e));
                        dataWriterPtr->write_part_by_elem_data(i, meta.id, data);
                    }
                    else {
                        auto dataNameMap = PFECNativeInterface::get_data_type_name_map();
                        throw _BETA_UTIL_MSG_EXCEPT2("Writing this type of data is not supported yet in the write_pfec_files method: " + dataNameMap[meta.type]);
                    }
                }
            } // End loop over partitions
        }

        // Close writers
        meshWriterPtr->close();
        if (needOrientations)
            orientWriterPtr->finalize_streams();
        if (groupWriterPtr != 0)
            groupWriterPtr->close();
        if (dataWriterPtr != 0)
            dataWriterPtr->finalize_streams();
            
        return;
    }
    */

    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::read_partition(MeshType & mesh, const std::string & meta_filepath, const PartitionIndex partition_id, const int time_step)
    {
        fea_io::MetaReader meta_reader;
        meta_reader.open(meta_filepath);
        auto mesh_filepath = meta_reader.get_mesh_path(time_step);

        MeshReader mesh_reader(mesh_filepath);
        auto mesh_info = mesh_reader.get_mesh_info();
        if (mesh_info.connectivity_scheme != fea_io::ConnectivityScheme::FlexMesh)
        {
            throw _UTIL_MSG_EXCEPT("The mesh was written with a different connectivity scheme.  A conversion has not been written.");
        }

        auto partition = mesh_reader.read_partition(partition_id);

        for (auto i : utility::indices(partition.node_global_ids))
        {
            const auto offset = mesh_info.dimensionality * i;
            Point p = { 0, 0, 0 };
            switch (mesh_info.dimensionality)
            {
            case 3:
                p[2] = partition.node_coords[offset + 2];
            case 2:
                p[1] = partition.node_coords[offset + 1];
            case 1:
                p[0] = partition.node_coords[offset];
            }

            auto n_ptr = mesh.add_node_by_coord(p);
            n_ptr->set_global_number(partition.node_global_ids[i]);
            n_ptr->set_local_number(static_cast<LocalIndex>(i));
            n_ptr->set_parent_partition(partition.node_parent_partitions[i]);
        }
        
        for (auto i : utility::indices(partition.element_global_ids))
        {
            std::vector<typename MeshType::NodeTypename *> connectivity;
            for (auto local_node_id : partition.element_connectivities[i])
            {
                connectivity.push_back(mesh.nodes()[local_node_id]);
            }

            auto e_ptr = mesh.elements().add();
            e_ptr->initialize(convert_element_type(partition.element_types[i]), connectivity);

            e_ptr->set_global_number(partition.element_global_ids[i]);
            e_ptr->set_local_number(static_cast<LocalIndex>(i));
            e_ptr->set_parent_partition(partition.element_parent_partitions[i]);
        }

        for (auto& subcomponent : partition.subcomponents)
        {
            auto find_iter = std::find_if(std::begin(mesh_info.components), std::end(mesh_info.components),
                [id = subcomponent.component_id](const fea_io::ComponentInfo<LocalIndex> & c)
                { return c.id == id; });

            if (find_iter == std::end(mesh_info.components))
            {
                throw _UTIL_MSG_EXCEPT("A subcomponent referenced a component id that did not exist in the mesh information.");
            }
            if (mesh_info.components[0].lcs_parameterization != fea_io::LCSParameterization::Quaternion)
            {
                throw _UTIL_MSG_EXCEPT("Local coordinate systems were parameterized with something other than Quaternions. A conversion needs to be added to read_partition.");
            }

            LCSAssignment component_lcs_assignment = convert_lcs_assignment((*find_iter).lcs_assignment);
            for (auto local_elem_id : subcomponent.elements_in_subcomponent)
            {
                mesh.elements()[local_elem_id]->set_component_id(subcomponent.component_id);
                mesh.elements()[local_elem_id]->allocate_lcs(component_lcs_assignment);
                // TODO assign lcs
            }
        }
    }

    template<typename MeshType>
    inline fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> FEAIOInterface<MeshType>::get_fea_io_mesh_info(
        const MeshType & mesh, const PartitionIndex n_partitions,
        const std::map<ComponentIndex, std::string> & component_names)
    {
        // Get the component IDs
        auto component_id_set = mesh.get_component_ids();
        std::vector<ComponentIndex> component_id_vector(component_id_set.begin(), component_id_set.end());

        // Create the meta info the mesh file
        fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> mesh_info;
        mesh_info.total_num_nodes = static_cast<GlobalIndex>(mesh.nodes().size());
        mesh_info.total_num_elements = static_cast<GlobalIndex>(mesh.elements().size());
        mesh_info.dimensionality = 3;
        mesh_info.extent = fea_io::Extent<FloatType>(min_coord(mesh.nodes()), max_coord(mesh.nodes()));
        mesh_info.connectivity_scheme = fea_io::ConnectivityScheme::FlexMesh;
        for (auto& component_id : component_id_vector)
        {
            auto elements_in_component = mesh.get_component(component_id);

            // Ensure the orientation type for all the elements is the same
            auto front_lcs_assignment = elements_in_component.front()->get_lcs_assignment();
            for (auto& e_ptr : elements_in_component)
            {
                if (e_ptr->get_lcs_assignment() != front_lcs_assignment)
                {
                    throw _UTIL_MSG_EXCEPT("The elements within the same component must have the same LCSAssignment to be written to the fea_io format.");
                }
            }

            // Figure out the name of the component
            auto comp_name_find_iter = component_names.find(component_id);
            std::string component_name;
            if (comp_name_find_iter == component_names.end())
            {
                component_name = "Component" + std::to_string(component_id);
            }
            else
            {
                component_name = comp_name_find_iter->second;
            }

            mesh_info.components.push_back({ component_id, component_name,
                convert_lcs_assignment(front_lcs_assignment), fea_io::LCSParameterization::Quaternion });
        }
        mesh_info.total_num_partitions = n_partitions;
        return mesh_info;
    }

    template<typename MeshType>
    inline std::vector<fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType>> FEAIOInterface<MeshType>::get_fea_io_partition_infos(MeshType & mesh, const std::vector<PartitionIndex>& partition_ids)
    {
        // Create the partition info for the mesh file
        std::vector<fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType>> partition_infos;
        for (auto& partition_id : partition_ids) {
            // Update local indices
            auto numbering_result = mesh.update_local_numbering(partition_id);
            std::set<NodePtr>& nodes_in_partition = numbering_result.first;
            std::set<ElementPtr>& elements_in_partition = numbering_result.second;

            // Figure out what each subcomponent contains in this partition
            std::vector<fea_io::SubcomponentInfo<LocalIndex>> subcomponent_infos;
            for (auto& component_id : mesh.get_component_ids()) {
                auto nodes_in_component = find_nodes_in_component<NodePtr>(elements_in_partition, component_id);
                auto elements_in_component = find_elements_in_component(elements_in_partition, component_id);

                subcomponent_infos.push_back({
                    component_id,
                    static_cast<LocalIndex>(nodes_in_component.size()),
                    static_cast<LocalIndex>(elements_in_component.size()),
                    0 });
            }

            LocalIndex total_num_connectivity_entries = 0;
            for (auto& e_ptr : elements_in_partition)
                total_num_connectivity_entries += e_ptr->get_num_nodes();

            fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType> partition_info;
            partition_info.partition_id = partition_id;
            partition_info.num_nodes = static_cast<LocalIndex>(nodes_in_partition.size());
            partition_info.num_elements = static_cast<LocalIndex>(elements_in_partition.size());
            partition_info.total_num_connectivity_entries = total_num_connectivity_entries;
            partition_info.extent = fea_io::Extent<FloatType>(min_coord(nodes_in_partition), max_coord(nodes_in_partition));
            partition_info.num_subcomponents = static_cast<LocalIndex>(subcomponent_infos.size());
            partition_info.subcomponents = subcomponent_infos;
        }
        return partition_infos;
    }

    template<typename MeshType>
    inline flex_mesh::LCSAssignment FEAIOInterface<MeshType>::convert_lcs_assignment(fea_io::LCSAssignment lcs_assignment)
    {
        if (lcs_assignment == fea_io::LCSAssignment::Component)
            return LCSAssignment::Component;
        else if (lcs_assignment == fea_io::LCSAssignment::Element)
            return LCSAssignment::Element;
        else if (lcs_assignment == fea_io::LCSAssignment::Node)
            return LCSAssignment::Node;
        else if (lcs_assignment == fea_io::LCSAssignment::None)
            return LCSAssignment::None;
        else
            throw _UTIL_MSG_EXCEPT("Cannot convert from fea_io::LCSAssignment.");
    }

    template<typename MeshType>
    inline fea_io::LCSAssignment FEAIOInterface<MeshType>::convert_lcs_assignment(flex_mesh::LCSAssignment lcs_assignment)
    {
        if (lcs_assignment == flex_mesh::LCSAssignment::Component)
            return fea_io::LCSAssignment::Component;
        else if (lcs_assignment == flex_mesh::LCSAssignment::Element)
            return fea_io::LCSAssignment::Element;
        else if (lcs_assignment == flex_mesh::LCSAssignment::Node)
            return fea_io::LCSAssignment::Node;
        else if (lcs_assignment == flex_mesh::LCSAssignment::None)
            return fea_io::LCSAssignment::None;
        else
            throw _UTIL_MSG_EXCEPT("Cannot convert from fea_io::LCSAssignment.");
    }

    template<typename MeshType>
    inline flex_mesh::ElementType FEAIOInterface<MeshType>::convert_element_type(fea_io::ElementType type)
    {
        // 1D
        if (type == fea_io::ElementType::L1D2N)
            return flex_mesh::ElementType::L1D2N;
        else if (type == fea_io::ElementType::L1D3N)
            return flex_mesh::ElementType::L1D3N;
        // 2D
        else if (type == fea_io::ElementType::T2D3N)
            return flex_mesh::ElementType::T2D3N;
        else if (type == fea_io::ElementType::Q2D4N)
            return flex_mesh::ElementType::Q2D4N;
        else if (type == fea_io::ElementType::C2D4N)
            return flex_mesh::ElementType::C2D4N;
        else if (type == fea_io::ElementType::T2D6N)
            return flex_mesh::ElementType::T2D6N;
        else if (type == fea_io::ElementType::Q2D8N)
            return flex_mesh::ElementType::Q2D8N;
        else if (type == fea_io::ElementType::C2D6N)
            return flex_mesh::ElementType::C2D6N;
        // 3D Linear
        else if (type == fea_io::ElementType::T3D4N)
            return flex_mesh::ElementType::T3D4N;
        else if (type == fea_io::ElementType::P3D5N)
            return flex_mesh::ElementType::P3D5N;
        else if (type == fea_io::ElementType::W3D6N)
            return flex_mesh::ElementType::W3D6N;
        else if (type == fea_io::ElementType::H3D8N)
            return flex_mesh::ElementType::H3D8N;
        else if (type == fea_io::ElementType::Q3D4N)
            return flex_mesh::ElementType::Q3D4N;
        else if (type == fea_io::ElementType::C3D8N)
            return flex_mesh::ElementType::C3D8N;
        // 3D Quadratic
        else if (type == fea_io::ElementType::T3D10N)
            return flex_mesh::ElementType::T3D10N;
        else if (type == fea_io::ElementType::P3D13N)
            return flex_mesh::ElementType::P3D13N;
        else if (type == fea_io::ElementType::W3D15N)
            return flex_mesh::ElementType::W3D15N;
        else if (type == fea_io::ElementType::H3D20N)
            return flex_mesh::ElementType::H3D20N;
        else if (type == fea_io::ElementType::Q3D8N)
            return flex_mesh::ElementType::Q3D8N;
        else if (type == fea_io::ElementType::C3D16N)
            return flex_mesh::ElementType::C3D16N;
        throw _UTIL_MSG_EXCEPT("No conversion available.");
    }

    template<typename MeshType>
    inline fea_io::ElementType FEAIOInterface<MeshType>::convert_element_type(flex_mesh::ElementType type)
    {
        // 1D
        if (type == flex_mesh::ElementType::L1D2N)
            return fea_io::ElementType::L1D2N;
        else if (type == flex_mesh::ElementType::L1D3N)
            return fea_io::ElementType::L1D3N;
        // 2D
        else if (type == flex_mesh::ElementType::T2D3N)
            return fea_io::ElementType::T2D3N;
        else if (type == flex_mesh::ElementType::Q2D4N)
            return fea_io::ElementType::Q2D4N;
        else if (type == flex_mesh::ElementType::C2D4N)
            return fea_io::ElementType::C2D4N;
        else if (type == flex_mesh::ElementType::T2D6N)
            return fea_io::ElementType::T2D6N;
        else if (type == flex_mesh::ElementType::Q2D8N)
            return fea_io::ElementType::Q2D8N;
        else if (type == flex_mesh::ElementType::C2D6N)
            return fea_io::ElementType::C2D6N;
        // 3D Linear
        else if (type == flex_mesh::ElementType::T3D4N)
            return fea_io::ElementType::T3D4N;
        else if (type == flex_mesh::ElementType::P3D5N)
            return fea_io::ElementType::P3D5N;
        else if (type == flex_mesh::ElementType::W3D6N)
            return fea_io::ElementType::W3D6N;
        else if (type == flex_mesh::ElementType::H3D8N)
            return fea_io::ElementType::H3D8N;
        else if (type == flex_mesh::ElementType::Q3D4N)
            return fea_io::ElementType::Q3D4N;
        else if (type == flex_mesh::ElementType::C3D8N)
            return fea_io::ElementType::C3D8N;
        // 3D Quadratic
        else if (type == flex_mesh::ElementType::T3D10N)
            return fea_io::ElementType::T3D10N;
        else if (type == flex_mesh::ElementType::P3D13N)
            return fea_io::ElementType::P3D13N;
        else if (type == flex_mesh::ElementType::W3D15N)
            return fea_io::ElementType::W3D15N;
        else if (type == flex_mesh::ElementType::H3D20N)
            return fea_io::ElementType::H3D20N;
        else if (type == flex_mesh::ElementType::Q3D8N)
            return fea_io::ElementType::Q3D8N;
        else if (type == flex_mesh::ElementType::C3D16N)
            return fea_io::ElementType::C3D16N;
        else
            throw _UTIL_MSG_EXCEPT("Failed to converge element type to fea_io type.");
    }

    template<typename MeshType>
    template<typename NodeContainer, typename ElementContainer>
    inline typename flex_mesh::FEAIOInterface<MeshType>::Partition flex_mesh::FEAIOInterface<MeshType>::get_fea_io_partition(
        NodeContainer & nodes_in_partition, const ElementContainer & elements_in_partition, const PartitionIndex partition_id)
    {
        fea_io::Partition<GlobalIndex, LocalIndex, FloatType> partition;
        partition.partition_id = partition_id;

        // Populate node information
        LocalIndex local_num = 0;
        for (auto& n_ptr : nodes_in_partition)
        {
            n_ptr->set_local_number(local_num);
            partition.node_global_ids.push_back(n_ptr->get_global_number());
            partition.node_parent_partitions.push_back(n_ptr->get_parent_partition());
            for (auto i = 0; i < 3; i++)
            {
                partition.node_coords.push_back((*n_ptr)[i]);
            }
            local_num++;
        }

        // Populate element information
        for (auto& e_ptr : elements_in_partition)
        {
            partition.element_global_ids.push_back(e_ptr->get_global_number());
            partition.element_parent_partitions.push_back(e_ptr->get_parent_partition());
            partition.element_types.push_back(convert_element_type(e_ptr->get_type()));
            std::vector<LocalIndex> local_connectivity;
            for (auto& n_ptr : *e_ptr)
            {
                local_connectivity.push_back(n_ptr->get_local_number());
            }
            partition.element_connectivities.push_back(local_connectivity);
        }

        // Find what components are in this partition
        std::set<ComponentIndex> component_ids_in_partition;
        for (auto& e_ptr : elements_in_partition)
        {
            component_ids_in_partition.insert(e_ptr->get_component_id());
        }

        // Figure out what each subcomponent contains in this partition
        std::vector<fea_io::Subcomponent<LocalIndex, FloatType>> subcomponents;
        for (auto& component_id : component_ids_in_partition)
        {
            auto nodes_in_subcomponent_pair = find_nodes_in_component<NodeTypename>(elements_in_partition, component_id);
            auto elements_in_subcomponent_pair = find_elements_in_component<ElementTypename>(elements_in_partition, component_id);
            std::vector<NodePtr>& nodes_in_subcomponent = nodes_in_subcomponent_pair.second;
            std::vector<ElementPtr>& elements_in_subcomponent = elements_in_subcomponent_pair.second;

            fea_io::Subcomponent<LocalIndex, FloatType> subcomponent;
            subcomponent.component_id = component_id;

            for (auto& i : nodes_in_subcomponent_pair.first)
            {
                subcomponent.nodes_in_subcomponent.push_back(i);
            }
            for (auto& i : elements_in_subcomponent_pair.first)
            {
                subcomponent.elements_in_subcomponent.push_back(i);
            }
            for (auto& e_ptr : elements_in_subcomponent)
            {
                auto& lcs = e_ptr->get_lcs();
                for (auto i = 0; i < lcs.size(); i++)
                {
                    const Quaternion& quatern = lcs[i];
                    subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.w());
                    subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.x());
                    subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.y());
                    subcomponent.lcs_parameters_in_subcomponent.push_back(quatern.z());
                }
            }
            subcomponents.push_back(subcomponent);
        }
        partition.subcomponents = subcomponents;

        return partition;
    }

    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::add_ghost_layers(const LocalIndex n_ghost_levels, const ElementAdjacencyList<ElementTypename> adjacency_list, std::set<NodeTypename>& nodes, std::set<ElementTypename>& elements)
    {
        std::vector<GlobalIndex> elements_included_in_set;
        for (auto& e_ptr : elements)
        {
            elements_included_in_set.push_back(e_ptr->get_global_number());
        }
        adjacency_list.skip_in_neighbor_search(elements_included_in_set);

        for (auto level : utility::range_from_zero(n_ghost_levels))
        {
            std::set<ElementTypename> ghost_layer_elements;
            for (auto& e_ptr : elements)
            {
                for (auto& neighbor_ptr : adjacency_list.get_neighbors(e_ptr->get_global_number()))
                {
                    ghost_layer_elements.insert(neighbor_ptr);
                }
            }

            for (auto& ghost_element : ghost_layer_elements)
            {
                elements_included_in_set.push_back(ghost_element->get_global_number());
                elements.insert(ghost_element);
                for (auto& n_ptr : ghost_element)
                {
                    nodes.insert(n_ptr);
                }
            }
        }
    }

    template<typename MeshType>
    inline void FEAIOInterface<MeshType>::write_kernel(MeshType & partition_mesh, const PartitionIndex & partition_id, fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer)
    {
        auto partition = get_fea_io_partition(partition_mesh.nodes(), partition_mesh.elements(), partition_id);
        mesh_writer.write_partition(partition);
    }
}

#endif