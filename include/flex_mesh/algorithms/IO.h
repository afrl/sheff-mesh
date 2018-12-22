#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Exceptions.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"
#include "Search.h"

#ifdef SWIG
#define USE_FEAIO
#endif

#ifdef USE_FEAIO
#include "fea_io/meta/MetaWriter.h"
#include "fea_io/mesh/MeshWriter.h"
#include "fea_io/data/DataWriter.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>

namespace flex_mesh
{
	/// Write a VTK-type unstructured mesh
	template<typename MeshType>
	static void write_vtk(MeshType& mesh,
		const std::string& file_path,
		const bool write_partition_ids = true,
		const bool write_lcs = true);

#ifdef USE_FEAIO

	/// This collection of mesh algorithms pertains to reading/writing with the fea_io format
	template<typename MeshType>
	class FEAIOInterface
	{
	public:
		using NodePtr = typename MeshType::NodeTypename*;
		using ElementPtr = typename MeshType::ElementTypename*;
		using NodeDataCallback = std::function<std::vector<double>(std::string, NodePtr)>;
		using ElementDataCallback = std::function<std::vector<double>(std::string, ElementPtr)>;
		using NodeByElemDataCallback = std::function<std::vector<std::vector<double>>(std::string, ElementPtr)>;

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
			const std::string& base_file_path,
			const std::string& meta_name);

		static void write(
			MeshType& mesh,
			const std::string& base_file_path,
			const std::string& meta_name,
			const std::map<ComponentIndex, std::string>& component_names,
			std::vector<fea_io::ArrayInfo> array_info,
			NodeDataCallback node_data_callback,
			ElementDataCallback element_data_callback,
			NodeByElemDataCallback node_by_elem_data_callback);

		static void read_partition(MeshType& mesh,
			const std::string &MetaFileName,
			const int Partition = 0,
			const int TimeStep = 0);

		/// Computes the information for the mesh needed for writing with fea_Io
		static fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> get_fea_io_mesh_info(
			const MeshType& mesh, const std::map<ComponentIndex, std::string>& component_names);
		/// Computes the partition information needed for writing in parallel with fea_io
		/// Note: currently, this interface just writes in serial.
		static std::vector<fea_io::PartitionInfo<GlobalIndex, LocalIndex, FloatType>> get_fea_io_partition_infos(
			MeshType& mesh, const std::vector<PartitionIndex>& partition_ids);
		static fea_io::Partition<GlobalIndex, LocalIndex, FloatType> get_fea_io_partition(
			MeshType& mesh, const PartitionIndex partition_id);
	private:
		static flex_mesh::LCSAssignment convert_lcs_assignment(fea_io::LCSAssignment lcs_assignment);
		static fea_io::LCSAssignment convert_lcs_assignment(flex_mesh::LCSAssignment lcs_assignment);
		static flex_mesh::ElementType convert_element_type(fea_io::ElementType type);
		static fea_io::ElementType convert_element_type(flex_mesh::ElementType type);
	};
#endif


	/////////////////
	// Definitions //
	/////////////////

	template<typename MeshType>
	inline void write_vtk(MeshType & mesh, const std::string & file_path, const bool write_partition_ids, const bool write_lcs)
	{
		std::ofstream fout(file_path.c_str());
		if (!fout)
			throw _FM_MSG_EXCEPT("Could not create file " + file_path);

		mesh.update_global_numbering();

		fout << "# vtk DataFile Version 2.0" << std::endl;
		fout << "Mesh from BetaMesh utility" << std::endl;
		fout << "ASCII" << std::endl;
		fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
		fout << "POINTS " << mesh.nodes().size() << " double" << std::endl;

		for (auto& n_ptr : mesh.nodes()) {
			fout << n_ptr->get_coord()[0] << " " << n_ptr->get_coord()[1] << " " << n_ptr->get_coord()[2] << "" << std::endl;
		}

		int total_num_cell_conn_entries = 0;
		for (auto& e_ptr : mesh.elements()) {
			total_num_cell_conn_entries++;
			total_num_cell_conn_entries += e_ptr->get_num_nodes();
		}

		fout << "CELLS " << mesh.elements().size() << " " << total_num_cell_conn_entries << std::endl;

		for (auto& e_ptr : mesh.elements()) {
			auto vtk_permutation = get_vtk_permutation(e_ptr->get_type());
			if (vtk_permutation.size() > 0) {
				auto& connectivity = e_ptr->get_conn();
				// Note: connectivity may have uint8_t type, which will be interpretted as a char
				// by ostream, resulting in an ASCII character instead of the integer. "+0" forces
				// it to be interpretted as an integer.
				fout << connectivity.size() + 0;
				for (unsigned int i = 0; i < connectivity.size(); ++i) {
					fout << " " << connectivity[vtk_permutation[i]]->get_global_number();
				}
			}
			else {
				// How do we handle polygons and polyhedrons?
				// For now, skip unsupported types
				// Should there be a warning?
			}
			fout << std::endl;
		}

		fout << "CELL_TYPES " << mesh.elements().size() << std::endl;

		for (auto& e_ptr : mesh.elements()) {
			fout << get_vtk_element_type(e_ptr->get_type()) << std::endl;
		}

		// Write the component id's as cell data
		// Note: the cell data has to all be together
		fout << "CELL_DATA " << mesh.elements().size() << std::endl;
		fout << "SCALARS ComponentIDs int 1\n";
		fout << "LOOKUP_TABLE default\n";
		for (auto& e_ptr : mesh.elements()) {
			fout << e_ptr->get_component_id() << std::endl;
		}

		// Write orientation data
		if (write_lcs) {
			fout << "VECTORS LocalMaterialX float\n";
			for (auto& e_ptr : mesh.elements()) {
				Vector3D axis_in_global = e_ptr->get_lcs(0)._transformVector(Vector3D::Unit(0));
				for (auto i = 0; i < axis_in_global.size(); i++) fout << axis_in_global(i) << " ";
				fout << std::endl;
			}
			fout << "VECTORS LocalMaterialY float\n";
			for (auto& e_ptr : mesh.elements()) {
				Vector3D axis_in_global = e_ptr->get_lcs(0)._transformVector(Vector3D::Unit(1));
				for (auto i = 0; i < axis_in_global.size(); i++) fout << axis_in_global(i) << " ";
				fout << std::endl;
			}
			fout << "VECTORS LocalMaterialZ float\n";
			for (auto& e_ptr : mesh.elements()) {
				Vector3D axis_in_global = e_ptr->get_lcs(0)._transformVector(Vector3D::Unit(2));
				for (auto i = 0; i < axis_in_global.size(); i++) fout << axis_in_global(i) << " ";
				fout << std::endl;
			}
		}

		// Write the partition id's if asked for
		if (write_partition_ids) {
			fout << "SCALARS PartitionIDs int 1\n";
			fout << "LOOKUP_TABLE default\n";
			for (auto& e_ptr : mesh.elements()) {
				fout << e_ptr->get_parent_partition() << std::endl;
			}

			fout << "POINT_DATA " << mesh.nodes().size() << std::endl;
			fout << "SCALARS PartitionIDs int 1\n";
			fout << "LOOKUP_TABLE default\n";
			for (auto& n_ptr : mesh.nodes()) {
				fout << n_ptr->get_parent_partition() << std::endl;
			}
		}

		fout.close();
	}

#ifdef USE_FEAIO
	template<typename MeshType>
	inline void FEAIOInterface<MeshType>::write(MeshType & mesh, const std::string & base_file_path, const std::string& meta_name, const std::map<ComponentIndex, std::string>& component_names, std::vector<fea_io::ArrayInfo> array_info, NodeDataCallback node_data_callback, ElementDataCallback element_data_callback, NodeByElemDataCallback node_by_elem_data_callback)
	{
		// Get the number of partitions first
		PartitionIndex num_partitions = mesh.get_num_partitions();

		std::vector<PartitionIndex> partition_ids(num_partitions);
		// Fills vector with 0..num_partitions-1
		std::iota(partition_ids.begin(), partition_ids.end(), 0);

		// Start a new meta file
		fea_io::MetaWriter meta_writer;
		meta_writer.open(base_file_path, meta_name);
		std::string mesh_filepath = meta_writer.add_mesh_path(0, 0, meta_name);
		//writer.add_data_path(0, 0, meta_name);
		//writer.add_group_path(0, 0, meta_name);

		fea_io::MeshWriter<GlobalIndex, LocalIndex, FloatType> mesh_writer;
		auto mesh_info = get_fea_io_mesh_info(mesh, component_names);
		mesh_writer.initialize(mesh_filepath, mesh_info);

		for (auto& partition_id : partition_ids) {
			auto partition = get_fea_io_partition(mesh, partition_id);
			mesh_writer.write_partition(partition);
		}

		/*
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
					throw _BETA_MSG_EXCEPT2("Writing this type of data is not supported yet in the write_pfec_files method: " + dataNameMap[meta.type]);
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
						throw _BETA_MSG_EXCEPT2("Writing this type of data is not supported yet in the write_pfec_files method: " + dataNameMap[meta.type]);
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
			*/
		return;
	}

	template<typename MeshType>
	inline fea_io::MeshInfo<GlobalIndex, LocalIndex, FloatType> FEAIOInterface<MeshType>::get_fea_io_mesh_info(const MeshType & mesh, const std::map<ComponentIndex, std::string> & component_names)
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
		for (auto& component_id : component_id_vector) {
			auto elements_in_component = mesh.get_component(component_id);

			// Ensure the orientation type for all the elements is the same
			auto front_lcs_assignment = elements_in_component.front()->get_lcs_assignment();
			for (auto& e_ptr : elements_in_component) {
				if (e_ptr->get_lcs_assignment() != front_lcs_assignment) {
					throw _MSG_EXCEPT("The elements within the same component must have the same LCSAssignment to be written to the fea_io format.");
				}
			}

			// Figure out the name of the component
			auto comp_name_find_iter = component_names.find(component_id);
			std::string component_name;
			if (comp_name_find_iter == component_names.end())
				component_name = "Component" + std::to_string(component_id);
			else
				component_name = comp_name_find_iter->second;

			mesh_info.components.push_back({ component_id, component_name,
				convert_lcs_assignment(front_lcs_assignment), fea_io::LCSParameterization::Quaternion });
		}
		PartitionIndex num_partitions = mesh.get_num_partitions();
		mesh_info.total_num_partitions = num_partitions;
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
			for (auto& component_id : component_id_vector) {
				auto nodes_in_component = find_nodes_in_component<NodePtr>(elements_in_partition, component_id);
				auto elements_in_component = find_elements_in_component(elements_in_partition, component_id);

				subcomponent_infos.push_back({
					component_id,
					static_cast<LocalIndex>(nodes_in_component.size()),
					static_cast<LocalIndex>(elements_in_component.size())
					});
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
	inline fea_io::Partition<GlobalIndex, LocalIndex, FloatType> FEAIOInterface<MeshType>::get_fea_io_partition(MeshType & mesh, const PartitionIndex partition_id)
	{
		fea_io::Partition<GlobalIndex, LocalIndex, FloatType> partition;
		partition.partition_id = partition_id;

		// Update local indices
		auto numbering_result = mesh.update_local_numbering(partition_id);
		std::set<NodePtr>& nodes_in_partition = numbering_result.first;
		std::set<ElementPtr>& elements_in_partition = numbering_result.second;

		// Populate node information
		for (auto& n_ptr : nodes_in_partition) {
			partition.node_global_ids.push_back(n_ptr->get_global_number());
			partition.node_parent_partitions.push_back(n_ptr->get_parent_partition());
			for (auto i = 0; i < 3; i++)
				partition.node_coords.push_back((*n_ptr)[i]);
		}

		// Populate element information
		for (auto& e_ptr : elements_in_partition) {
			partition.element_global_ids.push_back(e_ptr->get_global_number());
			partition.element_parent_partitions.push_back(e_ptr->get_parent_partition());
			partition.element_types.push_back(convert_element_type(e_ptr->get_type()));
			partition.element_connectivities.push_back(e_ptr->get_local_conn());
		}

		// Find what components are in this partition
		std::set<ComponentIndex> component_ids_in_partition;
		for (auto& e_ptr : elements_in_partition)
			component_ids_in_partition.insert(e_ptr->get_component_id());

		// Figure out what each subcomponent contains in this partition
		std::vector<fea_io::Subcomponent<LocalIndex, FloatType>> subcomponents;
		for (auto& component_id : component_ids_in_partition) {
			auto nodes_in_subcomponent = find_nodes_in_component<NodePtr>(elements_in_partition, component_id);
			auto elements_in_subcomponent = find_elements_in_component(elements_in_partition, component_id);

			fea_io::Subcomponent<LocalIndex, FloatType> subcomponent;
			subcomponent.component_id = component_id;
			for (auto& n_ptr : nodes_in_subcomponent)
				subcomponent.nodes_in_subcomponent.push_back(n_ptr->get_local_number());
			for (auto& e_ptr : elements_in_subcomponent) {
				subcomponent.elements_in_subcomponent.push_back(e_ptr->get_local_number());

				auto& lcs = e_ptr->get_lcs();
				for (auto i = 0; i < lcs.size(); i++) {
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
			throw _MSG_EXCEPT("Cannot convert from fea_io::LCSAssignment.");
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
			throw _MSG_EXCEPT("Cannot convert from fea_io::LCSAssignment.");
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
			throw _MSG_EXCEPT("Failed to converge element type to fea_io type.");
	}

#endif

}