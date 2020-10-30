#pragma once

#include "pfecBaseWriter.h"
#include "pfecWriterInterfaceObjects.h"
#include "pfecUtility.h"

#include <iterator>

namespace PFECNativeInterface {

    /// TODO add description of format and class
    /// Limitations to note:
    /// - number of nodes or elemnts in a partition is limited by the size of an PFEC_INT
    /// - the global indexing type can be changed in the Defines.h
    /// - the number of partitions is limited by the size of an PFEC_INT
    class MeshWriter : public BaseWriter
    {
    public:

        MeshWriter();
        /// Destructor, which calls close to cleanup the stream and variables
        ~MeshWriter();
        /// Write the header information to the stream, which includes size information and two lookup tables to jump to the data for each partition.
        /// The lookup tables just point to the beginning of the file until update_lookup_tables is called, so be sure to call that method before
        /// calling clear or the destructor is called.
		/// Note: partitions split up the mesh in a non-physical way that makes it convienient for computation, but domains are physical
		/// divisions of the mesh that define different materails or regions that should be kept independant.  When averaging stresses,
		/// stresses inside the same domain are averaged while those between domains are not.
        /// Parameters:
        ///   numGlobalNodes - the global number of nodes in the entire mesh
        ///   numGlobalElems - the global number of elements in the entire mesh
        ///   numDims - the number of dimensions for the problem (1 = 1D, 2 = 2D, 3 = 3D), which must match the number of coordinate components.
		///   domains - the domain numbers that exist in the whole mesh
        ///   globalNumParts - the total number of paritions that the mesh is broken into, but this file may not contain all partitions
        ///   localPartIDs - the partition numbers that will be contained in this file
        void write_header_info(PFEC_IND numGlobalNodes, PFEC_IND numGlobalElems, PFEC_INT numDims, std::vector<PFEC_INT> domains,
            PFEC_INT globalNumParts, std::vector<PFEC_INT> localPartIDs);
        /// Writes the nodal information for a particular partition.  This includes the global indices and the coordinates of each node.
        /// Parameters:
        ///   partID - the partition number
        ///   nodeIndices - a vector of the global indices of the nodes
		///   partitions - the partition that owns this node (some will not be owned)
        ///   coords - a vector of the coordinates of the nodes, where each coordinate is a vector of length of the number of dimensions
        /// Returns the starting position of the stream when this data is writen, which should be used to update the lookup tables.
        void write_part_nodal(const PFEC_INT partID, const std::vector<PFEC_IND> &nodeIndices, const std::vector<PFEC_INT> &partitions, const std::vector< std::vector<double> > &coords);
		void write_part_nodal(const PFEC_INT partID, const NodeInfo & ni);
        /// Templated version of write_part_nodal that uses a generic iterator, allow you to avoid allocating new memory if you
        /// already store the nodal information in a custom class.
        template <typename NodeAdapterType>
        void write_part_nodal2(PFEC_INT partID,
							   NodeAdapterType begin,
							   NodeAdapterType end);
        /// Writes the element information for a particular partition.  This includes the global index, the element type, the number of nodes
        /// in the connectivity, and the actual connectivity.
        /// Parameters:
        ///   partID - the partition number
        ///   globalIndices - a vector of the global indices of the elements
        ///   elemGeomTypes - a vector of the element geometry types, which is an integer.  (See definition of element formulation types in Defines.h)
        ///   localConnectivity - a vector of the connectivity of the elements
        void write_part_element(const PFEC_INT partID,
			const std::vector<PFEC_IND> &globalIndices,
			const std::vector<PFEC_INT> &partitions,
			const std::vector<PFEC_INT> &domains,
			const std::vector<ElemGeomType> &elemGeomTypes,
			const std::vector< std::vector<PFEC_INT>> &localConnectivity);
		void write_part_element(const PFEC_INT partID,const ElementInfo & ei);
        /// Templated version of write_part_element that uses a generic iterator, allow you to avoid allocating new memory if you
        /// already store the element information in a custom class.
        template <typename ElementAdapterType>
        void write_part_element2(PFEC_INT partID,
								 ElementAdapterType begin,
								 ElementAdapterType end);
        /// Updates the lookup table, closes the stream, and resets all variables to the default values.
        virtual void close();

    protected:
        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

        // Positions to save
        std::streampos node_lookup_pos_;
        std::streampos elem_lookup_pos_;

        // Important mesh information to store
        PFEC_IND num_global_nodes_;
        PFEC_IND num_global_elems_;
        PFEC_INT num_dims_;
        PFEC_INT num_global_parts_;

        // Variables needed to track changes.  Save positions as write_part_nodal or write_part_element
        // is called, and use these positions to update the lookup tables at the beginning of the destructor
        // being called.
        std::vector<PFEC_INT> updated_node_part_ids_;
        std::vector<PFEC_INT> updated_elem_part_ids_;
        std::vector<std::streampos> updated_node_part_pos_;
        std::vector<std::streampos> updated_elem_part_pos_;

        /// Updates the lookup tables with the stream positions for the node or element
        /// information that has been written since opening the stream.
        void update_lookup_tables();
        /// Updates one of the lookup tables
        void update_lookup_table(std::streampos tableStart,
                                 std::vector<PFEC_INT> &updatedIDs,
                                 std::vector<std::streampos> &updatedPos);
        /// Checks some common things for the methods in the class
        void common_error_checks(PFEC_INT partID, std::vector<PFEC_INT> &updateVec);
    };

    //////////////////////////////////////
    // Definitions of templated methods //
    //////////////////////////////////////

	template<typename NodeAdapterType>
	inline void PFECNativeInterface::MeshWriter::write_part_nodal2(PFEC_INT partID, NodeAdapterType begin, NodeAdapterType end)
	{
        common_error_checks(partID, updated_node_part_ids_);

        // Save the position for the lookup table
        std::streampos startPos = out_stream_.tellp();

        // Track how many were written
        PFEC_INT numNodesInPart = 0;

        // Iterate over the nodes and write the data
        out_stream_.write((char*)&partID, sizeof(PFEC_INT));
        std::streampos numNodesPos = out_stream_.tellp(); // Save position to update this later
        out_stream_.write((char*)&numNodesInPart, sizeof(PFEC_INT));
        for (auto it = begin; it != end; ++it) {
            PFEC_IND ind = it.get_global_ind();
			PFEC_INT owningPart = it.get_owning_partition();
            PFEC_INT len = it.get_coord_len();
            const double* coordStart = it.get_coords();
            out_stream_.write((char*)&ind, sizeof(PFEC_IND));
			out_stream_.write((char*)&owningPart, sizeof(PFEC_INT));
            out_stream_.write((char*)coordStart, len*sizeof(double));
            numNodesInPart++;
        }

        // Update how many nodes are in this partition
        std::streampos endPos = out_stream_.tellp();
        out_stream_.seekp(numNodesPos);
        out_stream_.write((char*)&numNodesInPart, sizeof(PFEC_INT));
        out_stream_.seekp(endPos);

        // Update the change log
        updated_node_part_ids_.push_back(partID);
        updated_node_part_pos_.push_back(startPos);
    }

	template<typename ElementAdapterType>
	inline void PFECNativeInterface::MeshWriter::write_part_element2(PFEC_INT partID, ElementAdapterType begin, ElementAdapterType end)
	{
		common_error_checks(partID, updated_elem_part_ids_);

        // Save the position for the lookup table
        std::streampos startPos = out_stream_.tellp();

        // Track how many were written
        PFEC_INT numElemsInPart = 0;

        // Write the data
        out_stream_.write((char*)&partID, sizeof(PFEC_INT));
        std::streampos numElemsPos = out_stream_.tellp(); // Save position to update this later
        out_stream_.write((char*)&numElemsInPart, sizeof(PFEC_INT));
        for (auto it = begin; it != end; ++it) {
            PFEC_IND ind = it.get_global_ind();
            out_stream_.write((char*)&ind, sizeof(PFEC_IND));
			PFEC_INT owningPart = it.get_owning_partition();
			out_stream_.write((char*)&owningPart, sizeof(PFEC_INT));
			PFEC_INT domainID = it.get_domain_id();
			out_stream_.write((char*)&domainID, sizeof(PFEC_INT));
            ElemGeomType gt = it.get_geom_type();
            out_stream_.write((char*)&gt, sizeof(ElemGeomType));
            PFEC_INT numNodeInConn = it.get_local_conn_len();
            PFEC_INT * localConn = it.get_local_conn();
            out_stream_.write((char*)localConn, numNodeInConn*sizeof(PFEC_INT));
            numElemsInPart++;
        }

        // Update how many nodes are in this partition
        std::streampos endPos = out_stream_.tellp();
        out_stream_.seekp(numElemsPos);
        out_stream_.write((char*)&numElemsInPart, sizeof(PFEC_INT));
        out_stream_.seekp(endPos);

        // Update the change log
        updated_elem_part_ids_.push_back(partID);
        updated_elem_part_pos_.push_back(startPos);
    }
};