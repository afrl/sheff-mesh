#include "pfecMeshWriter.h"
#include "pfecUtility.h"

using namespace PFECNativeInterface;

PFECNativeInterface::MeshWriter::MeshWriter() : BaseWriter()
{
}

PFECNativeInterface::MeshWriter::~MeshWriter()
{
    close();
}

void PFECNativeInterface::MeshWriter::write_header_info(PFEC_IND numGlobalNodes, PFEC_IND numGlobalElems,
									   PFEC_INT numDims, std::vector<PFEC_INT> domains,
                                       PFEC_INT globalNumParts, std::vector<PFEC_INT> localPartNums)
{
    write_basic_header();

    // Store the important info
    num_global_nodes_ = numGlobalNodes;
    num_global_elems_ = numGlobalElems;
    num_dims_ = numDims;
    num_global_parts_ = globalNumParts;
    local_part_ids_ = localPartNums;
    num_local_parts_ = static_cast<PFEC_INT>(localPartNums.size());
    for (PFEC_INT i = 0; i < num_local_parts_; i++)
        part_id_map_[local_part_ids_[i]] = i;
    PFEC_INT numDomains = static_cast<PFEC_INT>(domains.size());

    // Write the basic information
    out_stream_.write((char*)&numGlobalNodes, sizeof(PFEC_IND));
    out_stream_.write((char*)&numGlobalElems, sizeof(PFEC_IND));
    out_stream_.write((char*)&numDims, sizeof(PFEC_INT));
	out_stream_.write((char*)&numDomains, sizeof(PFEC_INT));
	for (PFEC_INT i = 0; i < numDomains; i++)
		out_stream_.write((char*)&domains[i], sizeof(PFEC_INT));
    out_stream_.write((char*)&globalNumParts, sizeof(PFEC_INT));
    out_stream_.write((char*)&num_local_parts_, sizeof(PFEC_INT));

    // Write empty lookup tables and save the positions for later
	std::vector<std::streamoff> offsets;
    node_lookup_pos_ = write_part_lookup_table(localPartNums, offsets);
    elem_lookup_pos_ = write_part_lookup_table(localPartNums, offsets);

    header_end_ = out_stream_.tellp();
}

void PFECNativeInterface::MeshWriter::write_part_nodal(const PFEC_INT partID, const std::vector<PFEC_IND> &nodeIndices, const std::vector<PFEC_INT> &partitions, const std::vector< std::vector<double> > &coords)
{
    common_error_checks(partID, updated_node_part_ids_);
    // Check sizes are the same
    if (nodeIndices.size() != coords.size() || nodeIndices.size() != partitions.size())
        throw _MSG_EXCEPT_C(ERR_PARAMETERS_DIFF_SIZE);

    // Save the position for the lookup table
    std::streampos startPos = out_stream_.tellp();

    // Write the data
    PFEC_INT numNodes = static_cast<PFEC_INT>(nodeIndices.size());
    out_stream_.write((char*)&partID, sizeof(PFEC_INT));
    out_stream_.write((char*)&numNodes, sizeof(PFEC_INT));
    for (PFEC_INT i = 0; i < numNodes; i++) {
        out_stream_.write((char*)&(nodeIndices[i]), sizeof(PFEC_IND));
		out_stream_.write((char*)&(partitions[i]), sizeof(PFEC_INT));
        out_stream_.write((char*)&(coords[i][0]), num_dims_*sizeof(double));
    }

    // Update the change log
    updated_node_part_ids_.push_back(partID);
    updated_node_part_pos_.push_back(startPos);
}

void PFECNativeInterface::MeshWriter::write_part_nodal(const PFEC_INT partID, const NodeInfo & ni)
{
	this->write_part_nodal(partID, ni.nodeIndices, ni.partitions, ni.coords);
}

void PFECNativeInterface::MeshWriter::write_part_element(const PFEC_INT partID, const std::vector<PFEC_IND> &globalIndices, const std::vector<PFEC_INT> &partitions, const std::vector<PFEC_INT> &domains,
	const std::vector<ElemGeomType> &elemGeomTypes, const std::vector<std::vector<PFEC_INT>> &localConnectivity)
{
    common_error_checks(partID, updated_elem_part_ids_);
    // Check sizes are the same
	if ((globalIndices.size() != elemGeomTypes.size()) ||
		(globalIndices.size() != localConnectivity.size()) ||
		globalIndices.size() != partitions.size() ||
		globalIndices.size() != domains.size())
        throw _MSG_EXCEPT_C(ERR_PARAMETERS_DIFF_SIZE);

    // Save the position for the lookup table
    std::streampos startPos = out_stream_.tellp();

    // Write the data
    PFEC_INT numElems = static_cast<PFEC_INT>(globalIndices.size());
    out_stream_.write((char*)&partID, sizeof(PFEC_INT));
    out_stream_.write((char*)&numElems, sizeof(PFEC_INT));
    for (PFEC_INT i = 0; i < numElems; i++) {
        out_stream_.write((char*)&(globalIndices[i]), sizeof(PFEC_IND));
		out_stream_.write((char*)&(partitions[i]), sizeof(PFEC_INT));
		out_stream_.write((char*)&(domains[i]), sizeof(PFEC_INT));
        out_stream_.write((char*)&(elemGeomTypes[i]), sizeof(ElemGeomType));
        PFEC_INT numNodeInConn = static_cast<PFEC_INT>(localConnectivity[i].size());
        out_stream_.write((char*)&(localConnectivity[i][0]), numNodeInConn*sizeof(PFEC_INT));
    }

    // Update the change log
    updated_elem_part_ids_.push_back(partID);
    updated_elem_part_pos_.push_back(startPos);
}

void PFECNativeInterface::MeshWriter::write_part_element(const PFEC_INT partID, const ElementInfo & ei)
{
	write_part_element(partID,
		ei.globalIndices,
		ei.partitions,
		ei.domains,
		ei.elemGeomTypes,
		ei.localConnectivity);
}

void PFECNativeInterface::MeshWriter::close()
{
    // Update lookup tables
    if (out_stream_.is_open())
        update_lookup_tables();

    // Reset variables
    node_lookup_pos_ = 0;
    elem_lookup_pos_ = 0;
    num_global_nodes_ = 0;
    num_global_elems_ = 0;
    num_dims_ = 0;
    updated_node_part_ids_.clear();
    updated_elem_part_ids_.clear();
    updated_node_part_pos_.clear();
    updated_elem_part_pos_.clear();

    BaseWriter::close();
}

void PFECNativeInterface::MeshWriter::update_lookup_tables()
{
    update_lookup_table(node_lookup_pos_, updated_node_part_ids_, updated_node_part_pos_);
    update_lookup_table(elem_lookup_pos_, updated_elem_part_ids_, updated_elem_part_pos_);
}

void PFECNativeInterface::MeshWriter::update_lookup_table(std::streampos tableStart, std::vector<PFEC_INT>& updatedIDs, std::vector<std::streampos>& updatedPos)
{
    // Make a map of the part indices where the key is the part# and the value is the index in the other vectors
    std::map<PFEC_INT, PFEC_INT> updatedPartIDMap;
    for (PFEC_INT i = 0; i < updatedIDs.size(); i++) updatedPartIDMap[updatedIDs[i]] = i;

    // Seek to the nodal table and rewrite it
    out_stream_.seekp(tableStart);
    // Loop over the partitions in this file
    for (PFEC_INT i = 0; i < local_part_ids_.size(); i++) {
        PFEC_INT partID = local_part_ids_[i];
        // If the nodes for this partition was updated, then update the table
        if (updatedPartIDMap.find(partID) != updatedPartIDMap.end()) {
            PFEC_INT indInVec = updatedPartIDMap[i];
            out_stream_.write((char*)&(updatedIDs[indInVec]), sizeof(PFEC_INT));
            // Compute offset from header and write that to the file
            std::streamoff so = updatedPos[indInVec] - header_end_;
            out_stream_.write((char*)&(so), sizeof(std::streamoff));
        }
    }
}

void PFECNativeInterface::MeshWriter::common_error_checks(PFEC_INT partID, std::vector<PFEC_INT>& updateVec)
{
    // Check part is within bounds
	if (!in_vec(BaseIO::local_part_ids_, partID)) // This throws a compiler error for some reason
        throw _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    // Make sure the part hasn't been updated yet
    if (in_vec(updateVec, partID))
        throw _MSG_EXCEPT_C(ERR_PART_ALREADY_WRITTEN);
    // Check stream is still open
    if (!out_stream_.is_open())
        throw _MSG_EXCEPT_C(ERR_STREAM_ALREADY_CLOSED);
}
