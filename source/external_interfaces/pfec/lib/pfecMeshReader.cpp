#include "pfecMeshReader.h"

using namespace PFECNativeInterface;

PFECNativeInterface::MeshReader::MeshReader() : BaseReader()
{
    num_global_nodes_ = 0;
    num_global_elems_ = 0;
    num_dims_ = 0;
    num_global_parts_ = 0;
}
void PFECNativeInterface::MeshReader::get_part_sizes(const PFEC_INT partID, PFEC_INT & numNodes, PFEC_INT & numElements, PFEC_INT & numDims)
{
    common_error_checks(partID);

    move_to_partition(partID, node_lookup_table_);
    in_stream_.read((char*)&numNodes, sizeof(PFEC_INT));

    move_to_partition(partID, elem_lookup_table_);
    in_stream_.read((char*)&numElements, sizeof(PFEC_INT));

    numDims = num_dims_;
}
void MeshReader::read_part_nodal(PFEC_INT partID, std::vector<PFEC_IND>& nodeIndices, std::vector<PFEC_INT> &partitions, std::vector<std::vector<double>>& coords)
{
    common_error_checks(partID);

    move_to_partition(partID, node_lookup_table_);

    // Clear the vectors, in case there was previous data left in them
    nodeIndices.clear();
	partitions.clear();
    coords.clear();

    // Read some of the header info and resize the vectors accordingly data
    PFEC_INT numNodes = 0;
    in_stream_.read((char*)&numNodes, sizeof(PFEC_INT));
    nodeIndices.resize(numNodes);
	partitions.resize(numNodes);
    coords.resize(numNodes, std::vector<double>(num_dims_, 0.0));

    // Read the data
    for (PFEC_INT i = 0; i < numNodes; i++) {
        in_stream_.read((char*)&(nodeIndices[i]), sizeof(PFEC_IND));
		in_stream_.read((char*)&(partitions[i]), sizeof(PFEC_INT));
        in_stream_.read((char*)&(coords[i][0]), num_dims_*sizeof(double));
    }
}

NodeInfo PFECNativeInterface::MeshReader::read_part_nodal(PFEC_INT partID)
{
	NodeInfo ni;
	read_part_nodal(partID, ni.nodeIndices, ni.partitions, ni.coords);
	return ni;
}

void PFECNativeInterface::MeshReader::read_part_element(PFEC_INT partID, std::vector<PFEC_IND>& globalIndices, std::vector<PFEC_INT>& partitions, std::vector<PFEC_INT>& domains, std::vector<ElemGeomType>& elemGeomTypes, std::vector<std::vector<PFEC_INT>>& localConnectivity)
{
	common_error_checks(partID);

	move_to_partition(partID, elem_lookup_table_);

	// Clear the vectors, in case there was previous data left in them
	globalIndices.clear();
	partitions.clear();
	elemGeomTypes.clear();
	localConnectivity.clear();

	// Read some of the header info and resize the vectors accordingly data
	PFEC_INT numElems = 0;
	in_stream_.read((char*)&numElems, sizeof(PFEC_INT));
	globalIndices.resize(numElems);
	partitions.resize(numElems);
	domains.resize(numElems);
	elemGeomTypes.resize(numElems);
	localConnectivity.resize(numElems, std::vector<PFEC_INT>());

	// Get the map that relates the type to the number of nodes
	std::map<ElemGeomType, PFEC_INT> connMap = get_conn_num_node_map();

	// Read the data
	for (PFEC_INT i = 0; i < numElems; i++) {
		in_stream_.read((char*)&(globalIndices[i]), sizeof(PFEC_IND));
		in_stream_.read((char*)&(partitions[i]), sizeof(PFEC_INT));
		in_stream_.read((char*)&(domains[i]), sizeof(PFEC_INT));
		in_stream_.read((char*)&(elemGeomTypes[i]), sizeof(ElemGeomType));
		PFEC_INT numNodesInConn = connMap[elemGeomTypes[i]];
		localConnectivity[i].resize(numNodesInConn, 0);
		in_stream_.read((char*)&(localConnectivity[i][0]), numNodesInConn * sizeof(PFEC_INT));
	}
}

ElementInfo PFECNativeInterface::MeshReader::read_part_element(PFEC_INT partID)
{
	ElementInfo ei;
	read_part_element(partID,
		ei.globalIndices,
		ei.partitions,
		ei.domains,
		ei.elemGeomTypes,
		ei.localConnectivity);
	return ei;
}

void MeshReader::close()
{
    BaseReader::close();
    node_lookup_table_.clear();
    elem_lookup_table_.clear();
    num_global_nodes_ = 0;
    num_global_elems_ = 0;
    num_dims_ = 0;
    num_global_parts_ = 0;
}

void PFECNativeInterface::MeshReader::common_error_checks(PFEC_INT partID)
{
    // Check part is within bounds
    if (part_id_map_.count(partID) == 0)
        throw _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    // Check stream is still open
    if (!in_stream_.is_open())
        throw _MSG_EXCEPT_C(ERR_STREAM_ALREADY_CLOSED);
}

void MeshReader::read_header_info()
{
    BaseReader::read_header_info();

    // Read the basic information
    in_stream_.read((char*)&num_global_nodes_, sizeof(PFEC_IND));
    in_stream_.read((char*)&num_global_elems_, sizeof(PFEC_IND));
    in_stream_.read((char*)&num_dims_, sizeof(PFEC_INT));
	PFEC_INT numDomains = 0;
	in_stream_.read((char*)&numDomains, sizeof(PFEC_INT));
	domains_.resize(numDomains, 0);
	for (PFEC_INT i = 0; i < numDomains; i++)
		in_stream_.read((char*)&domains_[i], sizeof(PFEC_INT));
    in_stream_.read((char*)&num_global_parts_, sizeof(PFEC_INT));
    in_stream_.read((char*)&num_local_parts_, sizeof(PFEC_INT));

    // Read nodal and element lookup tables
    read_part_lookup_table(node_lookup_table_);
    read_part_lookup_table(elem_lookup_table_);
}