#include "pfecGroupReader.h"

PFECNativeInterface::GroupReader::GroupReader() : BaseReader()
{
}

void PFECNativeInterface::GroupReader::get_group_info(std::vector<std::string>& nodeGroupNames, std::vector<std::string>& elemGroupNames)
{
    nodeGroupNames = node_group_names_;
    elemGroupNames = elem_group_names_;
}

void PFECNativeInterface::GroupReader::get_part_node_group(PFEC_INT partID, std::string nodeGroupName, std::vector<PFEC_IND>& indices, bool &localNumbering)
{
    move_to_group_info(partID, nodeGroupName, node_group_name_to_index_map_, node_group_offsets_);

    in_stream_.read((char*)&localNumbering, sizeof(localNumbering));
    PFEC_IND numNodesInGroup = 0;
    in_stream_.read((char*)&numNodesInGroup, sizeof(PFEC_IND));
    indices.resize(numNodesInGroup);
    in_stream_.read((char*)&indices[0], numNodesInGroup*sizeof(PFEC_IND));
}

void PFECNativeInterface::GroupReader::get_part_elem_group(PFEC_INT partID, std::string elemGroupName, std::vector<PFEC_IND>& indices, bool &localNumbering)
{
    move_to_group_info(partID, elemGroupName, elem_group_name_to_index_map_, elem_group_offsets_);

    in_stream_.read((char*)&localNumbering, sizeof(localNumbering));
    PFEC_IND numElemsInGroup = 0;
    in_stream_.read((char*)&numElemsInGroup, sizeof(PFEC_IND));
    indices.resize(numElemsInGroup);
    in_stream_.read((char*)&indices[0], numElemsInGroup*sizeof(PFEC_IND));
}

void PFECNativeInterface::GroupReader::read_header_info()
{
    BaseReader::read_header_info();

    // Read number of partitions in this file
    in_stream_.read((char*)&num_local_parts_, sizeof(num_local_parts_));
    // Read number of node groups
    PFEC_INT numNodeGroups = 0;
    in_stream_.read((char*)&numNodeGroups, sizeof(PFEC_INT));
    // Read number of element groups
    PFEC_INT numElemGroups = 0;
    in_stream_.read((char*)&numElemGroups, sizeof(PFEC_INT));

    // Read partition list
    local_part_ids_.resize(num_local_parts_);
    in_stream_.read((char*)&(local_part_ids_[0]), num_local_parts_*sizeof(PFEC_INT));
    for (PFEC_INT i = 0; i < num_local_parts_; i++)
        part_id_map_[local_part_ids_[i]] = i;

    // Read group names
    read_names(numNodeGroups, node_group_names_);
    read_names(numElemGroups, elem_group_names_);
    for (PFEC_INT i = 0; i < numNodeGroups; i++)
		node_group_name_to_index_map_[node_group_names_[i]] = i;
    for (PFEC_INT i = 0; i < numElemGroups; i++)
		elem_group_name_to_index_map_[elem_group_names_[i]] = i;

    
    // Read lookup table
    read_lookup_table(numNodeGroups, node_group_offsets_);
    read_lookup_table(numElemGroups, elem_group_offsets_);
}

void PFECNativeInterface::GroupReader::read_lookup_table(PFEC_INT numGroups, std::vector<std::vector<std::streamoff>>& offsets)
{
    offsets.resize(numGroups, std::vector<std::streamoff>(num_local_parts_, -1));
    for (PFEC_INT i = 0; i < numGroups; i++)
        in_stream_.read((char*)&(offsets[i][0]), num_local_parts_*sizeof(std::streamoff));
}

void PFECNativeInterface::GroupReader::read_names(PFEC_INT num_entries, std::vector<std::string>& names)
{
    names.resize(num_entries);
    for (PFEC_INT i = 0; i < num_entries; i++) {
        names[i] = read_string();
    }
}

void PFECNativeInterface::GroupReader::move_to_group_info(PFEC_INT partID, std::string groupName, std::map<std::string, PFEC_INT>& nameToIndexMap, std::vector<std::vector<std::streamoff>>& offsets)
{
    // Check part is within bounds
    if (part_id_map_.count(partID) == 0)
        throw _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    // Check data ID is within bounds
    if (nameToIndexMap.count(groupName) == 0)
        throw _MSG_EXCEPT_C(ERR_GROUP_ID_OUT_OF_BOUND);

    // Move to the position of the group info
    PFEC_INT partIndex = part_id_map_[partID];
    PFEC_INT groupIndex = nameToIndexMap[groupName];
    in_stream_.seekg(header_end_ + offsets[groupIndex][partIndex]);

    // Read partition ID and chack we are in the right place
    PFEC_INT partIDInFile = 0;
    in_stream_.read((char*)&partIDInFile, sizeof(PFEC_INT));
    if (partIDInFile != partID)
        throw _MSG_EXCEPT_C(ERR_LOOKUP_TABLE_CORRUPT);
    // Read group ID and chack we are in the right place
    PFEC_INT groupIDInFile = 0;
	std::string nameInFile = read_string();
    if (nameInFile != groupName)
        throw _MSG_EXCEPT_C(ERR_LOOKUP_TABLE_CORRUPT);
}
