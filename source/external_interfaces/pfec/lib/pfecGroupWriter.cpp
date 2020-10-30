#include "pfecGroupWriter.h"

PFECNativeInterface::GroupWriter::GroupWriter() : BaseWriter()
{
    header_written_ = false;
    node_group_pos_ = 0;
    elem_group_pos_ = 0;
}

void PFECNativeInterface::GroupWriter::write_header_info(std::vector<PFEC_INT>& localPartIDs, std::vector<std::string>& nodeGroupNames, std::vector<std::string>& elemGroupNames)
{
    if (header_written_)
        _MSG_EXCEPT_C(ERR_HEADER_ALREADY_WRITTEN);
    if (!out_stream_.is_open())
        _MSG_EXCEPT_C(ERR_STREAM_ALREADY_CLOSED);

    BaseWriter::write_basic_header();

    num_local_parts_ = static_cast<PFEC_INT>(localPartIDs.size());
    local_part_ids_ = localPartIDs;
    PFEC_INT numNodeGroups = static_cast<PFEC_INT>(nodeGroupNames.size());
	node_group_names_ = nodeGroupNames;
    PFEC_INT numElemGroups = static_cast<PFEC_INT>(elemGroupNames.size());
	elem_group_names_ = elemGroupNames;

    // Write number of parts
    out_stream_.write((char*)&num_local_parts_, sizeof(num_local_parts_));
    // Write number of node groups
    out_stream_.write((char*)&numNodeGroups, sizeof(numNodeGroups));
    // Write number of element groups
    out_stream_.write((char*)&numElemGroups, sizeof(numElemGroups));

    // Write partition id list
    out_stream_.write((char*)&(local_part_ids_[0]), num_local_parts_*sizeof(PFEC_INT));
    for (PFEC_INT i = 0; i < num_local_parts_; i++)
        part_id_map_[local_part_ids_[i]] = i;

    // Write group names
    write_names(nodeGroupNames);
    write_names(elemGroupNames);
    for (PFEC_INT i = 0; i < numNodeGroups; i++)
		node_group_name_to_index_map_[node_group_names_[i]] = i;
    for (PFEC_INT i = 0; i < numElemGroups; i++)
		elem_group_name_to_index_map_[elem_group_names_[i]] = i;

    // Write empty lookup tables
    node_group_pos_ = out_stream_.tellp();
    node_group_offsets_.resize(numNodeGroups, std::vector<std::streamoff>(num_local_parts_, -1));
    write_lookup_table(node_group_offsets_);
    elem_group_pos_ = out_stream_.tellp();
    elem_group_offsets_.resize(numElemGroups, std::vector<std::streamoff>(num_local_parts_, -1));
    write_lookup_table(elem_group_offsets_);

    header_end_ = out_stream_.tellp();
}

void PFECNativeInterface::GroupWriter::write_part_node_group(const PFEC_INT partID, const std::string nodeGroupName, const std::vector<PFEC_IND>& indices, const bool localNumbering)
{
    if (!header_written_)
        _MSG_EXCEPT_C(ERR_HEADER_NOT_WRITTEN);
    if (!out_stream_.is_open())
        _MSG_EXCEPT_C(ERR_STREAM_ALREADY_CLOSED);
    if (part_id_map_.count(partID) == 0)
        _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    if (node_group_name_to_index_map_.count(nodeGroupName) == 0)
        _MSG_EXCEPT_C(ERR_GROUP_ID_OUT_OF_BOUND);
    if (node_group_offsets_[node_group_name_to_index_map_[nodeGroupName]][part_id_map_[partID]] >= 0)
        _MSG_EXCEPT_C(ERR_PART_ALREADY_WRITTEN);

    node_group_offsets_[node_group_name_to_index_map_[nodeGroupName]][part_id_map_[partID]] = out_stream_.tellp() - header_end_;
    out_stream_.write((char*)&partID, sizeof(partID));
	write_string(nodeGroupName);
    out_stream_.write((char*)&localNumbering, sizeof(localNumbering));
    PFEC_IND numNodesInGroup = indices.size();
    out_stream_.write((char*)&numNodesInGroup, sizeof(PFEC_IND));
    out_stream_.write((char*)&indices[0], numNodesInGroup*sizeof(PFEC_IND));
}

void PFECNativeInterface::GroupWriter::write_part_elem_group(const PFEC_INT partID, const std::string elemGroupName, const std::vector<PFEC_IND>& indices, const bool localNumbering)
{
    if (!header_written_)
        _MSG_EXCEPT_C(ERR_HEADER_NOT_WRITTEN);
    if (!out_stream_.is_open())
        _MSG_EXCEPT_C(ERR_STREAM_ALREADY_CLOSED);
    if (part_id_map_.count(partID) == 0)
        _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    if (elem_group_name_to_index_map_.count(elemGroupName) == 0)
        _MSG_EXCEPT_C(ERR_GROUP_ID_OUT_OF_BOUND);
    if (elem_group_offsets_[elem_group_name_to_index_map_[elemGroupName]][part_id_map_[partID]] >= 0)
        _MSG_EXCEPT_C(ERR_PART_ALREADY_WRITTEN);

    elem_group_offsets_[elem_group_name_to_index_map_[elemGroupName]][part_id_map_[partID]] = out_stream_.tellp() - header_end_;
    out_stream_.write((char*)&partID, sizeof(partID));
	write_string(elemGroupName);
    out_stream_.write((char*)&localNumbering, sizeof(localNumbering));
    PFEC_IND numElemsInGroup = indices.size();
    out_stream_.write((char*)&numElemsInGroup, sizeof(PFEC_IND));
    out_stream_.write((char*)&indices[0], numElemsInGroup*sizeof(PFEC_IND));
}

void PFECNativeInterface::GroupWriter::close()
{
    // Update lookup tables
    if (out_stream_.is_open())
        update_lookup_tables();

    // Reset variables
    header_written_ = false;
    node_group_names_.clear();
    elem_group_names_.clear();
    node_group_name_to_index_map_.clear();
    elem_group_name_to_index_map_.clear();
    node_group_offsets_.clear();
    elem_group_offsets_.clear();
    node_group_pos_ = 0;
    elem_group_pos_ = 0;

    BaseWriter::close();
}

void PFECNativeInterface::GroupWriter::write_lookup_table(std::vector<std::vector<std::streamoff>>& offsets)
{
    for (PFEC_INT i = 0; i < offsets.size(); i++)
        out_stream_.write((char*)&(offsets[i][0]), offsets[i].size()*sizeof(std::streamoff));
}

void PFECNativeInterface::GroupWriter::write_names(std::vector<std::string>& names)
{
    PFEC_INT numEntries = static_cast<PFEC_INT>(names.size());
    for (PFEC_INT i = 0; i < numEntries; i++) {
        write_string(names[i]);
    }
}

void PFECNativeInterface::GroupWriter::update_lookup_tables()
{
    out_stream_.seekp(node_group_pos_);
    write_lookup_table(node_group_offsets_);
    out_stream_.seekp(elem_group_pos_);
    write_lookup_table(elem_group_offsets_);
}
