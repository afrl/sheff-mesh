#include "pfecDataReader.h"

using namespace PFECNativeInterface;

PFECNativeInterface::DataReader::DataReader() : BaseReader()
{
    time = 0.0;
    num_data_blocks_ = 0;
}

void PFECNativeInterface::DataReader::get_data_info(std::vector<PFEC_INT>& dataIDs, std::vector<DataType>& dataTypes, std::vector<PFEC_INT>& dataNumComp, std::vector<std::string>& dataNames, std::vector< std::vector<std::string> > &dataCompNames) const
{
    dataIDs = data_ids_;
    dataTypes = data_types_;
    dataNumComp = data_num_comp_;
    dataNames = data_names_;
    dataCompNames = data_comp_names_;
}

std::vector<DataMetaData> PFECNativeInterface::DataReader::get_data_info() const
{
	std::vector<DataMetaData> meta;
	for (auto i = 0; i < data_ids_.size(); i++) {
		meta.push_back(DataMetaData(data_ids_[i], data_types_[i], data_num_comp_[i], data_names_[i], data_comp_names_[i]));
	}
	return meta;
}

void PFECNativeInterface::DataReader::read_part_unique_data(PFEC_INT partID, PFEC_INT dataID, std::vector< std::vector<double> > &data)
{
    common_error_checks(partID, dataID);
    DataType thisDataType = data_types_[data_id_map_[dataID]];
    if (thisDataType != DataType::Nodal && thisDataType != DataType::Element)
        _MSG_EXCEPT_C(ERR_USED_WRONG_DATA_READ_METHOD);

    move_to_data_block(partID, dataID);

    // Read the length of the data and resize the vector
    unsigned numComps = data_num_comp_[data_id_map_[dataID]];
    PFEC_INT len = 0;
    in_stream_.read((char*)&len, sizeof(len));
    data.resize(len, std::vector<double>(numComps, 0.0));
    for (PFEC_INT i = 0; i < len; i++) {
        in_stream_.read((char*)&(data[i][0]), numComps*sizeof(double));
    }
}

void PFECNativeInterface::DataReader::read_part_by_elem_data(PFEC_INT partID, PFEC_INT dataID, std::vector< std::vector< std::vector<double> > >& data)
{
    common_error_checks(partID, dataID);
    DataType thisDataType = data_types_[data_id_map_[dataID]];
    if (thisDataType != DataType::NodalByElement && thisDataType != DataType::QuadPoint)
        _MSG_EXCEPT_C(ERR_USED_WRONG_DATA_READ_METHOD);

    move_to_data_block(partID, dataID);

    PFEC_INT num_elems = 0;
    in_stream_.read((char*)&num_elems, sizeof(num_elems));
    data.resize(num_elems);
    unsigned numComps = data_num_comp_[data_id_map_[dataID]];
    for (PFEC_INT i = 0; i < num_elems; i++) {
        // Read the number of nodes/quads in this element
        PFEC_INT len = 0;
        in_stream_.read((char*)&len, sizeof(len));
        data[i].resize(len, std::vector<double>(numComps, 0.0));
        // Read the data for this element
        for (PFEC_INT j = 0; j < len; j++)
            in_stream_.read((char*)&(data[i][j][0]), numComps*sizeof(double));
    }
}

void PFECNativeInterface::DataReader::read_part_orientation_data(PFEC_INT partID, PFEC_INT dataID, std::vector<std::vector<std::array<double, 4>>>& quaterions)
{
    common_error_checks(partID, dataID);
    DataType thisDataType = data_types_[data_id_map_[dataID]];
    if (thisDataType != DataType::Orientation)
        _MSG_EXCEPT_C(ERR_USED_WRONG_DATA_READ_METHOD);

    move_to_data_block(partID, dataID);

    PFEC_INT num_elems = 0;
    in_stream_.read((char*)&num_elems, sizeof(num_elems));
    quaterions.resize(num_elems);
    for (PFEC_INT i = 0; i < num_elems; i++) {
        // Read the number of nodes/quads in this element
        PFEC_INT len = 0;
        in_stream_.read((char*)&len, sizeof(len));
		// Read the number of rotations
		PFEC_INT numQuaterionTerms = 0;
		in_stream_.read((char*)&numQuaterionTerms, sizeof(numQuaterionTerms));
		quaterions[i].resize(len);
        // Read angles
        for (PFEC_INT j = 0; j < len; j++)
            in_stream_.read((char*)&(quaterions[i][j][0]), numQuaterionTerms *sizeof(double));
    }
}

std::ifstream * PFECNativeInterface::DataReader::get_part_data_stream(PFEC_INT partID, PFEC_INT dataID)
{
    throw _MSG_EXCEPT_C(ERR_NOT_IMPLEMENTED);
}

void PFECNativeInterface::DataReader::close()
{
    time = 0.0;
    num_data_blocks_ = 0;
    data_ids_.clear();
    data_id_map_.clear();
    data_types_.clear();
    data_num_comp_.clear();
    data_names_.clear();
    data_comp_names_.clear();
    local_part_offsets_.clear();
}

std::map<PFEC_INT, std::streamoff> DataReader::read_part_header(PFEC_INT partID)
{
    // Read the vector of offsets (which will be relative to the end of the
    // lookup table for the partition)
    std::map<PFEC_INT, std::streamoff> partDataOffsets;
    for (PFEC_INT i = 0; i < num_data_blocks_; i++) {
        PFEC_INT dataID = 0;
        std::streamoff dataOff = 0;
        in_stream_.read((char*)&dataID, sizeof(PFEC_INT));
        in_stream_.read((char*)&dataOff, sizeof(std::streamoff));
        partDataOffsets[dataID] = dataOff;
    }
    return partDataOffsets;
}

void PFECNativeInterface::DataReader::move_to_data_block(PFEC_INT partID, PFEC_INT dataID)
{
    move_to_partition(partID, local_part_offsets_);

    // Read the partition header to get the data offsets
    std::map<PFEC_INT, std::streamoff> dataOffsets = read_part_header(partID);
    std::streampos endOfDataTable = in_stream_.tellg();

    // Move the stream to the position of the data within the partition
    in_stream_.seekg(endOfDataTable + dataOffsets[dataID]);

    // Check the data ID in the file
    PFEC_INT dataIDInFile = 0;
    in_stream_.read((char*)&dataIDInFile, sizeof(PFEC_INT));
    if (dataIDInFile != dataID)
        throw _MSG_EXCEPT_C(ERR_LOOKUP_TABLE_CORRUPT);
}

void PFECNativeInterface::DataReader::common_error_checks(PFEC_INT partID, PFEC_INT dataID)
{
    // Check part is within bounds
    if (part_id_map_.count(partID) == 0)
        throw _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    // Check data ID is within bounds
    if (data_id_map_.count(dataID) == 0)
        throw _MSG_EXCEPT_C(ERR_DATA_ID_OUT_OF_BOUND);
    // Check stream is still open
    if (!in_stream_.is_open())
        throw _MSG_EXCEPT_C(ERR_STREAM_ALREADY_CLOSED);
}

void PFECNativeInterface::DataReader::read_header_info()
{
    BaseReader::read_header_info();

    // Read time
    in_stream_.read((char*)&time, sizeof(double));
	// Read num partitions
	in_stream_.read((char*)&num_local_parts_, sizeof(PFEC_INT));
    // Read data information
    in_stream_.read((char*)&num_data_blocks_, sizeof(PFEC_INT));
    data_ids_.resize(num_data_blocks_, 0);
    data_types_.resize(num_data_blocks_, DataType::Nodal);
    data_num_comp_.resize(num_data_blocks_, 0);
    data_names_.resize(num_data_blocks_, "");
    data_comp_names_.resize(num_data_blocks_);
    for (PFEC_INT i = 0; i < num_data_blocks_; i++) {
        in_stream_.read((char*)&(data_ids_[i]), sizeof(PFEC_INT));
        in_stream_.read((char*)&(data_types_[i]), sizeof(DataType));
        in_stream_.read((char*)&(data_num_comp_[i]), sizeof(PFEC_INT));
        data_names_[i] = read_string();
        // Read names of the components (only if the components are > 0)
        if (data_num_comp_[i] > 0)
            for (PFEC_INT j = 0; j < data_num_comp_[i]; j++)
                data_comp_names_[i].push_back(read_string());
    }
    for (PFEC_INT i = 0; i < num_data_blocks_; i++)
        data_id_map_[data_ids_[i]] = i;

    // Read the partition lookup table
    read_part_lookup_table(local_part_offsets_);
}
