#include "pfecBaseReader.h"

using namespace PFECNativeInterface;

PFECNativeInterface::BaseReader::BaseReader() : BaseIO()
{
}

void PFECNativeInterface::BaseReader::open(std::string filePath)
{
    in_stream_.open(filePath, std::ios::binary);
    if (!in_stream_.is_open())
        throw _MSG_EXCEPT2_C(ERR_COULDNT_OPEN_STREAM, "Was trying to open " + filePath);
    read_header_info();
    header_end_ = in_stream_.tellg();
}

void PFECNativeInterface::BaseReader::close()
{
    in_stream_.close();
    BaseIO::close();
}

std::string PFECNativeInterface::BaseReader::read_string()
{
    PFEC_INT strSize = 0;
    in_stream_.read((char*)&strSize, sizeof(PFEC_INT));
    char* buff = new char[strSize + 1];
    in_stream_.read(buff, strSize*sizeof(char));
    buff[strSize] = '\0'; // Termination character
    std::string str(buff);
    delete[] buff;
    return str;
}

void PFECNativeInterface::BaseReader::read_part_lookup_table(std::vector<std::streamoff>& offsets)
{
    local_part_ids_.resize(num_local_parts_, 0);
    offsets.resize(num_local_parts_, 0);
    for (PFEC_INT i = 0; i < num_local_parts_; i++) {
        in_stream_.read((char*)&(local_part_ids_[i]), sizeof(PFEC_INT));
        in_stream_.read((char*)&(offsets[i]), sizeof(std::streamoff));
    }
    // Build the partition ID map to make it easy to get the index for partition ID
    for (PFEC_INT i = 0; i < num_local_parts_; i++)
        part_id_map_[local_part_ids_[i]] = i;
}

void PFECNativeInterface::BaseReader::move_to_partition(PFEC_INT partID, std::vector<std::streamoff> &offsets)
{
    // Move the stream to the correct partition
    PFEC_INT localPartInd = part_id_map_[partID];
    in_stream_.seekg(header_end_ + offsets[localPartInd]);
    // Read partition ID and chack we are in the right place
    PFEC_INT partIDInFile = 0;
    in_stream_.read((char*)&partIDInFile, sizeof(PFEC_INT));
    if (partIDInFile != partID)
        throw _MSG_EXCEPT_C(ERR_LOOKUP_TABLE_CORRUPT);
}

void PFECNativeInterface::BaseReader::read_header_info()
{
    // Read version of file format that this file follows
    PFEC_INT formatVer = -1;
    in_stream_.read((char*)&formatVer, sizeof(PFEC_INT));
    if (formatVer != PFEC_FORMAT_VER)
        throw _MSG_EXCEPT_C(ERR_FILE_FORMAT_MISMATCH);
}
