#include "pfecBaseWriter.h"

#include <iostream>

PFECNativeInterface::BaseWriter::BaseWriter()
  : BaseIO()
{
}

void PFECNativeInterface::BaseWriter::open(std::string filePath, bool appendMode)
{
    if (appendMode)
        out_stream_.open(filePath, std::ios::binary | std::fstream::in | std::fstream::out);
    else
        out_stream_.open(filePath, std::ios::binary | std::fstream::out);
    if (!out_stream_.is_open())
        throw _MSG_EXCEPT_C(ERR_COULDNT_OPEN_STREAM);
}

void PFECNativeInterface::BaseWriter::close()
{
    out_stream_.close();
    BaseIO::close();
}

void PFECNativeInterface::BaseWriter::write_basic_header()
{
    // Write the version of file format that this file follows
    out_stream_.write((char*)&PFEC_FORMAT_VER, sizeof(PFEC_INT));
    // std::cout << "Version " << PFEC_FORMAT_VER << std::endl;
}

std::streampos PFECNativeInterface::BaseWriter::write_part_lookup_table(
  std::vector<PFEC_INT>& localPartIDs,
  std::vector<std::streamoff>& offsets)
{
    if (offsets.size() != localPartIDs.size() && offsets.size() != 0)
        throw _MSG_EXCEPT_C(ERR_PARAMETERS_DIFF_SIZE);

    std::streampos defaultPos = -1;
    std::streampos tableStart = out_stream_.tellp();
    // std::cout << "Writing part start at " << tableStart << std::endl;
    for (PFEC_INT i = 0; i < localPartIDs.size(); i++)
    {
        out_stream_.write((char*)&(localPartIDs[i]), sizeof(PFEC_INT));
        if (offsets.size() != 0)
            out_stream_.write((char*)&(offsets[i]), sizeof(std::streamoff));
        else
            out_stream_.write((char*)&(defaultPos), sizeof(std::streamoff));
    }
    return tableStart;
}

void PFECNativeInterface::BaseWriter::write_string(const std::string& s)
{
    PFEC_INT nameSize = static_cast<PFEC_INT>(s.size());
    out_stream_.write((char*)&nameSize, sizeof(PFEC_INT));
    out_stream_.write((char*)&(s[0]), nameSize * sizeof(char));
}
