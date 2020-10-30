#include "pfecException.h"

#define _ERR_CASE_RETURN(x) return #x " occured in \n\t\t" + file + ":" + std::to_string(line) + " in \n\t\t" + ((cls == "") ? " function " + func : func + " method of " + cls);
#define _ERR_CASE(x) { case x: _ERR_CASE_RETURN(x) }

std::string PFECNativeInterface::MsgException::build_msg(const std::string file, const std::string cls, std::string func, int line)
{
    switch (code) {
        // Generic
        _ERR_CASE(ERR_NOT_IMPLEMENTED);
        // Index bounds
        _ERR_CASE(ERR_PART_ID_OUT_OF_BOUND);
        _ERR_CASE(ERR_NODE_ID_OUT_OF_BOUND);
        _ERR_CASE(ERR_ELEMENT_ID_OUT_OF_BOUND);
        _ERR_CASE(ERR_DATA_ID_OUT_OF_BOUND);
        _ERR_CASE(ERR_GROUP_ID_OUT_OF_BOUND);
        // Logic errors
        _ERR_CASE(ERR_COULDNT_OPEN_STREAM);
        _ERR_CASE(ERR_STREAM_ALREADY_CLOSED);
        _ERR_CASE(ERR_HEADER_NOT_WRITTEN);
        _ERR_CASE(ERR_HEADER_ALREADY_WRITTEN);
        _ERR_CASE(ERR_FILE_FORMAT_MISMATCH);
        _ERR_CASE(ERR_PART_ALREADY_WRITTEN);
		_ERR_CASE(ERR_PART_MISMATCH);
		_ERR_CASE(ERR_PART_WAS_NOT_WRITTEN);
        _ERR_CASE(ERR_LOOKUP_TABLE_CORRUPT);
        _ERR_CASE(ERR_META_FILE_ALREADY_HAS_ENTRY);
        _ERR_CASE(ERR_META_INVALID);
        _ERR_CASE(ERR_META_INVALID_TIMESTEP);
        _ERR_CASE(ERR_META_BAD_TIME_REF);
        _ERR_CASE(ERR_USED_WRONG_DATA_READ_METHOD);
        // Value errors
        _ERR_CASE(ERR_PARAMETERS_DIFF_SIZE);
        _ERR_CASE(ERR_PARAMETER_WRONG_SIZE);
        // Interface errors
        _ERR_CASE(ERR_INTERFACE);
		_ERR_CASE(ERR_NO_MPI);
    default:
        _ERR_CASE_RETURN(ERR_UNKNOWN)
    }
}
