#pragma once

#include <string>
#include <exception>
#include <typeinfo>

namespace PFECNativeInterface {

    enum BFErrorCodes {
        // Generic
        ERR_UNKNOWN,
        ERR_NOT_IMPLEMENTED,
        // Index bounds
        ERR_PART_ID_OUT_OF_BOUND,
        ERR_NODE_ID_OUT_OF_BOUND,
        ERR_ELEMENT_ID_OUT_OF_BOUND,
        ERR_DATA_ID_OUT_OF_BOUND,
        ERR_GROUP_ID_OUT_OF_BOUND,
        // Logic errors
        ERR_COULDNT_OPEN_STREAM,
        ERR_STREAM_ALREADY_CLOSED,
        ERR_HEADER_NOT_WRITTEN,
        ERR_HEADER_ALREADY_WRITTEN,
        ERR_FILE_FORMAT_MISMATCH,
        ERR_PART_ALREADY_WRITTEN,
		ERR_PART_MISMATCH,
		ERR_PART_WAS_NOT_WRITTEN,
        ERR_LOOKUP_TABLE_CORRUPT,
        ERR_META_FILE_ALREADY_HAS_ENTRY,
        ERR_META_INVALID,
		ERR_META_INVALID_TIMESTEP,
        ERR_META_BAD_TIME_REF,
        ERR_USED_WRONG_DATA_READ_METHOD,
        // Value errors
        ERR_PARAMETERS_DIFF_SIZE,
        ERR_PARAMETER_WRONG_SIZE,
        // interface errors
        ERR_INTERFACE,
		ERR_NO_MPI,
    };

// Use this in functions
#define _MSG_EXCEPT(x) PFECNativeInterface::MsgException(x, __FILE__, __func__, __LINE__, "")
#define _MSG_EXCEPT2(x, m) PFECNativeInterface::MsgException(x, __FILE__, __func__, __LINE__, m)
// Use this in class methods
#define _MSG_EXCEPT_C(x) PFECNativeInterface::MsgException(x, __FILE__, typeid(*this).name(), __func__, __LINE__, "")
#define _MSG_EXCEPT2_C(x, m) PFECNativeInterface::MsgException(x, __FILE__, typeid(*this).name(), __func__, __LINE__, m)

    class MsgException : public std::exception
    {
    public:
        MsgException(const BFErrorCodes errcode, const std::string file, const std::string func, int line, const std::string comments) {
            code = errcode;
            msg = build_msg(file, "", func, line);
			msg += "\n\nComments: " + comments + "\n";
        }

        MsgException(const BFErrorCodes errcode, const std::string file, const std::string cls, std::string func, int line, const std::string comments) {
            code = errcode;
            msg = build_msg(file, cls, func, line);
			msg += "\n\nComments: " + comments + "\n";
        }

        virtual const char* what() const throw() {
            return msg.c_str();
        }

    private:
        BFErrorCodes code;
        std::string msg;

        std::string build_msg(const std::string file, const std::string cls,  std::string func, int line);
    };
}