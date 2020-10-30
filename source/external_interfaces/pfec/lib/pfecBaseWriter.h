#pragma once

#include "pfecBaseIO.h"

namespace PFECNativeInterface
{
    /// A base class for all writers to inherit from (except the MetaWriter)
    /// Contains some basic methods and properties needed for all writers
    class BaseWriter : public BaseIO
    {
    public:

        BaseWriter();
		virtual ~BaseWriter() {}

        /// Opens stream
        virtual void open(std::string filePath, bool appendMode = false);
        /// Closes the stream and cleansup any internal variables.
        /// If overloaded, be sure to call parent's method.
        virtual void close();

    protected:
        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

        /// The stream to file being written to
        std::fstream out_stream_;

        /// Writes any constant information that is common to all file types
        void write_basic_header();
        /// Writes an empty lookup table by partition and returns the position
        /// of the beginning of the table, so it can be updated later.
        std::streampos write_part_lookup_table(std::vector<PFEC_INT> &localPartIDs,
                                               std::vector<std::streamoff> &offsets);
        void write_string(const std::string &s);
    };
};