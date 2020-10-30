#pragma once

#include "pfecBaseIO.h"

namespace PFECNativeInterface
{

    /// A base class for all readers to inherit from (except the MetaReader)
    /// Contains some basic methods and properties needed for all readers
    class BaseReader : public BaseIO
    {
    public:

        BaseReader();
        /// Opens stream and reads any header information (via a
        /// call to reader_header_info).
        virtual void open(std::string filePath);
        /// Closes stream and cleans any internal variables (if overloaded
        /// be sure to call parent's close).  Each class will take care of
        /// its own variables in the inheritance hierarchy.
        virtual void close();
        /// Returns a vector of the partition IDs that is contained within this file.
        std::vector<PFEC_INT> get_parts_in_file() { return local_part_ids_; };

    protected:
        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

        /// The stream to file being read
        std::ifstream in_stream_;

        /// Just reads a string from the binary stream. Assumes the stream is in the correct
        /// position.
        std::string read_string();
        /// Reads a lookup table that is based on partition, which is expected
        /// in the format of:
        /// <PFEC_INT partitionID><streamoff offset>[Repeated for each local partition]
        /// Note that the offset is relative to the end of the lookup table, not
        /// the beginning of the file. (Assumes num_local_parts_ has been set)
        void read_part_lookup_table(std::vector<std::streamoff> &offsets);
        /// Moves the stream to the start of a particular partition and reads
        /// the partID to check if it is right place.
        void move_to_partition(PFEC_INT partID, std::vector<std::streamoff> &offsets);

        /// Called after stream is open within the public open method
        /// This method definition reads part of the header universal to all readers.
        /// Handle reading any header information, lookup tables, etc. in this method.
        /// Note: if you overload, be sure to call parent's method.
        virtual void read_header_info();
    };

};