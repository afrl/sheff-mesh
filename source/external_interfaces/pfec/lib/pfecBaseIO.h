#pragma once

#include "pfecDefines.h"
#include "pfecException.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace PFECNativeInterface
{

    /// An abstract class that all reader/writers inherit from
    /// It contains methods and properties common to all reader/writers
    class BaseIO
    {
    public:
        BaseIO();
        /// Destructor calls cloase to make sure the stream is cleaned up
        virtual ~BaseIO();
        /// Closes stream and cleans up any internal variables,
        /// allowing another open to be called
        virtual void close();

    protected:
        /// The number of partitions in the file
		PFEC_INT num_local_parts_;
        /// The ids of the partitions in this file
        std::vector<PFEC_INT> local_part_ids_;
        /// A convenient map of the partition ID to the index of the vectors
        std::map<PFEC_INT, PFEC_INT> part_id_map_;
        /// Stream position of the end of the header information, which is basis
        /// for the offsets contained in the header lookup tables. Set at end of open method.
        std::streampos header_end_;
    };

};