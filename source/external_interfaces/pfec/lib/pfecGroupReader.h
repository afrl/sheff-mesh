#pragma once

#include "pfecBaseReader.h"

namespace PFECNativeInterface {

    /// TODO add description and example
    class GroupReader : public BaseReader
    {
    public:

        GroupReader();
        /// Gets the ids and names of the node and element groups
        void get_group_info(std::vector<std::string> &nodeGroupNames,
                            std::vector<std::string> &elemGroupNames);
        /// Get a list of the nodes in the respective partition and node group
        void get_part_node_group(const PFEC_INT partID, const std::string nodeGroupName, std::vector<PFEC_IND> &indices, bool &localNumbering);
        /// Get a list of the elements in the respective partition and element group
        void get_part_elem_group(const PFEC_INT partID, const std::string elemGroupName, std::vector<PFEC_IND> &indices, bool &localNumbering);

    private:
        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

        /// The names for the node groups
        std::vector<std::string> node_group_names_;
        /// The names for the element groups
        std::vector<std::string> elem_group_names_;
        // Convienient maps
		std::map<std::string, PFEC_INT> node_group_name_to_index_map_;
		std::map<std::string, PFEC_INT> elem_group_name_to_index_map_;
        /// Offsets in the stream for the node group data
        /// First index is group, second index is partition
        std::vector< std::vector<std::streamoff> > node_group_offsets_;
        /// Offsets in the stream for the node group data
        /// First index is group, second index is partition
        std::vector< std::vector<std::streamoff> > elem_group_offsets_;

        /// Assumes the stream is open and the get pointer is at the beginning of the file
        /// ready to parse the header information.  Only called when first opening a file.
        virtual void read_header_info();
        /// Reads the lookup table of offsets ordered by group then by partition
        void read_lookup_table(PFEC_INT numGroups,
                                    std::vector< std::vector<std::streamoff> > &offsets);
        /// A helper function used to read a list of strings
        void read_names(PFEC_INT num_entries,
                                std::vector<std::string> &names);
        void move_to_group_info(PFEC_INT partID,
                                std::string groupName,
                                std::map<std::string, PFEC_INT> &nameToIndexMap,
                                std::vector< std::vector<std::streamoff> > &offsets);
    };
};