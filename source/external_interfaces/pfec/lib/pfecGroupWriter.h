#pragma once

#include "pfecBaseWriter.h"

namespace PFECNativeInterface {

    class GroupWriter : public BaseWriter
    {
    public:

        GroupWriter();
        /// TODO
        void write_header_info(std::vector<PFEC_INT> &localPartIDs,
                               std::vector<std::string> &nodeGroupNames,
                               std::vector<std::string> &elemGroupNames);
        /// TODO
        void write_part_node_group(const PFEC_INT partID,
                                   const std::string nodeGroupName,
                                   const std::vector<PFEC_IND> &indices,
                                   const bool localNumbering = true);
        /// TODO
        void write_part_elem_group(const PFEC_INT partID,
                                   const std::string elemGroupName,
                                   const std::vector<PFEC_IND> &indices,
                                   const bool localNumbering = true);
        /// Updates the lookup table, closes the stream, and resets all variables to the default values.
        virtual void close();

    protected:
        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

        bool header_written_;
        /// The names of the node groups
        std::vector<std::string> node_group_names_;
        /// The names of the node groups
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
        /// Position of node group lookup table
        std::streampos node_group_pos_;
        /// Position of element group lookup table
        std::streampos elem_group_pos_;

        /// Writes and lookup table to the current position in the output stream
        /// given the offsets ordered by group id then by partition in the same
        /// order as the vectors
        void write_lookup_table(std::vector< std::vector<std::streamoff> > &offsets);
        /// Writes the names of a list of groups
        void write_names(std::vector<std::string> &names);
        /// Updates both lookup tables
        void update_lookup_tables();
    };
};