#pragma once

#include "pfecBaseReader.h"
#include <array>

namespace PFECNativeInterface {

    class DataReader : public BaseReader
    {
    public:

        /// Default constructor just sets the some variables to initial values
        DataReader();
        /// Get what type of data is contained within this file
        /// Orientation data is an exception, and the dataNumComp will be 0
        void get_data_info(std::vector<PFEC_INT> &dataIDs,
                           std::vector<DataType> &dataTypes,
                           std::vector<PFEC_INT> &dataNumComp,
                           std::vector<std::string> &dataNames,
                           std::vector< std::vector<std::string> > &dataCompNames) const;
		/// An alternative get_data_info method that returns a vector of DataMetaData objects
		/// instead of filling a bunch of vectors.
		std::vector<DataMetaData> get_data_info() const;
        /// Get the data for the case of a unique value at each node/element.  This method should be used for
        /// data types of Nodal and Element, see DataTypes in Defines.h.  The dataID will be checked against
        /// the header information to make sure the data type is Nodal or Element.
        /// Note: an empty vector will be returned if the number of components is 0 for the dataID, which
        /// implies it is unstructured data and is intended for restarts or something. (use get_part_data_stream_
        /// in this case).
        /// Parameters:
        ///   partID - id of the partition to read from
        ///   dataID - id of the data to read
        ///   data - the vector that will be filled with the data. (Note: it will be emptied if it contains anything).
        ///          It is a vector with each index corresponding to the node/element in the partition mesh with
        ///          the same order.  The length of the inner vector matches the number of components for the dataID.
        void read_part_unique_data(PFEC_INT partID,
                                   PFEC_INT dataID,
                                   std::vector< std::vector<double> > &data);
        /// Get the data for the case of each node/quadpoint has a value by element. This method should be used
        /// for data types of NodalByElement and QuadPoint, see DataTypes in Defines.h.  The dataID will be checked
        /// against the header information to make sure the data type is NodalByElement or QuadPoint.
        /// Note: an empty vector will be returned if the number of components is 0 for the dataID, which
        /// implies it is unstructured data and is intended for restarts or something. (use get_part_data_stream_
        /// in this case).
        /// Parameters:
        ///   partID - id of the partition to read from
        ///   dataID - id fo the data to read
        ///   data - the vector that will be filled with the data. (Note: it will be emptied if it contains anything).
        ///          It is a vector with each index corresponding to the element in the partition mesh with
        ///          the same order.  The next inner vector corresponds to the node/quadpoint in the element in the order
        ///          defined by the element connectivity (if nodal) or the quad point order in program that wrote the
        ///          information.  The length of the inner most vector matches the number of components for the dataID.
        void read_part_by_elem_data(PFEC_INT partID,
                                    PFEC_INT dataID,
                                    std::vector< std::vector< std::vector<double> > > &data);
        ///
        void read_part_orientation_data(PFEC_INT partID,
                                        PFEC_INT dataID,
										std::vector<std::vector<std::array<double, 4>>>& quaterions);
        // TODO implement read_part_by_elem_data using iterators
        /// TODO: Still needs to be implemented
        /// Returns a pointer to the beginning of the data for a particular partition ID and data ID
        /// allowing the programing using this library to process custom formatted data (like quad point data)
        std::ifstream* get_part_data_stream(PFEC_INT partID, PFEC_INT dataID);
        /// Will clean up all the memory allocated and close the mesh file.
        void close();

    protected:
        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

        // Header information
        double time;
        // Information of what data is contained in the data file
        PFEC_INT num_data_blocks_;
        std::vector<PFEC_INT> data_ids_;
        std::map<PFEC_INT, PFEC_INT> data_id_map_;
        std::vector<DataType> data_types_;
        std::vector<PFEC_INT> data_num_comp_;
        std::vector<std::string> data_names_;
        std::vector< std::vector<std::string> > data_comp_names_;
        /// Offsets (relative to end of header) to the beginning of a partition,
        /// which in itself contains another lookup table to get to the data.
        std::vector<std::streamoff> local_part_offsets_;

        /// Reads the header for a particular partition that indicates where data is stored
        std::map<PFEC_INT, std::streamoff> read_part_header(PFEC_INT partID);
        /// Moves the stream to a particular data block within a partition in the file
        void move_to_data_block(PFEC_INT partID, PFEC_INT dataID);
        /// Checks some common things for the methods in the class
        void common_error_checks(PFEC_INT partID, PFEC_INT dataID);

        ///Reads the header information (called by open)
        virtual void read_header_info();
    };

};