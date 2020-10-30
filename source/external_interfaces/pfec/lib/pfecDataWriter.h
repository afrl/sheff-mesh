#pragma once

#include "pfecBaseWriter.h"
#include <sstream>
#include <array>

namespace PFECNativeInterface {

    /// This class works fine, but it needs to be reworked to make it more intutive.
    /// A class used to write data files.  It functions a little different than the
    /// other writers.  There are primarily two ways to get data into the data file.
    /// Option 1: call initialize_streams, write the data to the streams directly
    ///           (but note that the calling program will be responsible to format
    ///           the data in the correct format), and call finalize_streams (which
    ///           is when the actual file is created and written to.  It should be
    ///           noted that for the case of multiple MPI processes trying to write
    ///           to the same file, this is the method that will eventually be evolved
    ///           to support MPI writes.
    /// Option 2: call initialize_streams but ignore the streams; use write_part_unique_data,
    ///           write_part_by_elem_data, and write_part_orientation_data methods;
    ///           call finalize after the writing is done.  This is the safer method
    ///           if you are worried about adhering to the file format.  The order
    ///           of which partition is not important, but all the data for each partition
    ///           must be written at once. Only one process can write to the file using this method.
    /// TODO: add examples for each option
    class DataWriter : public BaseWriter
    {
    public:

        /// Default constructor
        DataWriter(std::string filename);
		virtual ~DataWriter();

        /// Initializes the streams that Beta (or whatever) will write to.  The streams will stay in memory (RAM) until
        /// a finalize_streams call is made.  Whether architecture is distributed or centralized, this method is called.
        /// It is the reponsibility of the program using this library to ensure that each partition has the same data
        /// output, since technically, this library will write it without issue but the reader will have issues.
        /// Parameters:
        /// \param localPartIDs - the IDs of the partitions the parent process is responsible for writing. Note: in
        ///                  the case of a distributed program, each MPI process is responsible for one partition
        ///                  in general.  If some weird hybrid architecture is used, then perhaps the MPI process
        ///                  is responsible for multiple partitions.  These partitions should correspond directly
        ///                  to partitions used in the mesh.  The writer is loosely associated with the mesh and
        ///                  does not have access to the data needed to varify the partitions line up, so the
        ///                  program using this library must ensure the partitions are the same.
        /// \param dataMetaData contains all the neccessary metadata including the id, type, number of components
        ///  names, and names of components.  See DataMetaData.
        /// Returns a vector of size of localPartIndices with each index corresponding to the respective localPartIndices.
        /// Each partition has a vector of pointers to stringstreams (one for each data block to write).  For example, if
        /// the process is responsible for one partition and two data blocks (maybe displacements and stresses) then the
        /// return object will be a vector of length 1, which contians a vector of length 2.  The order of streams
        /// will correspond to the order of the dataIds, dataTypes, dataNumComp, and dataNames passed (which all must be
        /// consistent with each other).
        std::vector< std::vector<std::stringstream*> > initialize_streams(const std::vector<PFEC_INT> &localPartIDs,
                                                                          const std::vector<DataMetaData> &dataMetaData);

		/// TODO
		void initialize_for_serial_writing(const std::vector<PFEC_INT> &localPartIDs,
										   const std::vector<DataMetaData> &dataMetaData);
		///
		void start_serial_data_next_partition();

        // Methods used to write via passing to methods (option 2)

        /// TODO: description
        /// Parameters:
        ///   partID - id of the partition to read from
        ///   dataID - id fo the data to read
        ///   data - the vector that will be filled with the data. (Note: it will be emptied if it contains anything).
        ///          It is a vector with each index corresponding to the node/element in the partition mesh with
        ///          the same order.  The length of the inner vector matches the number of components for the dataID.
        void write_part_unique_data(PFEC_INT partID,
                                    PFEC_INT dataID,
                                    std::vector< std::vector<double> > &data);
        /// TODO: description
        /// This is basically the low level point version of write_part_unique_data
        /// \param data is organized by node then by component (assumes the same number of comps per node)
        void write_part_unique_data(PFEC_INT partID,
                                    PFEC_INT dataID,
                                    PFEC_INT numVecValues,
                                    double * data);
        /// TODO: description
        /// Parameters:
        ///   partID - id of the partition to read from
        ///   dataID - id fo the data to read
        ///   data - the vector that will be filled with the data. (Note: it will be emptied if it contains anything).
        ///          It is a vector with each index corresponding to the element in the partition mesh with
        ///          the same order.  The next inner vector corresponds to the node/quadpoint in the element in the order
        ///          defined by the element connectivity (if nodal) or the quad point order in program that wrote the
        ///          information.  The length of the inner most vector matches the number of components for the dataID.
        void write_part_by_elem_data(PFEC_INT partID,
                                     PFEC_INT dataID,
                                     std::vector< std::vector< std::vector<double> > > &data);
        /// TODO: description
        void write_part_orientation_data(PFEC_INT partID,
										 PFEC_INT dataID,
                                         std::vector< std::vector< std::array<double, 4> > > &quaterions);

        /// Used for a centralized analysis.  Finalizes all streams, see initialize_streams, and writes the data
        /// to the file.  The filename can be relative to the current working directory, or an absolute path can
        /// be given.
        /// Parameters:
        ///   filename - the name of the output file relative to the current working directory or an absolute file path.
        ///              Note: an extension is not assumed, so one must be provided.
        /// TODO: should an extension be assumed?
        void finalize_streams();
        /// Used for a distributed analysis.  Finalizes all streams, see initialize_streams, communicates the length
        /// of the data in the buffer streams, the rank 0 process will write the header information, then each process
        /// or group of processes will write their data to the file until every process is done.  This will be blocking,
        /// so the return will occur when all the processes are done.
        /// Parameters:
        ///   filename - the name of the output file relative to the current working directory or an absolute file path.
        ///              Note: an extension is not assumed, so one must be provided.
        ///   mpiRank - the rank fo this process
        ///   mpiCComm - the comm pool used to communicate with the other processors trying to output.
        ///   numProcsToWriteConcurrently - the number of processors allowed to write at the same time.  If -1, then
        ///                                 there is no limit, but this can have negative performance effects if there
        ///                                 is a lot more processors than physical hard drive headers writing.  To fix
        ///                                 this, limit the number of processors that can write concurrently to avoid
        ///                                 overloading the headers.
        void finalize_streams(const int mpiCComm, const bool outputStatus = false);
		void finalize_streams_improved(const int mpiCComm, const bool outputStatus = false);

		//////////////////
		// Misc helpers //
		//////////////////

		/// A super simple helper to get the DataMetaData struct to be used with orientation data.
		DataMetaData get_orientation_meta() const {
			DataMetaData orientMeta(0, PFECNativeInterface::DataType::Orientation, 4, "Quaterions", { "b0", "b1", "b2", "b3" });
			return orientMeta;
		};

		/// Not used, so it will through a not implemented error
		virtual void open(std::string filePath) {
			throw _MSG_EXCEPT_C(ERR_NOT_IMPLEMENTED);
		};
		/// Not used, so it will through a not implemented error
		virtual void close() {
			throw _MSG_EXCEPT_C(ERR_NOT_IMPLEMENTED);
		};

    protected:
		enum class OutputMode : char {
			BufferStreams,
			DirectToFile
		};

        //////////////////////////////////
        // Protected properties/methods //
        //////////////////////////////////

		std::string filename_;
		OutputMode mode_;
		/// This is only used if the mode is set to serial writing
		PFEC_INT curr_serial_part_index_;
        double time_;
        PFEC_INT num_data_blocks_;
        std::vector<DataMetaData> data_metadata_;
        std::map<PFEC_INT, PFEC_INT> data_id_map_;
        // The first index corresponds to the local index of the partition (0-num_local_parts_)
        // The second index corresponds to the local index of the data block
        std::vector< std::vector<std::stringstream*> > data_buffers_;
        /// Position of the partition lookup table
        std::streampos part_lookup_pos_;
        /// Stores the positions of the start of the partition data
        std::vector<std::streampos> parts_pos_;
        /// Stores the positions of the start of the lookup table for each parition
        std::vector<std::streampos> parts_lookup_start_pos_;
        /// Stores the positions of the end of the lookup table for each parition
        std::vector<std::streampos> parts_lookup_end_pos_;
        /// Stores the positions of the start of data blocks as they are written
        /// First index: partition ordered by local_part_ids_
        /// Second index: data ordered by data_ids_
        std::vector< std::vector<std::streampos> > data_pos_;

        // Helper methods

        /// Initializes alot of variables and vectors
        void initialize(const std::vector<PFEC_INT> &localPartIDs, const std::vector<DataMetaData> &dataMetaData);
        /// Writes the the header information that goes at the beginning
        /// of the data file
        void write_header();
        /// Writes a partition header, returns position of the end of the table
        void write_part_header(PFEC_INT partID);
        /// Writes a lookup table that goes in the partition headers.  The offsets in the
        /// tbale are relative to the end of the partition header.
        void write_part_data_lookup_table(std::vector<std::streamoff> &offsets);
        /// Updates all header information
        void update_lookup_tables();
        /// Checks some common things for the methods in the class
        void common_error_checks(PFEC_INT partID, PFEC_INT dataID);
        /// Resets all variables to default values and cleans up any pointers. Note: this is called by finalize_streams,
        /// sat which point the stringstream points will be invalid.
        void cleanup();
    };

};