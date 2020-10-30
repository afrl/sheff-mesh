#include "pfecDataWriter.h"
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <algorithm>
#include <iostream>
#include <chrono>
#include <thread>

using namespace PFECNativeInterface;

#ifdef USE_MPI
// A helper for MPI calls
/// All the MPI processes with this file open will sync to be at the same
/// position in the file
void sync_stream_position(MPI_Comm comm, int procWithData, std::fstream& stream, std::streamoff pos)
{
	unsigned long long posInt = pos;
	int mpiResult = MPI_Bcast(
		&posInt,
		1,
		MPI_UNSIGNED_LONG_LONG,
		procWithData,
		comm);
	if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT(ERR_INTERFACE);
	pos = posInt;
	stream.seekp(pos);
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
}
#endif

PFECNativeInterface::DataWriter::DataWriter(std::string filename) : BaseWriter()
{
	filename_ = filename;
	mode_ = OutputMode::BufferStreams;
    time_ = 0.0;
    num_data_blocks_ = 0;
    part_lookup_pos_ = 0;
}

PFECNativeInterface::DataWriter::~DataWriter()
{
	if (mode_ == OutputMode::DirectToFile) {
		update_lookup_tables();
		cleanup();
	}
}

std::vector< std::vector<std::stringstream*> > DataWriter::initialize_streams(const std::vector<PFEC_INT> &localPartIDs, const std::vector<DataMetaData> &dataMetaData)
{
    initialize(localPartIDs, dataMetaData);

    // Now allocate the string streams
    data_buffers_.resize(num_local_parts_);
    for (PFEC_INT i = 0; i < num_local_parts_; i++) {
        data_buffers_[i].resize(num_data_blocks_);
        for (PFEC_INT j = 0; j < num_data_blocks_; j++)
            data_buffers_[i][j] = new std::stringstream();
    }

    return data_buffers_;
}

void PFECNativeInterface::DataWriter::initialize_for_serial_writing(const std::vector<PFEC_INT>& localPartIDs, const std::vector<DataMetaData>& dataMetaData)
{
	initialize(localPartIDs, dataMetaData);

	BaseWriter::open(filename_);

	write_header();
	write_part_header(local_part_ids_[curr_serial_part_index_]);

	mode_ = OutputMode::DirectToFile;
}

void PFECNativeInterface::DataWriter::start_serial_data_next_partition()
{
	if (curr_serial_part_index_ < local_part_ids_.size() - 1) {
		curr_serial_part_index_++;
		write_part_header(local_part_ids_[curr_serial_part_index_]);
		return;
	}
	throw _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
}

void PFECNativeInterface::DataWriter::write_part_unique_data(PFEC_INT partID, PFEC_INT dataID, std::vector<std::vector<double>>& data)
{
    common_error_checks(partID, dataID);
    PFEC_INT partInd = part_id_map_[partID];
    PFEC_INT dataInd = data_id_map_[dataID];
    unsigned numComps = data_metadata_[data_id_map_[dataID]].num_comps;

	std::ostream* os_ptr;
	if (mode_ == OutputMode::BufferStreams) os_ptr = data_buffers_[partInd][dataInd];
	else if (mode_ == OutputMode::DirectToFile) {
		if (curr_serial_part_index_ != partInd) throw _MSG_EXCEPT_C(ERR_PART_MISMATCH);
		// Save position of the start of the data for calculating offets
		data_pos_[partInd][dataInd] = out_stream_.tellp();
		os_ptr = &out_stream_;
		// Write data id for error checking
		os_ptr->write((char*)&(dataID), sizeof(PFEC_INT));
	}
	else throw _MSG_EXCEPT_C(ERR_UNKNOWN);
    PFEC_INT len = static_cast<PFEC_INT>(data.size());
	os_ptr->write((char*)&len, sizeof(len));
    for (PFEC_INT i = 0; i < len; i++) {
        if (data[i].size() != numComps) throw _MSG_EXCEPT_C(ERR_PARAMETER_WRONG_SIZE);
		os_ptr->write((char*)&(data[i][0]), numComps*sizeof(double));
    }
}

void PFECNativeInterface::DataWriter::write_part_unique_data(PFEC_INT partID, PFEC_INT dataID, PFEC_INT numVecValues, double * data)
{
    common_error_checks(partID, dataID);
    PFEC_INT partInd = part_id_map_[partID];
    PFEC_INT dataInd = data_id_map_[dataID];
    unsigned numComps = data_metadata_[data_id_map_[dataID]].num_comps;

	std::ostream* os_ptr;
	if (mode_ == OutputMode::BufferStreams) os_ptr = data_buffers_[partInd][dataInd];
	else if (mode_ == OutputMode::DirectToFile) {
		if (curr_serial_part_index_ != partInd) throw _MSG_EXCEPT_C(ERR_PART_MISMATCH);
		// Save position of the start of the data for calculating offets
		data_pos_[partInd][dataInd] = out_stream_.tellp();
		os_ptr = &out_stream_;
		// Write data id for error checking
		os_ptr->write((char*)&(dataID), sizeof(PFEC_INT));
	}
	else throw _MSG_EXCEPT_C(ERR_UNKNOWN);
    os_ptr->write((char*)&numVecValues, sizeof(numVecValues));
    for (PFEC_INT i = 0; i < numVecValues; i++) {
        os_ptr->write((char*)&(data[i*numComps]), numComps*sizeof(double));
    }
}

void PFECNativeInterface::DataWriter::write_part_by_elem_data(PFEC_INT partID, PFEC_INT dataID, std::vector<std::vector<std::vector<double>>>& data)
{
    common_error_checks(partID, dataID);
    PFEC_INT partInd = part_id_map_[partID];
    PFEC_INT dataInd = data_id_map_[dataID];
    unsigned numComps = data_metadata_[data_id_map_[dataID]].num_comps;

	std::ostream* os_ptr;
	if (mode_ == OutputMode::BufferStreams) os_ptr = data_buffers_[partInd][dataInd];
	else if (mode_ == OutputMode::DirectToFile) {
		if (curr_serial_part_index_ != partInd) throw _MSG_EXCEPT_C(ERR_PART_MISMATCH);
		// Save position of the start of the data for calculating offets
		data_pos_[partInd][dataInd] = out_stream_.tellp();
		os_ptr = &out_stream_;
		// Write data id for error checking
		os_ptr->write((char*)&(dataID), sizeof(PFEC_INT));
	}
	else throw _MSG_EXCEPT_C(ERR_UNKNOWN);
    PFEC_INT num_elems = static_cast<PFEC_INT>(data.size());
    os_ptr->write((char*)&num_elems, sizeof(num_elems));
    for (PFEC_INT i = 0; i < num_elems; i++) {
        // Write the number of nodes/quads in this element
        PFEC_INT len = static_cast<PFEC_INT>(data[i].size());
        os_ptr->write((char*)&len, sizeof(len));
        // Write the data for this element
        for (PFEC_INT j = 0; j < len; j++)
            os_ptr->write((char*)&(data[i][j][0]), numComps*sizeof(double));
    }
}

void PFECNativeInterface::DataWriter::write_part_orientation_data(PFEC_INT partID, PFEC_INT dataID, std::vector<std::vector<std::array<double, 4>>>& quaterions)
{
    common_error_checks(partID, dataID);

    PFEC_INT partInd = part_id_map_[partID];
    PFEC_INT dataInd = data_id_map_[dataID];
    unsigned numComps = data_metadata_[data_id_map_[dataID]].num_comps;

	std::ostream* os_ptr;
	if (mode_ == OutputMode::BufferStreams) os_ptr = data_buffers_[partInd][dataInd];
	else if (mode_ == OutputMode::DirectToFile) {
		if (curr_serial_part_index_ != partInd) throw _MSG_EXCEPT_C(ERR_PART_MISMATCH);
		// Save position of the start of the data for calculating offets
		data_pos_[partInd][dataInd] = out_stream_.tellp();
		os_ptr = &out_stream_;
		// Write data id for error checking
		os_ptr->write((char*)&(dataID), sizeof(PFEC_INT));
	}
	else throw _MSG_EXCEPT_C(ERR_UNKNOWN);
    PFEC_INT num_elems = static_cast<PFEC_INT>(quaterions.size());
    os_ptr->write((char*)&num_elems, sizeof(num_elems));
    for (PFEC_INT i = 0; i < num_elems; i++) {
        // Write the number of orientations in this element
        PFEC_INT len = static_cast<PFEC_INT>(quaterions[i].size());
        os_ptr->write((char*)&len, sizeof(len));
        // Write the number of rotations
        PFEC_INT numQuaterionTerms = 4;
        os_ptr->write((char*)&numQuaterionTerms, sizeof(numQuaterionTerms));
        // Write quaterion info
        for (PFEC_INT j = 0; j < len; j++)
            os_ptr->write((char*)&(quaterions[i][j][0]), numQuaterionTerms*sizeof(double));
    }
}

void DataWriter::finalize_streams()
{
    // Calculate the total size of the file and allocate it (via system call?)

    BaseWriter::open(filename_);

    write_header();

    // For each partition, write the partition lookup table, then for each write the data bloacks for that partition
    for (PFEC_INT i = 0; i < num_local_parts_; i++) {
        write_part_header(local_part_ids_[i]);
        for (PFEC_INT j = 0; j < num_data_blocks_; j++) {
            // Save position of the start of the data for calculating offets
            data_pos_[i][j] = out_stream_.tellp();
            // Write data id for error checking
            out_stream_.write((char*)&(data_metadata_[j].id), sizeof(PFEC_INT));
            // Get number of bytes
            data_buffers_[i][j]->seekg(0, std::ios::end);
            std::streamoff size = data_buffers_[i][j]->tellg();
			if (size == std::streamoff(0)) throw _MSG_EXCEPT_C(ERR_PART_WAS_NOT_WRITTEN);
            // Set the get pointer to the beginning
            data_buffers_[i][j]->seekg(0, std::ios::beg);
            // Write via the buffer
            out_stream_ << data_buffers_[i][j]->rdbuf();
            //out_stream_.write((char*)&(data_buffers_[i][j]->str().c_str()), size); // This makes a fully copy in memory
            // See http://stackoverflow.com/questions/28663075/how-to-get-characters-out-of-stringstream-without-copy
            // And http://stackoverflow.com/questions/4064601/copy-data-from-fstream-to-stringstream-with-no-buffer
        }
    }

    update_lookup_tables();

    // Close the stream and clean up
    cleanup();
}

void PFECNativeInterface::DataWriter::finalize_streams(const int mpiCComm, const bool outputStatus)
{
#ifdef USE_MPI
	// Basic flow:
	// Calculate the total size of the file via MPI commands
	// Rank 0 of the group creates the file and writes the headers
	// Note: each group is intended to have one data file output
	// I have a problem here:
	//  Option 1) just loop over the processors and let each write serially within the group
	//  Option 2) figure out how to allocate a null file of the right size, then write in parallel
	//            but limit the amount of processors allowed to write at once to avoid bottlenecking

	int mpiResult;
	// Get the number of processes
	int commSize;
	mpiResult = MPI_Comm_size(mpiCComm, &commSize);
	if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

	// Get the rank of the process
	int commRank;
	mpiResult = MPI_Comm_rank(mpiCComm, &commRank);
	if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

	// Only let the first processor write the header information
	if (commRank == 0) {
		BaseWriter::open(filename_);
		write_header();
	}

	// Let the others open in append mode
	if (commRank != 0)
		BaseWriter::open(filename_, true);

	// Every processor needs to know where the header ends to they can write their own
	// lookup tables
	sync_stream_position(mpiCComm, 0, out_stream_, header_end_);
	if (outputStatus) {
		for (PFEC_INT i = 0; i < num_local_parts_; i++) {
			if (commRank == i) std::cout << "Current position for " << i << ": " << out_stream_.tellp() << std::endl;
			MPI_Barrier(mpiCComm);
		}
	}

	// For each partition, write the partition lookup table, then for each write the data bloacks for that partition
	// Note: this effectively limits it to once process writing at a time
	// TODO: work on how to get pools of processes writing at one time
	for (PFEC_INT i = 0; i < num_local_parts_; i++) {
		// Let the root write the partition header since it will be responsible
		// for updating the lookup tables later
		if (commRank == 0)
			write_part_header(i);

		// Wait until all processors reach this point
		mpiResult = MPI_Barrier(mpiCComm);
		if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

		// Get all processors to be at the header of the data, to get ready for the next processor to write its data
		sync_stream_position(mpiCComm, 0, out_stream_, out_stream_.tellp());

		for (PFEC_INT j = 0; j < num_data_blocks_; j++) {
			// Save position of the start of the data for calculating offets
			data_pos_[i][j] = out_stream_.tellp();

			// Only the processor with the data does the writing
			if (i == commRank) {
				if (outputStatus) std::cout << "Writing part " << commRank << " for data " << data_metadata_[j].id << " at " << out_stream_.tellp() << " for buffer " << &(*data_buffers_[i][j]) << std::endl;
				// Write data id for error checking
				out_stream_.write((char*)&(data_metadata_[j].id), sizeof(PFEC_INT));
				// Get number of bytes
				data_buffers_[i][j]->seekg(0, std::ios::end);
				std::streamoff size = data_buffers_[i][j]->tellg();
				if (size == std::streamoff(0)) {
					std::cout << "Currenlty outputting to " << filename_ << std::endl;
					std::cout << "P" << commRank << ": Partition " << i << " data was 0 for buffer " << j << " at " << &(*data_buffers_[i][j]) << std::endl;
					std::cout << "Raw data: " << std::endl << data_buffers_[i][j]->str() << std::endl;
					throw _MSG_EXCEPT_C(ERR_PART_WAS_NOT_WRITTEN);
				}
				// Set the get pointer to the beginning
				data_buffers_[i][j]->seekg(0, std::ios::beg);
				// Write via the buffer
				out_stream_ << data_buffers_[i][j]->rdbuf();
				//out_stream_.write((char*)&(data_buffers_[i][j]->str().c_str()), size); // This makes a fully copy in memory
				// See http://stackoverflow.com/questions/28663075/how-to-get-characters-out-of-stringstream-without-copy
				// And http://stackoverflow.com/questions/4064601/copy-data-from-fstream-to-stringstream-with-no-buffer
			}

			// Wait until all processors reach this point
			mpiResult = MPI_Barrier(mpiCComm);
			if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

			// Get all processors to be at the end of the data, to get ready for the next
			// processor to write its data
			sync_stream_position(mpiCComm, i, out_stream_, out_stream_.tellp());
			mpiResult = MPI_Barrier(mpiCComm);
			if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);
			/*
			for (PFEC_INT i = 0; i < num_local_parts_; i++) {
				std::string mess = "Current position for " + std::to_string(i) + ": " + std::to_string(out_stream_.tellp()) + "\n";
				if (commRank == i) {
					std::cout << mess;
				}
				mpiResult = MPI_Barrier(mpiCComm);
				if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}
			mpiResult = MPI_Barrier(mpiCComm);
			if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);
			*/
		}
	}

	// Let the root write update the lookup tables
	if (commRank == 0) {
		//for (auto i : parts_pos_)
		//	std::cout << "Pos " << i << std::endl;
		update_lookup_tables();
	}

	// Close the stream and clean up
	cleanup();
#else
	throw _MSG_EXCEPT_C(ERR_NO_MPI);
#endif
}

void PFECNativeInterface::DataWriter::finalize_streams_improved(const int mpiCComm, const bool outputStatus)
{
#ifdef USE_MPI
	// Basic flow:
	// Calculate the total size of the file via MPI commands
	// Rank 0 of the group creates the file and writes the headers
	// Note: each group is intended to have one data file output
	// I have a problem here:
	//  Option 1) just loop over the processors and let each write serially within the group
	//  Option 2) figure out how to allocate a null file of the right size, then write in parallel
	//            but limit the amount of processors allowed to write at once to avoid bottlenecking

	int mpiResult;
	// Get the number of processes
	int commSize;
	mpiResult = MPI_Comm_size(mpiCComm, &commSize);
	if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

	// Get the rank of the process
	int commRank;
	mpiResult = MPI_Comm_rank(mpiCComm, &commRank);
	if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

	// Compute how much space is needed for each block of data on this process
	std::vector<long long> dataBlockSizes;
	for (PFEC_INT j = 0; j < num_data_blocks_; j++) {
		data_buffers_[commRank][j]->seekg(0, std::ios::end);
		dataBlockSizes.push_back(static_cast<long long>(data_buffers_[commRank][j]->tellg()));
	}

	// Let the root process create the file with blank data first
	if (commRank == 0) {
		// Tell the root process how much data this proc needed
		std::vector<long long> dataBlockSizesUnrolled(commSize*num_data_blocks_);
		MPI_Gather(
			&dataBlockSizes[0],
			num_data_blocks_,
			MPI_LONG_LONG,
			&dataBlockSizesUnrolled[0],
			num_data_blocks_,
			MPI_LONG_LONG,
			0,
			mpiCComm);

		// Find the largest data block size and create a buffer to dump
		auto maxSize = *std::max_element(dataBlockSizesUnrolled.begin(), dataBlockSizesUnrolled.end());
		std::vector<char> blankData(maxSize, '0');

		// Write the file with blank data.
		BaseWriter::open(filename_);
		write_header();
		for (PFEC_INT i = 0; i < num_local_parts_; i++) {
			write_part_header(i);
			for (PFEC_INT j = 0; j < num_data_blocks_; j++) {
				// Save position of the start of the data for calculating offets
				data_pos_[i][j] = out_stream_.tellp();
				// Write data id for error checking
				out_stream_.write((char*)&(data_metadata_[j].id), sizeof(PFEC_INT));
				// Write the blank data
				out_stream_.write(&(blankData[0]), dataBlockSizesUnrolled[i*num_data_blocks_ + j]);
			}
		}

		update_lookup_tables();
	}
	else {
		// All other processes (besides rank 0) just sends the data
		MPI_Gather(
			&dataBlockSizes[0],
			num_data_blocks_,
			MPI_LONG_LONG,
			NULL,
			num_data_blocks_,
			MPI_LONG_LONG,
			0,
			mpiCComm);
	}

	// Wait until all processors reach this point
	mpiResult = MPI_Barrier(mpiCComm);
	if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);

	// Let the others open in append mode
	if (commRank != 0)
		BaseWriter::open(filename_, true);
	// Get the beginning of the file
	out_stream_.seekp(0, std::ios::beg);
	std::streampos fileBegin = out_stream_.tellp();

	// Scatter where in the files each process should write each data block
	std::vector<std::streamoff> myDataOffsets(num_data_blocks_);
	if (commRank == 0) {
		// Build an unrolled vector of the data positions in the file
		std::vector<std::streamoff> dataOffsetsUnrolled;
		for (PFEC_INT i = 0; i < num_local_parts_; i++)
			for (PFEC_INT j = 0; j < num_data_blocks_; j++)
				dataOffsetsUnrolled.push_back(data_pos_[i][j]);

		MPI_Scatter(
			reinterpret_cast<char*>(&dataOffsetsUnrolled[0]),
			num_data_blocks_ * sizeof(std::streamoff),
			MPI_CHAR,
			reinterpret_cast<char*>(&myDataOffsets[0]),
			num_data_blocks_ * sizeof(std::streamoff),
			MPI_CHAR,
			0,
			mpiCComm);
	}
	else {
		MPI_Scatter(
			NULL,
			num_data_blocks_ * sizeof(std::streamoff),
			MPI_CHAR,
			reinterpret_cast<char*>(&myDataOffsets[0]),
			num_data_blocks_ * sizeof(std::streamoff),
			MPI_CHAR,
			0,
			mpiCComm);
	}

	// Update the data positions using the offset information that was communicated via MPI
	for (PFEC_INT j = 0; j < num_data_blocks_; j++) {
		data_pos_[commRank][j] = fileBegin + myDataOffsets[j];
		//out_stream_.seekp(data_pos_[commRank][j]);
		//std::cout << out_stream_.tellp() << std::endl;
	}

	// Now let each process write the data
	for (PFEC_INT i = 0; i < num_local_parts_; i++)
	{
		// Only the processor with the data does the writing
		if (i == commRank)
		{
			for (PFEC_INT j = 0; j < num_data_blocks_; j++)
			{
				// Move the appropriate position in the file
				out_stream_.seekp(data_pos_[i][j]);
				if (outputStatus) std::cout << "Writing part " << commRank << " for data " << data_metadata_[j].id << " at " << out_stream_.tellp() << " for buffer " << &(*data_buffers_[i][j]) << std::endl;
				// Write data id for error checking
				out_stream_.write((char*)&(data_metadata_[j].id), sizeof(PFEC_INT));
				//std::cout << out_stream_.tellp() << std::endl;
				// Write via the buffer
				// Set the get pointer to the beginning
				data_buffers_[i][j]->seekg(0, std::ios::beg);
				std::copy(std::istreambuf_iterator<char>(*data_buffers_[i][j]),
					std::istreambuf_iterator<char>(),
					std::ostreambuf_iterator<char>(out_stream_));
				// See http://stackoverflow.com/questions/28663075/how-to-get-characters-out-of-stringstream-without-copy
				// And http://stackoverflow.com/questions/4064601/copy-data-from-fstream-to-stringstream-with-no-buffer
				//std::cout << out_stream_.tellp() << std::endl;
			}
		}

		// Wait until all processors reach this point
		mpiResult = MPI_Barrier(mpiCComm);
		if (mpiResult != MPI_SUCCESS) throw _MSG_EXCEPT_C(ERR_INTERFACE);
	}

	// Close the stream and clean up
	cleanup();
#else
	throw _MSG_EXCEPT_C(ERR_NO_MPI);
#endif
}

void PFECNativeInterface::DataWriter::initialize(const std::vector<PFEC_INT> &localPartIDs, const std::vector<DataMetaData> &dataMetaData)
{
    // Save the parameters needed later for writing the header and processing the information
    num_local_parts_ = static_cast<PFEC_INT>(localPartIDs.size());
    local_part_ids_ = localPartIDs;
	for (PFEC_INT i = 0; i < num_local_parts_; i++)
		part_id_map_[local_part_ids_[i]] = i;
    num_data_blocks_ = static_cast<PFEC_INT>(dataMetaData.size());
    data_metadata_ = dataMetaData;
    for (PFEC_INT i = 0; i < num_data_blocks_; i++)
        data_id_map_[dataMetaData[i].id] = i;

    // Resize position vectors
    parts_pos_.resize(num_local_parts_, 0);
    parts_lookup_start_pos_.resize(num_local_parts_, 0);
    parts_lookup_end_pos_.resize(num_local_parts_, 0);
    data_pos_.resize(num_local_parts_, std::vector<std::streampos>(num_data_blocks_, 0));
}

void PFECNativeInterface::DataWriter::write_header()
{
    write_basic_header();

    // Write header specific to data files
    out_stream_.write((char*)&time, sizeof(double));
	out_stream_.write((char*)&num_local_parts_, sizeof(PFEC_INT));
    out_stream_.write((char*)&num_data_blocks_, sizeof(PFEC_INT));
    for (PFEC_INT i = 0; i < num_data_blocks_; i++) {
        out_stream_.write((char*)&(data_metadata_[i].id), sizeof(PFEC_INT));
        out_stream_.write((char*)&(data_metadata_[i].type), sizeof(DataType));
        out_stream_.write((char*)&(data_metadata_[i].num_comps), sizeof(PFEC_INT));
        write_string(data_metadata_[i].name);
        if (data_metadata_[i].num_comps > 0)
            for (PFEC_INT j = 0; j < data_metadata_[i].num_comps; j++)
                write_string(data_metadata_[i].comp_names[j]);
    }
    // Write an empty partition lookup table that will be updated later.
    part_lookup_pos_ = out_stream_.tellp();
	std::vector<std::streamoff> emptyOffsets;
    write_part_lookup_table(local_part_ids_, emptyOffsets);

    header_end_ = out_stream_.tellp();
	//std::cout << "Wrote header to " << header_end_ << std::endl;
}

void PFECNativeInterface::DataWriter::write_part_header(PFEC_INT partID)
{
    // Save position of the beginning of the partition data for calculating offsets
    parts_pos_[part_id_map_[partID]] = out_stream_.tellp();
	//std::cout << "Writing part header " << partID << " at " << parts_pos_[part_id_map_[partID]] << std::endl;
    // Write partition ID for error checking
    out_stream_.write((char*)&partID, sizeof(partID));
    // Save position of the lookup table to update later
    parts_lookup_start_pos_[part_id_map_[partID]] = out_stream_.tellp();
    // Write an empty lookup table
	std::vector<std::streamoff> emptyOffsets;
    write_part_data_lookup_table(emptyOffsets);
    // Save position of the end of the lookup table to calcualte offsets later
    parts_lookup_end_pos_[part_id_map_[partID]] = out_stream_.tellp();
}

void PFECNativeInterface::DataWriter::write_part_data_lookup_table(std::vector<std::streamoff>& offsets)
{
    std::streamoff defaultOff = -1;
    for (PFEC_INT i = 0; i < num_data_blocks_; i++) {
        out_stream_.write((char*)&(data_metadata_[i].id), sizeof(PFEC_INT));
        if (offsets.size() != num_data_blocks_)
            out_stream_.write((char*)&defaultOff, sizeof(std::streamoff));
        else
            out_stream_.write((char*)&offsets[i], sizeof(std::streamoff));
    }
}

void PFECNativeInterface::DataWriter::update_lookup_tables()
{
    // Move to the partition lookup table and overwrite the table with the
    // offets relative to the end of the main header.
	//std::cout << "Part Lookup pos: " << part_lookup_pos_ << std::endl;
    out_stream_.seekp(part_lookup_pos_);
    std::vector<std::streamoff> partsOff(parts_pos_.size());
	for (PFEC_INT i = 0; i < partsOff.size(); i++) {
		if (parts_pos_[i] == std::streampos(0)) throw _MSG_EXCEPT_C(ERR_PART_WAS_NOT_WRITTEN);
		partsOff[i] = parts_pos_[i] - header_end_;
		//std::cout << "Partition " << i << " offset: " << partsOff[i] << std::endl;
	}
    write_part_lookup_table(local_part_ids_, partsOff);

    // Move to the position of each lookup table within the partition and
    // overwrite it with the correct offsets
    for (PFEC_INT i = 0; i < num_local_parts_; i++) {
		//std::cout << "Part lookup start pos: " << parts_lookup_start_pos_[i] << std::endl;
        out_stream_.seekp(parts_lookup_start_pos_[i]);
        std::vector<std::streamoff> dataOffs(num_data_blocks_);
        for (PFEC_INT j = 0; j < num_data_blocks_; j++)
            dataOffs[j] = data_pos_[i][j] != std::streampos(0) ? data_pos_[i][j] - parts_lookup_end_pos_[i] : -1;
        write_part_data_lookup_table(dataOffs);
    }
}

void PFECNativeInterface::DataWriter::common_error_checks(PFEC_INT partID, PFEC_INT dataID)
{
    // Check part is within bounds
    if (part_id_map_.count(partID) == 0)
        throw _MSG_EXCEPT_C(ERR_PART_ID_OUT_OF_BOUND);
    // Check data ID is within bounds
    if (data_id_map_.count(dataID) == 0)
        throw _MSG_EXCEPT_C(ERR_DATA_ID_OUT_OF_BOUND);
}

void PFECNativeInterface::DataWriter::cleanup()
{
	filename_ = "";
    time_ = 0.0;
    num_data_blocks_ = 0;
    data_metadata_.clear();
    // Clean up buffers
    for (PFEC_INT i = 0; i < data_buffers_.size(); i++) {
        for (PFEC_INT j = 0; j < data_buffers_[i].size(); j++)
            delete data_buffers_[i][j];
        data_buffers_[i].clear();
    }
    data_buffers_.clear();
    part_lookup_pos_ = 0;
    parts_pos_.clear();
    parts_lookup_start_pos_.clear();
    parts_lookup_end_pos_.clear();
    data_pos_.clear();

    BaseWriter::close();
}

