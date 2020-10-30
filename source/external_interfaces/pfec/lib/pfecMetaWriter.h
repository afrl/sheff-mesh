#pragma once

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "pfecMetaReader.h"
#include "pfecMeshWriter.h"
#include "pfecDataWriter.h"
#include "pfecGroupWriter.h"

namespace PFECNativeInterface {

    /// TODO
    class MetaWriter : public MetaReader
    {
    public:
		MetaWriter();
        MetaWriter(std::string filePath, std::string analysisName);
        /// Destructor writes any final data and closes the stream
        ~MetaWriter();
		void open(std::string filePath, std::string analysisName);
		/// Set the analysis name of the meta file.  The name can be set at any point.
		void set_analysis_name(const std::string analysisName) { analysis_name_ = analysisName; };
        /// Starts a new time step for the writer.  All writers created via one of the get_X_writer
        /// methods will be associated with this time step.  The time step is closed when a new one
        /// is created or the destructor gets called due to the MetaWriter getting cleaned up.
        void start_new_time_step(int timeStepID, double timeVal);
        /// Returns a unique pointer to a mesh writer object for the specified time step
        std::unique_ptr<PFECNativeInterface::MeshWriter> get_mesh_writer(std::string fileRelPath);
        /// Returns a unique pointer to a data writer that should be used for writing orientation information for the specifed time step
        /// Though it is not enforced, only the orientation data will be read most likely if something else is included
        std::unique_ptr<PFECNativeInterface::DataWriter> get_orient_writer(std::string fileRelPath);
        /// Returns a unique pointer to a data writer object for the specified time step
        std::unique_ptr<PFECNativeInterface::DataWriter> get_data_writer(std::string fileRelPath);
        /// Returns a unique pointer to a group writer object for the specified time step
        std::unique_ptr<PFECNativeInterface::GroupWriter> get_group_writer(std::string fileRelPath);

		/// Set the path of the first time step mesh file
		void set_initial_mesh_path(const std::string meshPath) { time_step_metas_.front().mesh_file_path_ = meshPath; }
		void set_initial_orient_path(const std::string orientPath) { time_step_metas_.front().orient_file_path_ = orientPath; }
		/// A helper method for prepending a directory to the file paths in the TimeStepMeta obejcts
		void prepend_dir_to_paths(const std::string dir);
		/// Clears all existing data files from the TimeStepMeta's
		/// This is used when creating a writer based on the template of another meta file, but you want
		/// to throw away the data before writing new data.
		void clear_existing_data_references();
    private:
		/// Indicates if a timestep was being written
		bool time_step_started_;
		/// Stores the index of the current time step meta that is being written to
		PFEC_INT curr_time_meta_index_;

	private:
        /// Writes the JSON formatted meta file to meta_filepath_
        void write_json_meta();

		TimeStepMeta& get_time_step_meta(const int metaTimeIndex);
    };

};